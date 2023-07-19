

// c
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


// c++ 
#include <cmath> // floor function
#include <vector>
#include <iostream>
#include <fstream>
#include <istream>
#include <sstream>
#include <string>
#include <iomanip> 
#include <ctime> //<MT apr19> for struct tm (timeformat => calc summertime)
#include <algorithm>  // for using std::sort

using namespace std;

// own
#define PI acos(-1.)
#define SQR(x) ((x)*(x))

#include "Statistics.h"
#include "RandomUtils.h" // contains, e.g.,  myRand()
#include "InOut.h"


// constants

static const int MAXSTR=500;
static const double lat0=37.9804;  //latitude of last drone1 veh at first time
static const double lon0=23.7353;  //longitude 
static const double radiusEarth=6371000.8;
static const double DT=0.04; // sampling time (corresp. to 25fps)

// factors by comparing the simple cos(lat) formula with a
// coordinate system calculator for lat0 and long0, dlon=0.01 (dec degrees)
// with calculated dx/dy using fact_lon=cos
// (ellipsoids introduce much smaller error => using spherical Earth)
// notice that this calculator also has a latitude->Northern dependence
// on the longitude: Artifact of planification by Gauss-Krueger => ignore

static const double fact_lon=878.5/876.46*cos(lat0*PI/180.);


// rotate xy coords to make the left-hand or right-hand part of the
// investigation area parallel to the (rectangular) grid of streets

static const double rotAngle1=0.999; // right lower part (e.g., drone d1)
static const double rotAngle2=0; // left upper part (e.g., drone d9)



//#####################################################
// data structures
//#####################################################

// for the road axis

struct RoadData{
  double rotAngle;  // all the following with respect to rotAngle
  vector<double> x;
  vector<double> y;
  vector<double> coshead; //!! not just heading! smooth and intp will fail!
  vector<double> sinhead;
};


// a single vehicle trajectory in the time (not arclength) domain
// e.g., for the logical trajectory data (along/perp to road axis)

struct Trajectory {
  int ID;
  int type;
  vector<double> time;
  vector<double> x;
  vector<double> y;
  vector<double> heading;
  vector<double> vx;
  vector<double> vy;
  vector<double> ax;
  vector<double> ay;
};

// a vehicle state used in the logical coordinates (along/perp to road axis)
// a vector<vehicle> represents a snapshot of all vehicles
// (which be sorted with respect to one of the vehicle's components)
// while a trajectory represents a single vehicle at all times

struct vehicle {
  int ID;
  int type;
  double time;
  double x;
  double y;
  int lane;
  double heading;
  double vx;
  double vy;
  double ax;
  double ay;
  
  // index of the original vector of trajectories
  // for debugging to check if the time snapshot is correct
  int iveh; 
};


// compare function for sorting a vector of vehicles according to vehicle.x

bool compare_by_x(const vehicle &a, const vehicle &b){
  return a.x < b.x;
}

// compare function for sorting a vector of vectors
// according to element [1] of subvector
// (vector<vector<double>> lanesFinal;)

bool compare_by_secondEntry(const vector<double> &a, const vector<double> &b){
  return a[1] < b[1];
}

// vehicle properties

struct VehicleProperties{
  string type_str;
  int type_ID;
  double length;
  double width;
};
    


  
// collection of global parameters

struct TuningParameters{
  
  // WhatToDo=0:  initial global parameters

  double rotAngle;
  bool  filterForCars; // filters for cars (best lane-keeping behaviour)
  bool  filterForOrigin; // filters for origin location
  double origin_ymin; 
  double origin_ymax; 
  bool  filterForDest;
  double dest_ymin;
  double dest_ymax;
  bool filterForXY;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  int dit;
  
  // WhatToDo=1:  for representation data_s and corresp output

  double ds;
  double dsSmooth_heading; //[m] for data_s and corresp output 

  // WhatToDo=2:  heatmap and corresp output
  
  double gridsize;

  // WhatToDo=3:  for purely data-based road and lane finding/matching
  
  double distBetweenCuts;   //[m], x distance for identifying lanes
  double dsSmooth_grid;    //[m] for finding the y maxima of the heat
  double maxHeadMismatch; //[rad] cutoff for identifying contiguous lanes
  double dsSmooth_lane;     //[m] Smoothing of final lane y coords
  double wLane_min;       // [m] averaging for robust finding of heatpeak
  double maxSinHead;  // heatmap only for sin(heading)<maxSinHead (clear lanes)
  int ncutmin;        // only create a contig. lane if >ncutmin heatpeaks
  
  // WhatToDo=4: purely data-based road and lane finding/matching

  double wLane;    // lanewidth[m] 

  // WhatToDo=5: find neighbors and direct CF leaders

  double dxmax;   // max long offset for qualifying as leader/follower
  double dymax;   //max lateral offset for qualifying as neighbor
  double TL_dxmin;  // consider TL as virtual leader if dx<TL_dxmin
  double TL_bsafe;  // ... or if bkin>TL_bsafe
  double TL_bcrit; // Just pass the yellow TL if bkin>TL_bcit

};


void init(int WhatToDo, TuningParameters & param, char* paramFileName,
	  vector<VehicleProperties> & vehPropVec, char* vehProp_fname){
  
  InOut inout;
  double paramArray[1024];
  int typeID[100];
  double length[100];
  double width[100];

  int n_array;


  // vehicle properties
  
  inout.get_col(vehProp_fname, 2, n_array, typeID);
  inout.get_col(vehProp_fname, 3, n_array, length);
  inout.get_col(vehProp_fname, 4, n_array, width);
  for(int i=0; i<n_array; i++){
    VehicleProperties vehProp;
    vehProp.type_ID=typeID[i];
    vehProp.length=length[i];
    vehProp.width=width[i];
    vehProp.type_str=(typeID[i]==0) ? "motorcycle"
      : (typeID[i]==1) ? "car"
      : (typeID[i]==2) ? "medveh"
      : (typeID[i]==3) ? "truck"
      : (typeID[i]==4) ? "taxi"
      : (typeID[i]==5) ? "bus" : "trafficLight";
    vehPropVec.push_back(vehProp);
  }

  cout<<"init: vehPropVec[1].type_str="<<vehPropVec[1].type_str<<endl;


  
  // global parameters

  inout.get_array (paramFileName, n_array, paramArray);
 
  param.rotAngle= (WhatToDo<=2) ? 0 : paramArray[0];

  bool filter= (WhatToDo<=2) ? paramArray[1]
    : (WhatToDo==3) ? paramArray[2] : 0;
  param.filterForCars=(filter<1e-10) ? false : true;

  filter=(WhatToDo<=2) ? paramArray[3]
    : (WhatToDo==3) ? paramArray[4] : 0;
  param.filterForOrigin=(filter<1e-10) ? false : true;

  param.origin_ymin=(param.rotAngle==0) ? paramArray[5] : paramArray[6];
  param.origin_ymax=(param.rotAngle==0) ? paramArray[7] : paramArray[8];

  filter=(WhatToDo<=2) ? paramArray[9]
    : (WhatToDo==3) ? paramArray[10] : 0;
  param.filterForDest=(filter<1e-10) ? false : true;
  
  param.dest_ymin=(param.rotAngle==0) ? paramArray[11] : paramArray[12];
  param.dest_ymax=(param.rotAngle==0) ? paramArray[13] : paramArray[14];

  filter=(WhatToDo<=2) ? paramArray[15]
    : (WhatToDo==3) ? paramArray[16] : paramArray[17];
  param.filterForXY=(filter<1e-10) ? false : true;
  
  param.xmin=(param.rotAngle==0) ? paramArray[18]
    : (WhatToDo<=4) ? paramArray[19] : paramArray[20];
  param.xmax=(param.rotAngle==0) ? paramArray[21]
    : (WhatToDo<=4) ? paramArray[22] : paramArray[23];
  param.ymin=(param.rotAngle==0) ? paramArray[24] : paramArray[25];
  param.ymax=(param.rotAngle==0) ? paramArray[26] : paramArray[27];

  param.dit=int(paramArray[28]);

// WhatToDo=1: for representation data_s and corresp output

  param.ds              =paramArray[29];
  param.dsSmooth_heading=paramArray[30];

// WhatToDo=2: eadmap and corresp output

  param.gridsize=paramArray[31];

// WhatToDo=3: for purely data-based road and lane finding/matching
  
  param.distBetweenCuts=paramArray[32];
  param.dsSmooth_grid  =paramArray[33];
  param.maxHeadMismatch=paramArray[34];
  param.dsSmooth_lane  =paramArray[35];
  param.wLane_min      =paramArray[36];
  param.maxSinHead     =paramArray[37];
  param.ncutmin    =int(paramArray[38]);
  
  // WhatToDo=4 and 5: logical trajectories

  param.wLane=paramArray[39];
  
  // WhatToDo=5: find neighbors and direct CF leaders

  param.dxmax   =paramArray[40];  
  param.dymax   =paramArray[41];  
  param.TL_dxmin=paramArray[42];
  param.TL_bsafe=paramArray[43];
  param.TL_bcrit=paramArray[44];


  if(true){
    for(int i=0; i<n_array; i++){
      cout<<"i="<<i<<" paramArray[i]="<<paramArray[i]<<endl;
    }
    cout<<"\n\nparam.rotAngle="<<param.rotAngle
	<<"\nparam.filterForCars="<<param.filterForCars
	<<"\n\nparam.filterForOrigin="<<param.filterForOrigin
	<<"\nparam.origin_ymin="<<param.origin_ymin
	<<"\nparam.origin_ymax="<<param.origin_ymax
	<<"\nparam.filterForDest="<<param.filterForDest
	<<"\nparam.dest_ymin="<<param.dest_ymin
	<<"\nparam.dest_ymax="<<param.dest_ymax
	<<"\n\nparam.filterForXY="<<param.filterForXY
	<<"\nparam.xmin="<<param.xmin
	<<"\nparam.xmax="<<param.xmax
	<<"\nparam.ymin="<<param.ymin
	<<"\nparam.ymax="<<param.ymax
	<<"\n\nparam.dit="<<param.dit
	<<"\n\nparam.ds="<<param.ds
	<<"\nparam.dsSmooth_heading="<<param.dsSmooth_heading
	<<"\n\nparam.gridsize="<<param.gridsize
	<<"\n\nparam.distBetweenCuts="<<param.distBetweenCuts
	<<"\nparam.dsSmooth_grid="<<param.dsSmooth_grid
        <<"\nparam.maxHeadMismatch="<<param.maxHeadMismatch
        <<"\nparam.dsSmooth_lane="<<param.dsSmooth_lane 
        <<"\nparam.wLane_min="<< param.wLane_min     
        <<"\nparam.maxSinHead="<<param.maxSinHead     
        <<"\nparam.ncutmin="<< param.ncutmin   
	<<"\n\nparam.wLane="<<param.wLane
	<<"\n\nparam.dxmax="<<param.dxmax
	<<"\nparam.dymax="<<param.dymax
	<<"\nparam.TL_dxmin="<<param.TL_dxmin
	<<"\nparam.TL_bsafe="<<param.TL_bsafe
	<<"\nparam.TL_bcrit="<<param.TL_bcrit
	<<endl<<endl;
  }
}


//#####################################################
// get relative Cartesian North and East coordinates
// from the WGS84 coordinates  (input in degrees)
//#####################################################
double get_North(double lat){
  return radiusEarth*(lat-lat0)*PI/180.;
}

double get_East(double lon){
  return fact_lon*radiusEarth*(lon-lon0)*PI/180.;
}

//#####################################################
// interpolation routine
// argument is fractional index x, extrapolates as constant values
//#####################################################

double intp(const vector<double>ydata, double x){
  if(false){
    cout<<"intp: fract index x="<<x<<" ydata.size()="<<ydata.size()
	<<endl;
  }
  if(x<=0){return ydata[0];}
  else if( int(x)>=int(ydata.size()-1)){return ydata[ydata.size()-1];}
  else{
    int ix=int(x);
    double r=x-ix;
    return (1-r)*ydata[ix]+r*ydata[ix+1];
  }
}



//#####################################################
// transform string vehicle types into integer
// types=[Motorcycle, Car, Medium Vehicle, Heavy Vehicle, Taxi, Bus, redTL]
//#####################################################

int str_type2i(string str_type){
  int index=
    (str_type==" Motorcycle") ? 0 :
    (str_type==" Car") ? 1 :
    (str_type==" Medium Vehicle") ? 2 :
    (str_type==" Heavy Vehicle") ? 3:
    (str_type==" Taxi") ? 4:
    (str_type==" Bus") ? 5 :
    (str_type==" redTrafficLight") ? 6 : 100;  // 100 possibly new type
  return index;
}
    


//#####################################################
// parse one line of data containing some metadata
  // and a complete vehicle trajectory 
//#####################################################
  
void parseOneLine(const string line,  int & index, string & str_type,
		  vector<double> & data_one_line){

  stringstream ss(line);
  string substr;

    // get metadata
    
  if(ss.good()){
      getline( ss, substr, ';' );
      index=stoi(substr);
  }
  if(ss.good()){
      getline( ss, substr, ';' );
      str_type=substr; // !! the strings begin with a blanc!
  }

    // skip the metadata traveled_d; avg_speed;

  for(int j=0; j<2; j++){
      if(ss.good()){  
	getline( ss, substr, ';' );
      }
  }

    // store the line data in the double vector

  int j=0;
  while(ss.good()){  
    getline( ss, substr, ';' );
 
      //! NOT substr.empty() or substr.size() because substr
      // still contains line ending (c<r> etc)
      
    if(ss.good()){
      data_one_line.push_back( stod(substr)); // !! NICHT stof (single prec)
      if(false){
	cout<<"parseOneLine: j="<<j
	    <<" data_one_line[j]="<<data_one_line[j]<<endl;
      }
    }
    j++;
  }
  

}//parseOneLine


/*#####################################################
transforms one original trajectory stored in the structure data 
(not rotated, only translated to be relative to  lat0, lon0)
to the logical coordinates of the (directional) road which has the 
data structure of oneLane, i.e., one top-level element of lanesFinal

@param dit:         only every dit's data line recorded (multiples of 4 ms)
@param data:        original data from the csv files w/o interpretation
@param data_s:      distance-based data with heading info
@param IDvecFiltered:  filtered original IDvecs in data_s (orig: 0,1,2..)
@param road:        vector containing tuples {axis_x,axis_y,heat,heading} 
                    for regular x gridpoints (axis must not deviate 
                    more than, say, 10 degrees to the rotated x axis!)
@param rotAngle:    approximate angle [rad] of the road axis 
                    w/respect to traj_orig; must be the same as in data_s

@param xmin,xmax,   filter rectangle in the rotated coords.   
@param ymin,ymax:   (Must include the road axis)

@action:            trajectory as a vector of structs
                    each struct has the member variables
                    {ID, type, t, x, y, vx, vy, ax, ay, along, alat}
                    vx, vy make use of the headings in data_s 
                    (rotated @ rotangle, so they are differential headings)
                    (along, alat directly from the data (rot invariant)
                     and ax=|a|*cos(differential heading) etc

@output:            file proj.trajLogical

//#####################################################
*/



//#####################################################
// get data structure "RoadData" from file
//#####################################################

RoadData getRoadAxisFromFile(string fname_in, int laneIndex){

  RoadData road_out;
  road_out.rotAngle=-9999;
  
  cout <<"getRoadAxisFromFile: reading "<< fname_in<<endl;
  
  ifstream  infile (fname_in);
  if(!infile){
     cerr << "Error opening file " << fname_in << " for reading" << endl;
     exit(-1);
  }

  string line;
  while (getline(infile, line)){


    // skip comments or empty lines, parse the rest

    if((line.size()>0)&& ((isdigit(line[0]))||(line[0]=='-'))){
      stringstream ss(line);
      string substr;

      // extract stored value for rotAngle
      
      if (line.find("rotAngle") != string::npos) {
	if(ss.good()){
          //getline( ss, substr, '\t' ); // or ss>>rotAngle directly
	  //rotAngle=stod(substr);
	  ss>>road_out.rotAngle;
	}
      }

      // populate the x, y, heading vectors of road_out
      
      else{
	vector<double> numbers;
	while(ss.good()){
	  double number;
	  ss>>number;
	  numbers.push_back(number);
	}
	if(numbers.size()!=5){
	  cerr<<"getRoadAxisFromFile: warning: data line contains "
	      <<numbers.size()<<"!=5 numbers. Doing nothing"<<endl;
	}
	else{
	  if(int(numbers[0])==laneIndex){ // filter only for right laneIndex
	    road_out.x.push_back(numbers[1]);
	    road_out.y.push_back(numbers[2]);
	    road_out.coshead.push_back(cos(numbers[3]));
	    road_out.sinhead.push_back(sin(numbers[3]));
	    //cout<<"numbers[1]="<<numbers[1]<<endl;
	  }
	}
      }
    }
  }
  return road_out;
}

      
//#####################################################
// get logical trajectories from file
//#####################################################

vector<Trajectory> getLogicalTrajsFromFile(string fname_in){

  vector<Trajectory> trajs_out;
  Trajectory traj;
  int IDold=-9999; // initialize outside of valid range N+
  int ID;   // IDs used to determine new trajectory

  cout <<"getLogicalTrajsFromFile: reading "<< fname_in<<endl;
  
  ifstream  infile (fname_in);
  if(!infile){
     cerr << "Error opening file " << fname_in << " for reading" << endl;
     exit(-1);
  }

  string line;
  while (getline(infile, line)){

    bool isTL=(line[0]=='-');
    // skip comments or empty lines, parse the rest

    //cout<<"line[0]="<<line[0]<<" isTL="<<( (isTL) ? "true" : "false")<<endl;
    if((line.size()>0)&& (isdigit(line[0])||isTL) ){
      stringstream ss(line);
      string substr;
      vector<double> numbers;
      while(ss.good()){
	  double number;
	  ss>>number;
	  numbers.push_back(number);
      }
      if(numbers.size()!=10){
	  cerr<<"getLogicalTrajsFromFile: warning: data line contains "
	      <<numbers.size()<<"!=10 numbers. Doing nothing"<<endl;
      }
      else{
	ID=int(numbers[0]);
	if(ID==IDold){ // same trajectory
	  traj.ID=ID;
	  traj.type=int(numbers[1]);
	  traj.time.push_back(numbers[2]);
	  traj.x.push_back(numbers[3]);
	  traj.y.push_back(numbers[4]);
	  traj.heading.push_back(numbers[5]);
	  traj.vx.push_back(numbers[6]);
	  traj.vy.push_back(numbers[7]);
	  traj.ax.push_back(numbers[8]);
	  traj.ay.push_back(numbers[9]);
	}
	else{// new trajectory
	  if(traj.time.size()>0){//catch init situation
	    trajs_out.push_back(traj);
	    if(false){
	      cout<<"getLogicalTrajsFrom File: ID="<<IDold
		  <<" type="<<traj.type
		  <<" #dataPoints="<<traj.time.size()<<endl;
	    }
	  }
	  
	  IDold=ID;
	  traj={};  // reset all components of this struct if new traj

	}
      }
     
    }
    //cout<<"getLogicalTrajsFromFile: traj.ID="<<traj.ID<<endl;
  }

  if(false){
    for(unsigned iveh=0; iveh<trajs_out.size(); iveh++){
      if(trajs_out[iveh].type==6){
        cout<<"red TL at x="<<trajs_out[iveh].x[0]<<"found"
	    <<" y="<<trajs_out[iveh].y[0]
	    <<endl;
      }
    }
  }

  return trajs_out;
}//getLogicalTrajsFromFile



/*
//#####################################################
transformation to logical coordinates
based on the y coordinates of the near horizontal road axis with equidistant
x coordinates {xmin, xmin+dx, ...}
argument road has tuples {x,y,heat,coshead,sinhead}
@return: x_logical: always starting near zero and increasing 
                    (rotated if |road heading|>PI/2)
         y_logical: vertical distance to road axis y
                    (to the left in driving direction is positive)
         heading_road: just copied from input, approx PI if reverse
//#####################################################
*/

  
vector<double> logicalCoords(double x, double y, RoadData road){

  int ncut=road.x.size();  

  if(ncut<2){
    cerr<<"logicalCoords: error: road.x.size()="<<ncut
	<<" too small."<<endl;
    return {0,0,0};
  }
  

  double xmin=road.x[0];
  double xmax=road.x[road.x.size()-1];
  double dx=road.x[1]-road.x[0];
   
  double ixPred=(x-xmin)/dx;
  double yAxisPred=intp(road.y,ixPred);
  double cosheading_road=intp(road.coshead,ixPred);
  double sinheading_road=intp(road.sinhead,ixPred);
  double heading_road=atan2(sinheading_road, cosheading_road);
  double yDist2Axis=y-yAxisPred;
  double y_logical=yDist2Axis*cosheading_road;
  double x_logical=x+yDist2Axis*sinheading_road;
  bool rotateBy180=(fabs(heading_road)>0.5*PI);
  x_logical=(rotateBy180) ? xmax-x_logical : x_logical-xmin;
  if(rotateBy180){y_logical -=1;}

  if(false){
    cout<<"logicalCoords: x="<<x<<" y="<<y
	<<"\n ixPred="<<ixPred
	<<" yAxisPred="<<yAxisPred
	<<" yDist2Axis="<<yDist2Axis
	<<" heading_road="<<heading_road
	<<"\n x_logical="<<x_logical
	<<" y_logical="<<y_logical
	<<endl;
  }
  
  return {x_logical, y_logical, heading_road};
}
    
    
/*
//#####################################################
Get all logical trajectories of a data set 
based on the lane-axis data of lane index roadAxisLane
(lane index as stored in the roadData object created with WhatToDo=3,
file .lanes)
and filtered for heading near 0 or pi (reverseDirection=0 or 1, respectively)
The filter is based on the middle point of a given trajectory
depending on value of filterForXY34 in .param file, also filtered
(in the rotated system!) horizontally and vertically

@return: Logical trajectories as vector<Trajectory>
         - trajectory filtered according to right direction 
           (determined at the middle of the respective trajectory)
         - trajectory filtered according to a minimum length: 
           traj.time.size()>min_dataPoints
         - trajectory points filtered to the rectangular region 
           (in the rotated frame)
           defined by .param if param.filterForXY34 is active
         - even if not active or too large (two parallel road same direction)
           restricted in y by local parameter maxdist_logy 
           (max distance from reference axis)
         - virtual trajectories of traffic lights added if the corresponding
           .trafficLights or .trafficLightsReverse file found
         - in case of reverse, the logical frame is rotated by PI
         - regardless of reverse and rotation, x_logical always increases and 
           starts near 0, y=0 at the roadAxisLane and positive to the left

Basis: function logicalCoords(xrot,yrot,road)

//#####################################################
*/


vector<Trajectory> calcLogicalTrajs(TuningParameters param,
				    const vector<vector<double>> data,
				    const vector<vector<double>> data_s,
				    const vector<int>IDvecFiltered,
				    const vector<int> typeVec,
				    const RoadData road,
				    int roadAxisLane,
				    bool reverseDirection,
				    char* fnameTL){
					   
  vector<Trajectory> trajs_out;
  int min_dataPoints=10;
  double maxdist_logy=15;

  //cout<<"calcLogicalTrajs: road.x.size()="<<road.x.size()<<endl;

  double tmin=1e6;
  double tmax=-1e6;
  
  // calcLogicalTraj (2): rotate traj_orig by rotAngle and filter out
  // traj points outside filter rectangle

  const double cosAngle=cos(param.rotAngle);
  const double sinAngle=sin(param.rotAngle);

  // loop over the already filtered trajectories in traj_s but extraction
  // from the original data
  
  for(unsigned iveh_s=0; iveh_s<data_s.size(); iveh_s++){
    Trajectory traj;
    vector<double> vecx;
    vector<double> vecy;
    int iveh=IDvecFiltered[iveh_s];

    traj.ID=iveh;
    traj.type=typeVec[iveh];

    int is=0; // index of data_s[iveh_s] corresponding to time
    
    for(unsigned it=0; it<data[iveh].size()/6; it+=param.dit){
      double time=data[iveh][5+6*it];
      tmin=min(time,tmin);
      tmax=max(time,tmax);
      double x=get_East(data[iveh][1+6*it]);
      double y=get_North(data[iveh][0+6*it]);
      double xrot=cosAngle*x-sinAngle*y;
      double yrot=sinAngle*x+cosAngle*y;
      double heading_s=0;

      //######################################################
      // do only something if filters passed
      //######################################################
      
      if( (!param.filterForXY)||((xrot>=param.xmin)&&(xrot<=param.xmax)
				 &&(yrot>=param.ymin)&&(yrot<=param.ymax))){

	// find index in traj_s[iveh_s] corresponding to time

	while( (unsigned(4*is+0)<data_s[iveh_s].size())
	       &&(data_s[iveh_s][4*is+0]<time)){
	  is++;
	}
	if(unsigned(4*is+3)<data_s[iveh_s].size()){
	  heading_s=data_s[iveh_s][4*is+3];
	}

	
	// do the actual transform to logical coordinates
	// also takes care of reverse heading of road
	
	vector<double>logCoords
	  =logicalCoords(xrot,yrot,road);

	// populate vectors of traj if not too far from road axis

	if(fabs(logCoords[1])<maxdist_logy){
	  double heading_log=heading_s-logCoords[2]; // rel to road axis
	  double coshead=cos(heading_log);
	  double sinhead=sin(heading_log);
	  double speed=data[iveh][2+6*it]/3.6;
	  double acc=sqrt(SQR(data[iveh][3+6*it])+SQR(data[iveh][4+6*it]));
	  traj.time.push_back(time);  // component time of traj
	  traj.x.push_back(logCoords[0]);
	  traj.y.push_back(logCoords[1]);
	  traj.heading.push_back(heading_log);
	  traj.vx.push_back(speed*coshead);
	  traj.vy.push_back(speed*sinhead);
	  traj.ax.push_back(acc*coshead);
	  traj.ay.push_back(acc*sinhead);
	
	  if(false){
	    unsigned its=traj.time.size();
	    cout<<"compare: "
	      <<" alon="<<data[iveh][3+6*it]<<" ax="<<traj.ax[its-1]
	      <<" alat="<<data[iveh][4+6*it]<<" ay="<<traj.ay[its-1]
	      <<endl;
	  }
	} // maxdist_logy distance from road axis
	
      } // filter

    } // it loop


    // push_back to traj_out only if right direction and >min_dataPoints

    if(int(traj.time.size())>min_dataPoints){
      int itcenter=int(0.5*traj.heading.size());
      //cout<<"traj.heading.size()="<<traj.heading.size()<<endl;
      bool rightDirection=(reverseDirection)
        ? (fabs(traj.heading[itcenter])- PI<param.maxHeadMismatch)
        : (fabs(traj.heading[itcenter])<param.maxHeadMismatch);
    
      if(rightDirection){
        trajs_out.push_back(traj);
      }
    }
 
    if(false){
      //if(false&&(traj.time.size()>10)){
	//  if((iveh_s>=35)&&(iveh_s<39)){
      cout<<"\ncalcLogicalTrajs Test traj iveh_s= "<<iveh_s
	  <<" corresponding to iveh="<<iveh
	  <<" #data points="<<traj.time.size()
	  <<endl;
      if(false){
        for(unsigned itf=0; itf<traj.time.size(); itf+=25){
	  cout<<" time="<<traj.time[itf]
	      <<" x_log="<<traj.x[itf]
	      <<" y_log="<<traj.y[itf]
	      <<" heading_log="<<traj.heading[itf]
	      <<" vx="<<traj.vx[itf]
	      <<" vy="<<traj.vy[itf]
	      <<" ay="<<traj.ay[itf]
	      <<endl;
	}
	cout<<endl;
      }
    }
    
  }// iveh_s loop


  //=========================================================
  //calcLogicalTrajs: add traffic lights as virtual vehicles
  //=========================================================

  InOut inout;
  if(!inout.fileExists(fnameTL)){
    cerr<<"calcLogicalTrajs: warning: no file "<<fnameTL<<" found"
	<<" will not add any traffic lights"<<endl;
  }
  else{
    int nTL=inout.getNumberOfDataLines(fnameTL);

    double xTL[nTL];
    double Tcycle[nTL];
    double phaseRedBegin[nTL];
    double duration[nTL];
    inout.get_col(fnameTL, 1, nTL, xTL);
    inout.get_col(fnameTL, 2, nTL, Tcycle);
    inout.get_col(fnameTL, 3, nTL, phaseRedBegin);
    inout.get_col(fnameTL, 4, nTL, duration);

    double dt=DT*param.dit;

    int ID=0;
  
    for(int iTL=0; iTL<nTL; iTL++){
      double tfirstRed=Tcycle[iTL]*int(tmin/Tcycle[iTL])+phaseRedBegin[iTL];
      for(int ip=0; ip<int((tmax-tmin)/Tcycle[iTL]); ip++){ // phase
        for(int ilane=max(0, roadAxisLane-2); ilane<=roadAxisLane+2;
  	  ilane++){ // every lane gets its own TL
	  ID--; // need new negative ID for each TL and each red episode!
          Trajectory trajTL;
          trajTL.ID=ID;  // all ytraffic lights have "vehID"=6
          trajTL.type=6; // reserved for traffic lights
	  double y=(ilane-roadAxisLane)*param.wLane;
          for(int it=0;it<duration[iTL]/dt; it++){
  	    double time=tfirstRed+ip*Tcycle[iTL]+it*dt;
  	    trajTL.time.push_back(time);
  	    trajTL.x.push_back(xTL[iTL]);
  	    trajTL.y.push_back(y);
  	    trajTL.heading.push_back(0);
  	    trajTL.vx.push_back(0);
  	    trajTL.vy.push_back(0);
  	    trajTL.ax.push_back(0);
  	    trajTL.ay.push_back(0);
	  }
	  trajs_out.push_back(trajTL);
	}
      }
    }
  }
  return trajs_out;
}//calcLogicalTrajs


//#####################################################
// find vehicle ID nearest to a given spatiotemp point
// (to find interesting traj)
// difflane=lane to search minus reference lane of logTrajs
//#####################################################

int getIDnearestTo(double x, double t, int difflane, 
		   TuningParameters param,
		   const vector<Trajectory> logTrajs){
  
  double devx_min=1e10;
  int ID=-9999;
  double heading=-9999;
  double xveh=-9999;
  double yveh=-9999;
  for(int iveh=0; iveh<int(logTrajs.size()); iveh++){
    Trajectory traj=logTrajs[iveh];
    if(int(traj.time.size()>1)){
      double tmin=traj.time[0];
      double dt=traj.time[1]-traj.time[0];
      int it_interesting=round((t-tmin)/dt);
      if((it_interesting>=0)&&(it_interesting<int(traj.time.size()))){
	double devx=fabs(traj.x[it_interesting]-x);
	int dlane=int(round(traj.y[it_interesting]/param.wLane));
	if((dlane==difflane)&&(devx<devx_min)){
	  ID=traj.ID;
	  heading=traj.heading[it_interesting];
	  xveh=traj.x[it_interesting];
	  yveh=traj.y[it_interesting];
	  devx_min=devx;
	}
      }
    }
  }
  cout<<"getIDnearestTo(x="<<x<<", t="<<t
      <<"): ID="<<ID<<" deviation in x="<<devx_min
      <<" xveh="<<xveh<<" yveh="<<yveh
      <<" heading="<<heading
      <<endl;
  return(ID);
}

    



//#####################################################
// write to file proj.road<laneIndex>.traj
//#####################################################


void  writeLogicalTrajs(string projName, TuningParameters param,
			const vector<Trajectory> logTrajs,
			int laneIndex, int reverseDirection)
{

  char fname_out[2048];
  sprintf(fname_out,"%s.road%i.traj", projName.c_str(), laneIndex);
  // int(road.x[0]), int(road.y[0]));
  
  cout<<"\nWriting "<<fname_out<<endl;

  ofstream outfile(fname_out);

  outfile<<"#This file was produced by calling\n#extractTraj_pNEUMA "
  <<projName<<" <WhatToDo=4> <laneRefIndex="<<laneIndex<<">"
   <<" <reverse="<<reverseDirection<<">"
  <<"\n#from the road axis data file "<<projName<<".lanes"
	 <<" with lane index "<<laneIndex
	 <<"\n#parameters:"
	 <<"\n#param.rotAngle="<<param.rotAngle
	 <<"\n#param.xmin="<<param.xmin
	 <<"\n#param.xmax="<<param.xmax
	 <<"\n#param.ymin="<<param.ymin
	 <<"\n#param.ymax="<<param.ymax
	 <<"\n#Type: {0=motorcycle,1=car,2=medveh,3=truck,4=taxi,5=bus,6=TL}"
	 <<"\n#ID\ttype\ttime[s]\tx[m]\ty[m]\theading\tvx\tvy\tax\tay"
	 <<endl;
  outfile<<setprecision(3)<<fixed;
  for(unsigned iveh_out=0; iveh_out<logTrajs.size(); iveh_out++){
    for(unsigned it=0; it<logTrajs[iveh_out].time.size(); it++){
      outfile<<logTrajs[iveh_out].ID<<"\t"
	     <<logTrajs[iveh_out].type<<"\t"
	     <<logTrajs[iveh_out].time[it]<<"\t"
	     <<logTrajs[iveh_out].x[it]<<"\t"
	     <<logTrajs[iveh_out].y[it]<<"\t"
	     <<logTrajs[iveh_out].heading[it]<<"\t"
	     <<logTrajs[iveh_out].vx[it]<<"\t"
	     <<logTrajs[iveh_out].vy[it]<<"\t"
	     <<logTrajs[iveh_out].ax[it]<<"\t"
	     <<logTrajs[iveh_out].ay[it]
	     <<endl;
    }
    outfile<<endl;
  }

  


}//calcLogicalTrajs


		    

/*#####################################################
  get/write the local neighbourhood  of the veh with a certain ID_subj
  from the input logTrajs. Since logTrajs already is rotated if the reverse
  direction is true, no "isReverse" parameter needed: The first and last
  small section of logTrajs is ignored (see local params below)
  and the output is only started once the first leader is found

@param param:          global parameters 
                       (constraints xminRot56, xmaxRot56 not used)
@param vehPropVec:     determines encroachment 
                       using the .width and .length attr
@param onlyCF:         if true, only .FCdata file is written, 
                       no general neighborhood file .leaders or .followers
@param noMotoCFtarget: .FCdata file is onl written if not any motorc as target
@param ID_subj:        .FCdata and neighborhood for follower with ID_subj
@param logTrajs:       input: All logical trajectories (no orientation
                       needed because this already in logTrajs)

@param leaderCF_tseries: Reference just included to get access from outside
                         the function for debugging it in the top level

@return:               apart from above reference, 
                       just writes to files .FCdata (.leaders, .followers)
                       creates no data structure outside of function
#####################################################*/

void extractWriteNeighborhood(TuningParameters param,
			      string projName,
			      bool onlyCF, bool noMotoCFtarget, int ID_subj, 
			      int roadAxisLane,
			      const vector<VehicleProperties> vehPropVec,
			      const vector<Trajectory> logTrajs,
			      vector<vehicle> & leaderCF_tseries,
			      bool & writeCF)
{

  double dxStartIgnored=6.1; // 6.1 drop the first points of the 
                             // log trajectories for various reasons
  double dxEndIgnored=1.0;   // Same for the last part of the subject's traj

  double minDurationCF=10;  // [s]
  //!! vehPropVec[i].width: {moto,car,medVeh,truck,taxi,bus,redTL}

  int nveh=logTrajs.size();  // shortcut


  // no trajectory format because the vehicle ID can change
  // leaderCF_tseries taken to the calling environment
  
  vector<vector<vehicle>> relevantLeaders_tseries;
  vector<vector<vehicle>> relevantFollowers_tseries;



  
  // extractWriteNeighborhood (1):
  // check if the logical trajectories contain a trajectory for ID_subj
  
  int iveh_subj=0;
  int success=0;

  while( (!success)&&(iveh_subj<nveh)){
    success=(logTrajs[iveh_subj].ID==ID_subj);
    iveh_subj++;
  }

  if(!success){
    cerr<<"extractWriteNeighborhood: error: ID_subj "<<ID_subj
	<<" not in set of logical trajectories logTrajs"<<endl;
    exit(-1);
  }
  iveh_subj--;
  
  if(false){
    cout <<endl<<"in extractWriteNeighborhood: iveh_subj="<<iveh_subj
         <<" logTrajs[iveh_subj].ID="<<logTrajs[iveh_subj].ID<<endl;
  }


  // extractWriteNeighborhood (2):
  // get relative time offsets of the other trajectories
  // with respect to subject (notice: trajs are guaranteed to be longer
  // than a minimum of 10 data points,
  // "if(traj.time.size()>10){...}" so time[1] exists)

  vector<int> it_offset;
  double dt=logTrajs[iveh_subj].time[1]-logTrajs[iveh_subj].time[0];
  if(dt<1e-10){
    cerr<<"error: dt="<<dt<<" something went wrong in variable logTrajs\n";
    exit(-1);
  }
  int it0_subj=int(logTrajs[iveh_subj].time[0]/dt+0.5);
  for(int iveh=0; iveh<nveh; iveh++){
    //cout<<"...Neighborhood (2): iveh="<<iveh<<" logTrajs[iveh].type="<<logTrajs[iveh].type<<endl;
    it_offset.push_back(int(logTrajs[iveh].time[0]/dt+0.5)-it0_subj);
  }


  // extractWriteNeighborhood (3):
  // the actual time loop
  
  for(unsigned it_subj=0; it_subj<logTrajs[iveh_subj].time.size(); it_subj++){

    vector<vehicle> vehicles;
    vector<vehicle> leaders; 
    vector<vehicle> followers; 
    vehicle leaderCF;
    vehicle leaderCFnone;

    // subject shortcuts
    
    double time=logTrajs[iveh_subj].time[it_subj];
    double x=logTrajs[iveh_subj].x[it_subj];
    double vx=logTrajs[iveh_subj].vx[it_subj];
    double y=logTrajs[iveh_subj].y[it_subj];
    int type=logTrajs[iveh_subj].type;
    int ID=logTrajs[iveh_subj].ID;
    int lane=roadAxisLane+round(y/param.wLane);
    

    
    // extractWriteNeighborhood (4):
    // get same-time snapshots of all vehs except for the subject
    // (within the it_subj loop)

    
      // initialize with no immediate CF leader (CF file needs to be complete)
    
    leaderCFnone.time=time;
    leaderCFnone.x=x+1e4;
    leaderCFnone.y=y+1e4;
    leaderCFnone.vx=0;
    leaderCFnone.vy=0;
    leaderCFnone.ID=-9999;   // outside range, representing "no leader"
    leaderCFnone.type=-9999;   // outside range, representing "no leader"

    leaderCF=leaderCFnone;
  
    for(int iveh=0; iveh<nveh; iveh++){

      vehicle oneVehicle;
      
      int it=it_subj-it_offset[iveh];
      if( (it>=0)
	  &&(it<int(logTrajs[iveh].time.size()))
	  &&(iveh_subj!=iveh)){
	oneVehicle.iveh=iveh; // for debugging, only
	oneVehicle.ID  =logTrajs[iveh].ID;
	oneVehicle.type=logTrajs[iveh].type;
	oneVehicle.time=logTrajs[iveh].time[it];
	oneVehicle.x   =logTrajs[iveh].x[it];
	oneVehicle.y   =logTrajs[iveh].y[it];
	oneVehicle.lane=roadAxisLane+round(oneVehicle.y/param.wLane)-lane;
	oneVehicle.heading=logTrajs[iveh].heading[it];
	oneVehicle.vx   =logTrajs[iveh].vx[it];
	oneVehicle.vy   =logTrajs[iveh].vy[it];
	oneVehicle.ax   =logTrajs[iveh].ax[it];
	oneVehicle.ay   =logTrajs[iveh].ay[it];
	vehicles.push_back(oneVehicle);
      }
    }


     
    // sort the vehicles with respect to increasing longitudinal x

    sort(vehicles.begin(), vehicles.end(), compare_by_x);

    
    // extractWriteNeighborhood (5) (within the it_subj loop):
    // get relevant leaders and followers
    // filter long distance < param.dxmax (parameter)
    // lateral dist < 1.0*param.wLane (filter not lane based)

    for(unsigned i=0; i<vehicles.size(); i++){
      double dx=vehicles[i].x-x;
      double dy=vehicles[i].y-y;

      if( (abs(dy)<=param.dymax) && (abs(dx)<param.dxmax)
	  &&(vehicles[i].ID !=ID)){
	if(dx>0){
	  leaders.push_back(vehicles[i]);
	}
	else{
	  followers.push_back(vehicles[i]);
	}
      }
    }

 
    // select only relevant leaders/followers (lane based)

    vector<vehicle> relevantLeaders;
    bool vehL=false; // leading vehicle found
    bool vehLL=false; // leading left vehicle found
    bool vehLR=false; // leading right vehicle found

    vector<vehicle> relevantFollowers;
    bool vehF=false; // following vehicle found
    bool vehFL=false; // following left vehicle found
    bool vehFR=false; // following right vehicle found

    for(unsigned i=0; i<leaders.size(); i++){
      int dlane=round((leaders[i].y-y)/param.wLane); // !!relative to subj's y
      if(dlane==0){
	if (!vehL){
	  relevantLeaders.push_back(leaders[i]);
	}
	vehL=true;
      }
      
      if(dlane==1){ // left=positive
	if (!vehLL){
	  relevantLeaders.push_back(leaders[i]);
	}
	vehLL=true;
      }
      
       if(dlane==-1){ // right
	if (!vehLR){
	  relevantLeaders.push_back(leaders[i]);
	}
	vehLR=true;
      }
    }


    //!! ACHTUNG: BOESES FOUL mit unsigned:
    // followers.size()-1 ist riesige positive Zahl, wenn followers.size()=0
    // Auf integer casten, wenn Gefahr von negativen Indices besteht!

    for(int i=int(followers.size()-1); i>=0; i--){ // start with nearest
      int dlane=round((followers[i].y-y)/param.wLane);
      if(dlane==0){
	if (!vehF){
	  relevantFollowers.push_back(followers[i]);
	}
	vehF=true;
      }
      
      if(dlane==1){ // left=positive
	if (!vehFL){
	  relevantFollowers.push_back(followers[i]);
	}
	vehFL=true;
      }
      
       if(dlane==-1){ // right
	if (!vehFR){
	  relevantFollowers.push_back(followers[i]);
	}
	vehFR=true;
      }
    }

    
     
    // debug
    
    //bool debug=((it_subj>=100)&&(it_subj<110));
    //bool debug=((x>=100)&&(x<100.5));
    bool debug=false;

    if(debug){
      cout<<"\ncheck leaders and followers: subject x="<<x
	  <<" y="<<y
	  <<" t="<<time<<endl;
      for(unsigned i=0; i<vehicles.size(); i++){
	cout<<"all vehicles: i="<<i
	    <<" ID="<<vehicles[i].ID
	    <<" type="<<vehicles[i].type
	    <<" x="<<vehicles[i].x
	    <<" dx="<<vehicles[i].x-x
	    <<" dy="<<vehicles[i].y-y
	    <<endl;
      }
      for(unsigned i=0; i<leaders.size(); i++){
	cout<<"leaders for vehicle at x="<<x<<": i="<<i
	    <<" ID="<<leaders[i].ID
	    <<" dx="<<leaders[i].x-x
	    <<" dy="<<leaders[i].y-y
	    <<endl;
      }
      for(unsigned i=0; i<followers.size(); i++){
	cout<<"followers: i="<<i
	    <<" ID="<<followers[i].ID
  	    <<" dx="<<followers[i].x-x
	    <<" dy="<<followers[i].y-y
	    <<endl;
      }
      for(unsigned i=0; i<relevantLeaders.size(); i++){
	cout<<"relevant leaders: i="<<i
	    <<" ID="<<relevantLeaders[i].ID
	    <<" dx="<<relevantLeaders[i].x-x
	    <<" dy="<<relevantLeaders[i].y-y
	    <<endl;
      }
      for(unsigned i=0; i<relevantFollowers.size(); i++){
	cout<<"relevant followers: i="<<i
	    <<" ID="<<relevantFollowers[i].ID
	    <<" dx="<<relevantFollowers[i].x-x
	    <<" dy="<<relevantFollowers[i].y-y
	    <<endl;
      }

    }// debug
 
     
    // extractWriteNeighborhood (6) (within the it_subj loop):
    // select a single leaderCF
    // filter criterion: nearest laterally encroaching leader

    int iLeader=-9999;
    for(int i=0; (iLeader==-9999)&&(i<int(leaders.size())); i++){
      double wEncroach=0.5*(vehPropVec[type].width
			    +vehPropVec[leaders[i].type].width);
      if(fabs(y-leaders[i].y)<=wEncroach){
	leaderCF=leaders[i];
	iLeader=i;
      }
      if(false){
	cout<<" i="<<i<<"(y-leaders[i].y)="<<(y-leaders[i].y)
	    <<" wEncroach="<<wEncroach<<" iLeader="<<iLeader<<endl;
      }
    }
  

    //=======================================================
    // extractWriteNeighborhood (7) (within the it_subj loop):
    // special case: yellow/red traffic light further ahead of the leader
    // who just makes it through the (yellow) light
    //=======================================================

 
    // iLeader is index of the immediate leader in the vector of leaders
    // leaderCF initially initialized with "none": OK
    
    bool leaderIsTL=(leaderCF.type==6);  
    if( (!leaderIsTL) && (iLeader>=0) && (int(leaders.size())>iLeader+1)){
      for(unsigned i=iLeader+1; (!leaderIsTL)&&(i<leaders.size()); i++){
	if(false){
	  cout<<" before replacement 2: t="<<time<<" x="<<x
	      <<" i="<<i<<" leaders[i].type="<<leaders[i].type
	      <<" leaderCF.type="<<leaderCF.type<<endl;
	}
	if(leaders[i].type==6){ // is TL
	  leaderIsTL=true;

	  // need to calculate wEncroach again because new possible leader
	  
	  double wEncroach=0.5*(vehPropVec[type].width
			    +vehPropVec[leaders[i].type].width);
	  double distTL=leaders[i].x-x;
	  double distTL_presentLeader=leaders[i].x-leaderCF.x;
	  
	  // criteria for replacing actual leder by red traffic light:
	  // (1) distance distTL to traffic light <distTL_interact or < set minimum distance
	  // (2) actual leader's speed > own speed
	  // (3) leader's distance to TL < distTL_leader_interact (for nearly stopped leaders)

	  double distTL_interact=0.5*SQR(vx)/param.TL_bsafe;
	  double distTL_leader_interact=0.5*SQR(leaderCF.vx)/param.TL_bsafe;

	  if(false){
	    cout<<" before replacement 3: leaderCF.type="<<leaderCF.type<<endl;
	    cout<<"i="<<i<<" leaderIsTL!  t="<<time
		<<" (fabs(y-leaders[i].y))="<<(fabs(y-leaders[i].y))
		<<" wEncroach="<<wEncroach
		<<" vx="<<vx
		<<" leaders[i].vx-vx="<<leaders[iLeader].vx-vx
		<<" distTL="<<distTL
		<<" param.TL_dxmin="<<param.TL_dxmin
		<<endl;
	  }
	  if( (fabs(y-leaders[i].y)<=wEncroach)
	      && (leaderCF.vx>vx)
	      && (distTL<max(param.TL_dxmin, distTL_interact))
	      && (distTL_presentLeader<distTL_leader_interact) ){
	    leaderCF=leaders[i];
	    //cout<<" replaced actual leader by red TL!"<<endl;
	  }
	}
      }
    }
    
  
    //=======================================================
    // extractWriteNeighborhood (8) (within the it_subj loop):
    // the reverse special case: yellow light is too close and should
    // be passed => set it to green and try find the next leader 
    //=======================================================

    leaderIsTL=(leaderCF.type==6);
    double distTLcrit=0.5*SQR(vx)/param.TL_bcrit;
    double distTL=leaderCF.x-x;

    bool debug2=false;
    
    if(debug2&&leaderIsTL){
      cout<<"\nbefore setting too close TL to green: t="<<time<<" x="<<x
	      <<" leaderCF.type="<<leaderCF.type
	  <<" distTL="<<distTL
	  <<" distTLcrit="<<distTLcrit
	  <<endl;
    }

   
    if(leaderIsTL && (distTL<distTLcrit)){ // only safe way is to go on
      //if(false){  // deactivate this component

      if(debug2){cout<<" in reverting to green loop"<<endl;}
      bool foundAnotherLeader=false;
      for(unsigned i=iLeader+1;
	  leaderIsTL&&(!foundAnotherLeader)&&(i<leaders.size()); i++){
	leaderIsTL=(leaders[i].type==6);
	if(!leaderIsTL){
	  double wEncroach=0.5*(vehPropVec[type].width
				+vehPropVec[leaders[i].type].width);

	  if(fabs(y-leaders[i].y)<=wEncroach){
	    leaderCF=leaders[i];
	    foundAnotherLeader=true;
	  }
	}
      }

      if(!foundAnotherLeader){ // then no leader at all
        leaderCF=leaderCFnone;
      }
    }

    
			      
    // debug
    //if(debug){
    if(false){
    //if((it_subj>100)&&(it_subj<110)){
      cout<<"leaderCF for veh ID "<<ID_subj<<":"
	  <<" iveh_subj="<<iveh_subj
	  <<" type="<<leaderCF.type
	  <<" time="<<leaderCF.time
	// <<" time="<<time
	  <<" x="<<x
	  <<" y="<<y
	  <<"   IDlead="<<leaderCF.ID
	  <<" typeLead="<<leaderCF.type
	  <<" dx="<<leaderCF.x-x
	  <<" dy="<<leaderCF.y-y
	  <<endl;
    }


    //=======================================================
    // extractWriteNeighborhood (9) (within the it_subj loop):
    // if needed, complete the immediate CF leader.
    // Populate the vectors of
    // CF-leader and relevant leader and follower time series
    //=======================================================

    
    relevantLeaders_tseries.push_back(relevantLeaders);
    relevantFollowers_tseries.push_back(relevantFollowers);
    leaderCF_tseries.push_back(leaderCF);
    
  } // it_subj time loop
  


  //=======================================================
  // extractWriteNeighborhood (10): Find interesting trajectories
  //=======================================================

  if(true){
    
    // leaderCF_tseries: vector<vehicle>    [it].ID, [it].type. [it].time
    // subjVeh: Trajectory:                 .ID, .type, .time[it]

    Trajectory subjVeh=logTrajs[iveh_subj];

    double approachStop_xStart=0;
    double approachStop_tStart=0;
    double approachStop_vStart=0;
    double approachStop_sStart=0;
    double approachStop_sminStart=40;
    double approachStop_vlmax=3;
    double approachStop_smaxEnd=8;
    double approachStop_vmaxEnd=3;
    int countStop=0;
    bool inEpisodeStop=false;
    
    //cout<<"checking episodes for veh "<<subjVeh.ID<<endl;

    int nt=int(leaderCF_tseries.size());
    for(int it=0; it<nt; it++){
      int type=subjVeh.type;
      int typeLeader=leaderCF_tseries[it].type;
      int ID=subjVeh.ID;
      double dx=leaderCF_tseries[it].x-subjVeh.x[it];
      //double dy=leaderCF_tseries[it].y-subjVeh.y[it]; //!!! dy-based alert
      double gap=(typeLeader<0) ? dx  // type=-9999 => no leader, dx=1e6
        : dx-0.5*(vehPropVec[type].length
			    +vehPropVec[typeLeader].length);

      double vxl=leaderCF_tseries[it].vx;
      double vx=subjVeh.vx[it];

      // start potentially interesting episode
      
      if((gap>approachStop_sminStart)&&(vxl<=approachStop_vlmax)
	 &&(!inEpisodeStop)){
	inEpisodeStop=true;
	approachStop_xStart=subjVeh.x[it];
	//approachStop_yStart=subjVeh.y[it];
	approachStop_tStart=subjVeh.time[it];
	approachStop_sStart=gap;
	approachStop_vStart=vx;
      }

      // update counter

      if(inEpisodeStop){
	countStop++;
      }

      // close interesting episode and check if it really is interesting
      // !! (later possibly also restrict to unique leader)
      // if so, give signal onto the console

      if((inEpisodeStop)&&((vxl>approachStop_vlmax)||(it==nt-1))){
	//if(true){
	if((gap<approachStop_smaxEnd)&&(vx<approachStop_vmaxEnd)
	 &&(countStop*dt>5)
	   &&(subjVeh.x[it]-approachStop_xStart>0.5*approachStop_sminStart)){
	  cout<<"found interesting stop episode for veh "<<ID
	      <<" type"<<type<<" typeLeader="<<typeLeader
	      <<setprecision(3)
	      <<" tStart="<<approachStop_tStart
	      <<" xStart="<<approachStop_xStart
	      <<" xEnd="<<subjVeh.x[it]
	      <<" laneEnd="<<roadAxisLane+round(subjVeh.y[it]/param.wLane)
	    // start partially strange
	    // <<roadAxisLane+round(approachStop_yStart/param.wLane)
	      <<" sStart="<<approachStop_sStart
	      <<" vStart="<<approachStop_vStart
 	      <<" sEnd="<<gap
	      <<" vEnd="<<vx
	      <<" vlEnd="<<vxl
	      <<" duration="<<(countStop*dt)
	      <<endl;
	}
	inEpisodeStop=false;
	countStop=0;
      }
	
      // debug
      if(ID==2212){
	cout
	  <<"ID="<<ID<<" time="<<subjVeh.time[it]
	  <<" approachStop_tStart="<<approachStop_tStart
	  <<" countStop="<<countStop
	  <<" gap="<<gap
	  <<" v="<<vx
	  <<" vcrit-vxl="<<approachStop_vlmax-vxl
	  <<endl;

      }
    }

  }


  
  //=======================================================
  // extractWriteNeighborhood (11): write files
  // leaderCF_tseries is fully populated also if there is no leader
  //=======================================================

  char fname_out[2048];

  Trajectory subjVeh=logTrajs[iveh_subj];
  int nt=int(subjVeh.time.size()); // shortcut

  // restrict time indices to epsiodes with leaders since some
  // measurement problems at the beginning/and
  // (in between, a missing leader is no problem)
  // Sometimes, traj vehicles stop for a long
  // time at the beginning w/o reason. This is caught as well
  
  int it_CFlast=0;
  int it_CFfirst=0;

  if(subjVeh.type!=6){ // no traffic lights as subjets!

    bool leaderFound=false; // determine the first valid point
    for(int it=0; (!leaderFound)&&(it<nt); it++){
      double covered_dist=subjVeh.x[it]-subjVeh.x[0];
      leaderFound=(covered_dist>dxStartIgnored)&&(leaderCF_tseries[it].ID>0);
      if(leaderFound){it_CFfirst=it;}
    }
    
    leaderFound=false; // determine the last valid point
    for(int it=nt-1; (!leaderFound)&&(it>0); it--){
      double rest_dist=subjVeh.x[nt-1]-subjVeh.x[it];
      leaderFound=(rest_dist>dxEndIgnored)&&((leaderCF_tseries[it].ID>0));
      if(leaderFound){it_CFlast=it+1;}// "+1 because loop ends with <it_CFlast
    }
    

    if(false){
      cout<<"  ...Neighborhood: subjID="<<subjVeh.ID
	<<" first CF data point: leaderID="<<leaderCF_tseries[0].ID
	<<" t="<<subjVeh.time[0]
	<<" x="<<subjVeh.x[0]
	<<" v="<<subjVeh.vx[0]
	<<" last data point with true leader:"
	<<" leaderID="<<leaderCF_tseries[it_CFlast-1].ID
	<<" t="<<subjVeh.time[it_CFlast-1]
	<<" x="<<subjVeh.x[it_CFlast-1]
	<<" v="<<subjVeh.vx[it_CFlast-1]
	  <<""<<endl;
    }
  }

  // filter for .FCdata output. Not relevant for other neighbors
  
  writeCF=(subjVeh.type!=6); 
  if(noMotoCFtarget){
    for(int it=it_CFfirst; writeCF&&(it<it_CFlast); it++){
      if(leaderCF_tseries[it].type==0){
	writeCF=false;
      }
    }
  }

  writeCF=writeCF && (int(subjVeh.time.size())>1)
  && (it_CFlast-it_CFfirst>minDurationCF/(subjVeh.time[1]-subjVeh.time[0]));
	
  if(writeCF){ // no traffic lights as subjets!){

    sprintf(fname_out,"%s_road%i_veh%i.FCdata",
	    projName.c_str(), roadAxisLane, ID_subj);

    cout<<"Writing "<<fname_out<<endl;
    //cout<<"it_CFfirst="<<it_CFfirst<<" it_CFlast="<<it_CFlast<<endl;
    ofstream outfile(fname_out);

    outfile<<"#This file was produced by calling\n#extractTraj_pNEUMA " 
	 <<projName<<" (WhatToDo=) 5 "<<ID_subj<<" (param.dxmax=) "<<param.dxmax
	 <<"\n#based on logical trajectories for the"
	 <<"\n#road axis data file "<<projName<<".lanes"
	 <<" with roadAxisLane "<<roadAxisLane
	 <<"\n#Note1: y increasing to the left"
	 <<"\n#Note2: dx and dy always posLeader-posSubject"
	 <<"\n#Note3: s is bumper-to-bumper gap, traj at veh center assumed"
	 <<"\n#Type: {0=motorcycle,1=car,2=medveh,3=truck,4=taxi,5=bus,6=TL}"
         <<"\n#t[s]\tflwType\tx[m]\ty[m]\theading\tvx[m/s]\tvy[m/s]"
         <<"\t\tleadID\tleadTyp\tdx[m]\tdy[m]\tgap[m]\tlead_vx\tlead_vy"
 	 <<endl;
    outfile<<setprecision(3)<<fixed;
    for(int it=it_CFfirst; it<it_CFlast; it++){ // valid range
      double dx=leaderCF_tseries[it].x-subjVeh.x[it];
      double dy=leaderCF_tseries[it].y-subjVeh.y[it];
      int type=subjVeh.type;
      int typeLeader=leaderCF_tseries[it].type;
      double gap=(typeLeader<0) ? dx  // type=-9999 => no leader, dx=1e6
        : dx-0.5*(vehPropVec[type].length
			    +vehPropVec[typeLeader].length);

      outfile<<subjVeh.time[it]<<"\t"
	   <<subjVeh.type<<"\t"
	   <<subjVeh.x[it]<<"\t"
	   <<subjVeh.y[it]<<"\t"
	   <<subjVeh.heading[it]<<"\t"
	   <<subjVeh.vx[it]<<"\t"
	   <<subjVeh.vy[it]<<"\t\t"
	   <<leaderCF_tseries[it].ID<<"\t"
	   <<leaderCF_tseries[it].type<<"\t"
	   <<dx<<"\t"
	   <<dy<<"\t"
	   <<gap<<"\t"
	   <<leaderCF_tseries[it].vx<<"\t"
	   <<leaderCF_tseries[it].vy
	   <<endl;
    }
  }

  if(onlyCF){return;}
  
  
  sprintf(fname_out,"%s_road%i_veh%i.leaders",
	  projName.c_str(), roadAxisLane, ID_subj);
  //================================================================

  cout<<"Writing "<<fname_out<<endl;

  ofstream outfile2(fname_out);

  outfile2<<"#This file was produced by calling\n#extractTraj_pNEUMA " 
	 <<projName<<" <WhatToDo=5> "<<ID_subj<<" "<<param.dxmax
	 <<"\n#based on logical trajectories for the road axis data file "
	 <<projName<<".lanes"
	 <<" with roadAxisLane "<<roadAxisLane
	 <<"\n#t[s]\tsubjTyp\tx[m]\ty[m]"
	 <<"\t\tlead1ID\tlead1Ty\tl1_dx\tl1_dy\tlead2ID\t ... "
	  <<" (ordered for increasing x)"
	 <<endl;
 
  outfile2<<setprecision(3)<<fixed;
  for(int it=0; it<nt; it++){
    outfile2<<subjVeh.time[it]<<"\t"
	   <<subjVeh.type<<"\t"
	   <<subjVeh.x[it]<<"\t"
	   <<subjVeh.y[it]<<"\t\t";
    for(int iveh=0; iveh<int(relevantLeaders_tseries[it].size()); iveh++){
      vehicle leader=relevantLeaders_tseries[it][iveh];
      outfile2<<leader.ID<<"\t"
	     <<leader.type<<"\t"
	     <<leader.x-subjVeh.x[it]<<"\t"
	     <<leader.y-subjVeh.y[it]<<"\t";
    }
    for(int iveh=0; iveh<3-int(relevantLeaders_tseries[it].size()); iveh++){
      outfile2<<"na\tna\tna\tna\t";
    }
    outfile2<<endl;
  }



  sprintf(fname_out,"%s_road%i_veh%i.followers",
	  projName.c_str(), roadAxisLane, ID_subj);
  //================================================================

  cout<<"Writing "<<fname_out<<endl;

  ofstream outfile3(fname_out);

  outfile3<<"#This file was produced by calling\n#extractTraj_pNEUMA " 
	 <<projName<<" <WhatToDo=5> "<<ID_subj<<" "<<param.dxmax
	 <<"\n#based on logical trajectories for the road axis data file "
	 <<projName<<".lanes"
	 <<" with roadAxisLane "<<roadAxisLane
	 <<"\n#t[s]\tsubjTyp\tx[m]\ty[m]"
	 <<"\t\tflw1_ID\tflw1_Ty\tflw1_dx\tflw1_dy\tflw2_ID\t ..."
	  <<" (ordered for increasing x)"
	 <<endl;
 
  outfile3<<setprecision(3)<<fixed;
  for(int it=0; it<nt; it++){
    outfile3<<subjVeh.time[it]<<"\t"
	   <<subjVeh.type<<"\t"
	   <<subjVeh.x[it]<<"\t"
	   <<subjVeh.y[it]<<"\t\t";
    for(int iveh=0; iveh<int(relevantFollowers_tseries[it].size()); iveh++){
      vehicle follower=relevantFollowers_tseries[it][iveh];
      outfile3<<follower.ID<<"\t"
	     <<follower.type<<"\t"
	     <<follower.x-subjVeh.x[it]<<"\t"
	     <<follower.y-subjVeh.y[it]<<"\t";
    }
    for(int iveh=0; iveh<3-int(relevantFollowers_tseries[it].size()); iveh++){
      outfile3<<"na\tna\tna\tna\t";
    }
    outfile3<<endl;
  }

}//extractWriteNeighborhood



/*#####################################################
Write the complete trajectory data in a different format
(each vehicle at each time in its own row)
for better plotting possibilities. Uses the internal
complete trajectory representation "data"

@param projName:  drone data set, e.g., 20181024_d1_0900_0930
@param param:     set of global control parameters/filter/smoothing etc
@param typeVec: integer vehicle type  
                  {0=motorc, 1=car, 2=medVehicle, 3=heavyVeh, 4=taxi, 5=bus}
@param data:      internal data representation

@output:          file of name <projName>.rowPerVehTime_dit<dit>
//#####################################################
*/

void writeRowPerVehTime(string projName, TuningParameters param,
			const vector<int> IDvec,
			const vector<int> typeVec,
			const vector< vector<double> > data){

  cout<<endl
      <<"in writeRowPerVehTime: writing every "
      <<param.dit<<"th data line"<<endl;

  
  char sep='\t'; // column separator

  char fname_out[2048];
  sprintf(fname_out,"%s.rowPerVehTime_dit%i", projName.c_str(), param.dit);
 
  ofstream outfile(fname_out);

  outfile<<"# Original data in other format."
	 <<"\n# Generated by extractTraj_pNEUMA by the function"
	 <<" writeRowPerVehTime("
	 <<projName<<"..."<<endl
	 <<"# with parameters rotAngle="<<param.rotAngle
	 <<", dit="<<param.dit<<endl
	 <<"#\n# type: 0=motorcycle, 1=car, 2=medVehicle,"
	 <<" 3=heavyVeh, 4=taxi, 5=bus" <<endl
	 <<"vehID\ttype\tt[s]\tx[m]\ty[m]\tv[m/s]\taccLon\taccLat"
	 <<endl;

  for(unsigned iveh=0; iveh<data.size(); iveh++){

    for(unsigned it=0; it<data[iveh].size()/6; it++){
      if(it%param.dit==0){
	double x=get_East(data[iveh][1+6*it]);
	double y=get_North(data[iveh][0+6*it]);
        outfile<< IDvec[iveh]<<sep
	       << typeVec[iveh]<<sep
	       << data[iveh][5+6*it]<<sep
	       << setprecision(5)
	       << x<<sep
	       << y<<sep
	       << data[iveh][2+6*it]/3.6<<sep
	       << data[iveh][3+6*it]<<sep
	       << data[iveh][4+6*it]<<endl;
      }
    }
    outfile<<endl;
  }

  cout<<"writeRowPerVehTime: wrote "<<fname_out<<endl;
}// writeRowPerVehTime



/*#####################################################
change the data representation from time-based to distance-based,
i.e., long episodes of stopping will just result in a single data point
This is needed when, e.g., filteriung for heading or do 
lane identification/map matching.



Output format as in writeRowPerVehTime:
Each vehicle at each position on its own line

@param param:     set of global control parameters. For generation of a
                  lane axis, param.rotAngle should be such that the road
                  ends are parallel to the x axis
@param typeVec:   integer vehicle type  
                  {0=motorc, 1=car, 2=medVehicle, 3=heavyVeh, 4=taxi, 5=bus}
@param data:      internal data representation. Will be transformed to data_s

@internal:        Some finetuning vars at the beginning

@return/action:   populated data structure data_s and filtered IDs 
                  (for later use in other functions)
                  data_s[4*is+0]=time
                  data_s[4*is+1]=x
                  data_s[4*is+2]=y
                  data_s[4*is+3]=heading after rotAngle

@file output:          files <projName>.data_s
                        <projName>.heatmap
                        <projName>.peaks
                        <projName>.lanes
//#####################################################
*/

void transform2DistBased(TuningParameters param,
			 const vector<int> typeVec,
			 const vector< vector<double> > data,
			 vector< int> & IDvecFiltered,
			 vector<int> & typeVecFiltered,
			 vector< vector<double> > & data_s){
 
  
  const double cosAngle=cos(param.rotAngle);
  const double sinAngle=sin(param.rotAngle);
  
  Statistics stat;

  int nVeh=data.size();

  for(int iveh=0; iveh<nVeh; iveh++){
    double s=0;
    int is=0;

    double xold_dt=get_East(data[iveh][1+6*0]);
    double yold_dt=get_North(data[iveh][0+6*0]);

    double xlast=get_East(data[iveh][data[iveh].size()-5]);
    double ylast=get_North(data[iveh][data[iveh].size()-6]);

    double xrotold_dt=cosAngle*xold_dt-sinAngle*yold_dt;
    double yrotold_dt=sinAngle*xold_dt+cosAngle*yold_dt;
    double xrotold_ds=xrotold_dt;
    double yrotold_ds=yrotold_dt;

    //double xrotlast=cosAngle*xlast-sinAngle*ylast;
    double yrotlast=sinAngle*xlast+cosAngle*ylast;


    // Apply filters

    bool filterOriginPassed
      =((!param.filterForOrigin)
	|| ((yrotold_dt<param.origin_ymax)&&(yrotold_dt>param.origin_ymin)));

    bool filterDestPassed
      =((!param.filterForDest)
	|| ((yrotlast<param.dest_ymax)&&(yrotlast>param.dest_ymin)));

    bool filterTypePassed
      =((!param.filterForCars)
	|| (typeVec[iveh]==1));

    bool filtersPassed=(filterOriginPassed&&filterDestPassed
			&&filterTypePassed);
    
    if(filtersPassed){

      typeVecFiltered.push_back(typeVec[iveh]);
      IDvecFiltered.push_back(iveh);
    
      vector<double> vehVec; //[time1,x1,y1,heading1,time2, ...]
      vector<double> coshead;
      vector<double> sinhead;
    
      for(unsigned it=1; it<data[iveh].size()/6; it++){
        double x=get_East(data[iveh][1+6*it]);
        double y=get_North(data[iveh][0+6*it]);
        double xrot=cosAngle*x-sinAngle*y;
        double yrot=sinAngle*x+cosAngle*y;
        s+=sqrt(SQR(xrot-xrotold_dt)+SQR(yrot-yrotold_dt)); 

        if(false){
  	cout<<"it="<<it
  	    <<" xrot="<<xrot<<" xrotold_ds="<<xrotold_ds
	    <<" yrot="<<yrot<<" yrotold_ds="<<yrotold_ds
	    <<" s="<<s
	    <<" int(s/param.ds)="<<int(s/param.ds)<<" is="<<is
	    <<endl;
	}

        //if(int(s/param.ds)>is){ // !! 0/0 error if driving back!


	  double dsAct=sqrt(SQR(yrot-yrotold_ds)+SQR(xrot-xrotold_ds));
	if(dsAct>param.ds){
	  //if((dsAct>param.ds) && (yrot>=param.ymin) && (yrot<=param.ymax) ){
  	  double time=data[iveh][5+6*it];
  	  coshead.push_back((xrot-xrotold_ds)/dsAct);  
  	  sinhead.push_back((yrot-yrotold_ds)/dsAct);
  	  if(isnan(((xrot-xrotold_ds)/dsAct))){
  	  cerr<<"error: iveh="<<iveh<<" it="<<it<<" nt="<<data[iveh].size()/6
  	      <<" time="<<time
  	      <<" s="<<s<<" xrot="<<xrot<<" xrotold_ds="<<xrotold_ds
  	      <<" dsAct="<<dsAct
  	      <<endl; exit(-1);
	  }
  	  vehVec.push_back(time);
  	  vehVec.push_back(xrot);
  	  vehVec.push_back(yrot);
  	  vehVec.push_back(0.); // placeholder for heading; inserted later
  	  is=int(s/param.ds);
  	  xrotold_ds=xrot; // inside condition: xrotold_ds; outside: xrotold_dt
  	  yrotold_ds=yrot;
  
  
  	  if(false){
  	//if(iveh==nVeh/2){
  	  cout<<" analyzeDistBased: "
  	      <<" filtered iveh="<<iveh
  	      <<" type="<<typeVec[iveh]
  	      <<" time="<<time
  	      <<" is="<<is
  	      <<" s="<<s
  	      <<" xrot="<<xrot
  	      <<" yrot="<<yrot
  	      <<endl;
	  }
        }
        xrotold_dt=xrot; // inside condition: xrotold_ds; outside: xrotold_dt
        yrotold_dt=yrot;
        
      }
  
      // smooth coshead, sinhead, calc heading and insert into vehVec
  
      int ns=vehVec.size()/4;
      int dnSmooth=int(param.dsSmooth_heading/param.ds);
      vector<double> cosheadSmoothed=stat.smoothTriang(coshead,dnSmooth);
      vector<double> sinheadSmoothed=stat.smoothTriang(sinhead,dnSmooth);
      for(int is=0; is<ns; is++){
        vehVec[4*is+3]=atan2(sinheadSmoothed[is],cosheadSmoothed[is]);
        if(isnan(vehVec[4*is+3])){
  	cerr<<"error: is="<<is
  	    <<"  coshead[is]="<<coshead[is]
  	    <<"  cosheadSmoothed[is]="<<cosheadSmoothed[is]
  	    <<"  sinhead[is]="<<sinheadSmoothed[is]
  	    <<"  sinheadSmoothed[is]="<<sinheadSmoothed[is]
  	    <<"  headingSmoothed=vehVec[4*is+3]="<<vehVec[4*is+3]
  	    <<endl;
  	exit(-1);
        }
  	  
      }
    
      // debug
      
      if(false){
        
        cout<<"\nTest smoothing of headings: dnSmooth="<<dnSmooth<<endl;
        for(int is=0; is<min(ns,20); is++){
  	cout<<"is="<<is
  	    <<"  coshead[is]="<<coshead[is]
  	    <<"  cosheadSmoothed[is]="<<cosheadSmoothed[is]
  	    <<"  sinhead[is]="<<sinheadSmoothed[is]
  	    <<"  headingSmoothed=vehVec[4*is+3]="<<vehVec[4*is+3]
  	    <<endl;
        }
      }
  
      // set vehVec as element of data_s
  	  
      data_s.push_back(vehVec);

    } // filters
  
  } // iveh loop for input data

  cout<<"transform2DistBased finished: "
      <<" data.size()="<<data.size()
      <<" data_s.size()="<<data_s.size()
      <<" typeVecFiltered.size()="<<typeVecFiltered.size()
      <<endl;
  if(data_s.size()<10){
    cerr<<"transform2DistBased: error: data_s.size()="<<data_s.size()
	<<" too few trajectories passed the filters of"
	<<" <projName>.param<WhatToDo>"<<endl;
    exit(-1);
  }

  
} // transform2DistBased





  

// #############################################################
// write distance-based trajectories of data_s to file
// #############################################################
  
void writeDistBased(string projName, TuningParameters param,
		    const vector<int> IDvecFiltered,
		    const vector<int> typeVecFiltered,
		    const vector< vector<double> > data_s){

  

  char fname_out[2048];
  sprintf(fname_out,"%s.data_s", projName.c_str());

  cout<<"writing "<<fname_out<<endl;


  ostringstream title_oss;
  title_oss
    <<"# Trajectories distance based instead of time based."
    <<"\n# This file was produced by calling extractTraj_pNEUMA "
    <<projName
    <<" WhatToDo=1>"
    <<"\n#with the function writeDistBased with parameter param.ds="
    <<param.ds<<" param.dsSmooth_heading="<<param.dsSmooth_heading
    <<"\n#vehID\ttype\tt[s]\txrot[m]\tyrot[m]\theading";

  string titleString=title_oss.str();
  

  ofstream outfile(fname_out);

  outfile<<titleString<<endl;
  for(unsigned iveh=0; iveh<data_s.size(); iveh++){
    //cout<<" writing "<<fname_out<<": iveh="<<iveh<<endl;
    for(unsigned is=0; is<unsigned(data_s[iveh].size()/4); is++){ 
      outfile<<setprecision(5)
	     <<IDvecFiltered[iveh]<<"\t"
	     <<typeVecFiltered[iveh]<<"\t"
	     <<data_s[iveh][4*is+0]<<"\t"
	     <<data_s[iveh][4*is+1]<<"\t"
	     <<data_s[iveh][4*is+2]<<"\t"
	     <<data_s[iveh][4*is+3]<<endl;
    }
    outfile<<endl;
  }
}// writeDistBased


// #############################################################
// calculate heatmap from data_s (better suited than time-based data
// particularly for the "spin heatmap" headingMatrix
// for aplications with horizontal orientation (WhatToDo>=3)
// set excludeNonHoizontalOD=true but not for general WhatToDo=2
// because there is generally no horiz. orientation for this application
// #############################################################

void calculateHeatmap(TuningParameters param,
		      const vector< vector<double> > data_s,
		      bool excludeVerticalOD,
		      vector< vector<int> > & heatmap,
		      vector< vector<double> > & headingMatrix,
		      double & xminCompl, double & yminCompl){

							      
							      
  cout<<" calculateHeatmap(: create heatmap and matrix of headings"
      <<endl;

  double gridsize=param.gridsize;  // just a shortcut

  
  // get range of values of data_s
  
  double xminGrid=1e10;
  double xmaxGrid=-1e10;
  double yminGrid=1e10;
  double ymaxGrid=-1e10;

  for(unsigned iveh=0; iveh<data_s.size(); iveh++){
    for(unsigned is=0; is<unsigned(data_s[iveh].size()/4); is++){ 
      xminGrid=min(data_s[iveh][4*is+1],xminGrid);
      xmaxGrid=max(data_s[iveh][4*is+1],xmaxGrid);
      yminGrid=min(data_s[iveh][4*is+2],yminGrid);
      ymaxGrid=max(data_s[iveh][4*is+2],ymaxGrid);
    }
  }
 
  // omit the "fringe elements" of the grid with size<gridsize
  
  xminCompl=gridsize*int(xminGrid/gridsize); // reference for output
  yminCompl=gridsize*int(yminGrid/gridsize); // reference for output

  double xmaxCompl=gridsize*int(xmaxGrid/gridsize);
  double ymaxCompl=gridsize*int(ymaxGrid/gridsize);

  if(xminGrid>0){xminCompl+=gridsize;} // int(-0.6)=0 => no need for neg vals
  if(yminGrid>0){yminCompl+=gridsize;}
  if(xmaxGrid<0){xmaxCompl-=gridsize;}
  if(ymaxGrid<0){ymaxCompl-=gridsize;}

  // further crop by param dimensions (to avoid confusing boundary
  // points for determining the lanes from the heatpeaks)

  
  if(param.filterForXY&&(param.xmin>xminCompl)){
    xminCompl+=gridsize*int((param.xmin-xminCompl)/gridsize);}
  if(param.filterForXY&&(param.xmax<xmaxCompl)){
    xmaxCompl-=gridsize*int((xmaxCompl-param.xmax)/gridsize);}
  if(param.filterForXY&&(param.ymin>yminCompl)){
    yminCompl+=gridsize*int((param.ymin-yminCompl)/gridsize);}
  if(param.filterForXY&&(param.ymax<ymaxCompl)){
    ymaxCompl-=gridsize*int((ymaxCompl-param.ymax)/gridsize);}


  int nx=int((xmaxCompl-xminCompl+1e-10)/gridsize);
  int ny=int((ymaxCompl-yminCompl+1e-10)/gridsize);

  cout<<"heatmap: nx="<<nx<<" ny="<<ny<<endl;
  if((nx<10)||(ny<10)){
    cerr<<" heatmap: error: too few gridpoints nx="<<nx
	<<" or ny="<<ny<<endl;
    exit(-1);
  }

  
  // calculateHeatmap:
  // heatmaps with separate headingMatrix. Because of discontinuity of heading
  // need to average cos(heading) and sin(heading)
  // and calc avg heading with atan2 again
  // northern: heading>=0; southern: heading<0

  // note: uses the alredy rotated data matrix data_s
  
  
  vector< vector<double> > sumcoshead (nx, vector<double>(ny,0)); 
  vector< vector<double> > sumsinhead (nx, vector<double>(ny,0)); 
  heatmap.resize (nx, vector<int>(ny,0));              // init with zeros  
  headingMatrix.resize (nx, vector<double>(ny,-9999.)); // init with err vals

  for(unsigned iveh=0; iveh<data_s.size(); iveh++){

    for(unsigned is=0; is<data_s[iveh].size()/4; is++){
      double x=data_s[iveh][4*is+1];
      double y=data_s[iveh][4*is+2];
      double heading=data_s[iveh][4*is+3];

      // exclude the boundary data

      if((x>=xminCompl)&&(x<xmaxCompl)&&(y>=yminCompl)&&(y<ymaxCompl)){

	// int() ensures that cells are filled symmetrically
	  
	int ix=int((x-xminCompl)/gridsize);
	int iy=int((y-yminCompl)/gridsize);
	if((ix<0)||(ix>=nx)||(iy<0)||(iy>=ny)){
	    cerr<<" extractWriteHeatmap: at least one index out of bounds:"
		<<" ix="<<ix<<" nx="<<nx<<" iy="<<iy<<" ny="<<ny<<endl;
	}
	else{
	  double sinhead=sin(heading);
	  if( (!excludeVerticalOD) || (fabs(sinhead)<=param.maxSinHead)){
	    heatmap[ix][iy]+=1;
	    sumcoshead[ix][iy]+=cos(heading);
	    sumsinhead[ix][iy]+=sinhead;
	  }
	}
      }
    }
  }

  // calculateHeatmap: create and populate matrix of headings
  // from the sumcoshead and sumsinhead matrices
  // and initiate with out-of-bounds value as code for "no data"

  for(int ix=0; ix<nx; ix++){
    for(int iy=0; iy<ny; iy++){
      if(heatmap[ix][iy]>0){
	double coshead=sumcoshead[ix][iy]/heatmap[ix][iy];
	double sinhead=sumsinhead[ix][iy]/heatmap[ix][iy];
	headingMatrix[ix][iy]=atan2(sinhead,coshead);
	if(isnan(headingMatrix[ix][iy])){
	  cout<<"error: coshead="<<coshead<<" sinhead="<<sinhead
	      <<" heatmap[ix][iy]="<<heatmap[ix][iy]
	      <<" headingMatrix[ix][iy]="<<headingMatrix[ix][iy]
	      <<endl;
	}
      }
    }
  }

  cout<<"Finished calculaeHeatmap: gridsize="<<gridsize
      <<" Note: did not restrict range from param.xmin etc"
      <<"\n xminCompl="<<xminCompl
      <<" xmaxCompl="<<xmaxCompl
      <<" xminCompl+(nx-1)*gridsize="<<xminCompl+(nx-1)*gridsize
      <<"    nx="<<nx
      <<"\n yminCompl="<<yminCompl
      <<" ymaxCompl="<<ymaxCompl
      <<" yminCompl+(ny-1)*gridsize="<<yminCompl+(ny-1)*gridsize
      <<"    ny="<<ny
      <<endl;
  
}

// #############################################################
// write heatmap to proj.heatmap
// columns x, y, heading, heat=traffic
// excludeVerticalOD only changes fname. Written is the actual file
// calculated by calculateHeatmap
// #############################################################

void writeHeatmap(string projName, TuningParameters param, bool excludeVerticalOD,
		  const vector< vector<int> >  heatmap,
		  const vector< vector<double> >  headingMatrix,
		  double xminCompl, double yminCompl){

  double gridsize=param.gridsize;  // just a shortcut
  int nx=heatmap.size();
  int ny=heatmap[0].size();
  
  cout <<"writeHeatmap:"
       <<" xminCompl="<<xminCompl
       <<" xmaxCompl="<<xminCompl+nx*gridsize
       <<" yminCompl="<<yminCompl
       <<" ymaxCompl="<<yminCompl+ny*gridsize
       <<" gridpoints x=xminCompl+(ix+0.5)*gridsize etc"
       <<endl;

  char fname_out[2048];
  sprintf(fname_out,"%s.heatmap%i", projName.c_str(),
	  ((excludeVerticalOD) ? 3 : 2));

  cout<<"writing "<<fname_out<<endl;

 
  ofstream outfileHeat(fname_out);

  outfileHeat
    <<"#This file was produced by calling extractTraj_pNEUMA "
    <<projName<<" <WhatToDo=2>"
    <<"\n#with the function writeHeatmap"
    <<" with parameters"
    <<"\n#param ds="<<param.ds
    <<"\n#params gridsize="<<gridsize
    <<"\n#param.dsSmooth_heading="<<param.dsSmooth_heading
    <<"\n#xc[m]\tyc[m]\tavgHead\tsum\n";

  for(int ix=0; ix<nx; ix++){
    double x=xminCompl+(ix+0.5)*gridsize;
    for(int iy=0; iy<ny; iy++){
      double y=yminCompl+(iy+0.5)*gridsize;
      outfileHeat<<setprecision(6)
	     <<x<<"\t"
	     <<y<<"\t"
	     <<headingMatrix[ix][iy]<<"\t"
	     <<heatmap[ix][iy]
	     <<endl;
    }
    outfileHeat<<endl;
  }
   

  
}//writeHeatmap




  
//###############################################################
// Identify contiguous lanes
// from the heatpeaks of heatmap calculated by
// cutting the heatmap in y and x directions
// !!! WHATCH OUT! Lane identification very sensitive to rotation
// Essentially trust only the contiguous main lane(s) and take it
// for the whole directional road
// (just filtered in WhatToDo=4 as directional option)
//###############################################################


void identifyWriteLanes(string projName, TuningParameters param,
			const vector< vector<int>> heatmap,
			const vector< vector<double>> headingMatrix,
			double xminCompl, double yminCompl){

  
  Statistics stat;

  double gridsize=param.gridsize;  // just a shortcut
  int nx=heatmap.size();
  int ny=heatmap[0].size();
  cout<<"\nIdentifyWriteLanes: Identifying contiguous lanes from heatmap"
      <<" and headingMatrix"<<endl
      <<" size of heatmap: nx="<<nx<<" ny="<<ny<<endl;
  
  vector< vector<double> > lanes_cuty;
  int dn=max(1,int(param.distBetweenCuts/gridsize));
  int dnSmoothGrid=max(1,int(param.dsSmooth_grid/gridsize));
  for(int icut=0; icut<int(nx/dn-0.5); icut++){
    int ix=int(dn*(icut+0.5)); //!! Therefore the "-0.5" in upper limit icut
    double xc=xminCompl+(ix+0.5)*gridsize;
    //bool debug=(icut<3);
    bool debug=false;

    // ============================================================
    // identifyWriteLanes (1): prepare 1d-cut for finding the peaks
    // ============================================================
    
    vector<double> heatPeaks;  // vector of iy values corresp to peaks
    vector<double> heatPeaksConnected;
    vector<double>heatmap_1d_aggr(ny,0);
    vector<double>heatmap_1d_smooth;

    // for a better sample, add the neighboring x gridpoints
    // including the expected direction

    int avg_nx=5; // !! corresp. to half-width of the neighboring gridpoints

    int dix=min(avg_nx, min(avg_nx,nx-1-ix));  // 
    for(int ixs=ix-dix; ixs<=ix+dix; ixs++){
      for(int iy=0; iy<ny; iy++){
	double heading=headingMatrix[ix][iy];
	double diy=tan(heading)*(ix-ixs);
	if(dix<=6){diy=0;}// no calc heading at the extreme points

	// c++ bug: floor function returns double
	// some bias at the boundary. Don't care
	
	int iylow=max(0,int(floor(iy+diy)));  
	int iyhigh=min(ny-1,int(floor(iy+diy))+1);
	double r=iy+diy-iylow;
	heatmap_1d_aggr[iy]+=(1-r)*heatmap[ixs][iylow]+r*heatmap[ixs][iyhigh];
	heatmap_1d_aggr[iy]/=(2*dix+1.);

	if(false){
	  cout<<"(ix-ixs)="<<(ix-ixs)<<" iy="<<iy
	      <<" heading="<<heading<<" diy="<<diy
	      <<" iylow="<<iylow
	      <<" iyhigh="<<iyhigh
	      <<" r="<<r
	      <<endl;
	}
	
      }
    }
    
    heatmap_1d_smooth=stat.smoothTriang(heatmap_1d_aggr,dnSmoothGrid);


    // ============================================================
    // identifyWriteLanes (2): Actually find the peaks
    // filtered to vertical range param.origin_ymin, param.origin_ymax
    // (inside ic (icut) loop)
    // ============================================================

    //!!! add minimum heat
    // (in smoothed data, corresponds to approx ten times real minheat)
    double minheat=0.4; //!!! corresponds to minheatReal approx 15*minheat
    for(int iy=1; iy<ny-1; iy++){

    // get three neighbor heats to find maxima and calc there
    // location by a parabola defined by the three points
    
      double hm=heatmap_1d_smooth[iy-1];
      double h0=heatmap_1d_smooth[iy];
      double hp=heatmap_1d_smooth[iy+1];

      // peak criterion including minimum heat
      // "+1e-6" because otherwise sometimes double comparison error
    
      //if( (h0>hm+1e-6)&&(h0>hp+1e-6) ){ 
      if( (h0>=minheat)&&(h0>hm+1e-6)&&(h0>hp+1e-6) ){
	cout<<"h0="<<h0<<endl;
        double diypeak=0.5*(hp-hm)/(2*h0-hp-hm);

	// in addition to basic fractional di shift diypeak,
	// use weighted average with weightings controlled by diypeak
	// to become more robust
	
	int nHalfAvg=int(0.5*param.wLane_min/gridsize);
	double wFirst=(diypeak<0) ? 1 : 1-2*diypeak;
	double wLast=(diypeak>=0) ? 1 : 1-(-2*diypeak);
	double diypeak2=0;
	double heatsum=0;
	for(int iya=max(iy-nHalfAvg,0);
	    iya<=min(iy+nHalfAvg,int(heatmap_1d_smooth.size()));
	    iya++){
	  
	  double w=(iya==iy-nHalfAvg) ? wFirst
	    : (iya==iy+nHalfAvg) ? wLast : 1;
	  
	  diypeak2+=w*(iya-iy)*heatmap_1d_smooth[iya];
	  heatsum+=w*heatmap_1d_smooth[iya];
	  if(false){
	    cout<<"heatpeaks: diypeak="<<diypeak
		<<" iya-iy="<<iya-iy<<" w="<<w<<endl;
	  }
	}
	diypeak2/=heatsum;
	double iypeak=iy+diypeak2;
	//cout<<"diypeak="<<diypeak<<"diypeak2="<<diypeak2<<endl<<endl;
	double ypeak=yminCompl+(iypeak+0.5)*gridsize;

	
        // filter for y range (addtl to filter y origin and dest)
      
        if((ypeak>param.ymin)&&(ypeak<param.ymax)){ 
	  heatPeaks.push_back(iypeak);

	  if(false){
	  //if(debug){
	  cout<<"peak found:"
	      <<" ix="<<ix
	      <<" iy="<<iy
	      <<" xpeak=xc="<<xc
	      <<" ypeak="<<ypeak
	      <<" h1dsmooth[iy-1]="<<hm
	      <<" h1dsmooth[iy]="  <<h0
	      <<" h1dsmooth[iy+1]="<<hp
	      <<" h[ix][iy-1]="<<heatmap[ix][iy-1]
	      <<" h[ix][iy]="<<heatmap[ix][iy]
	      <<" h[ix][iy+1]="<<heatmap[ix][iy+1]
	      <<endl;
	  }
	  
	}
      }
    }

    if(false){
      cout<<"icut="<<icut<<" ix="<<ix<<" xc="<<xc
          <<" heatPeaks.size()="<<heatPeaks.size()
	//<<" heatPeaksConnected.size()="<<heatPeaksConnected.size()
 	//  <<" lanes_cuty[icut-1].size()="<<lanes_cuty[icut-1].size()
	  <<endl;
    }

    
    // ============================================================
    // identifyWriteLanes (3): Identify contiguous lanes by
    // connecting peaks of neighboring cuts if the heading criterion
    // is fulfilled
    // ============================================================
  
    // initially, there is no connection to test/set up

    if(icut==0){
      heatPeaksConnected=heatPeaks;
    }

    // afterwards, initialize length of heatPeaksConnected to number
    // of previous peaks and fill with -1 for "not connected"

    else{
       heatPeaksConnected=vector<double>(lanes_cuty[icut-1].size(),-1);

      // "pair" the old peaks from the previous cut with the new peaks 
      // in the expected direction
      // pair only if minimum is found (not possible if length of old
      // peaks vector=0) and the minimum is smaller than max_iy_mismatch
      // the old peaks that found no partner
      // notice: polygamy is not excluded; 2:1 matches possible
      // (remaining -1 elements in heatPeaksConnected)
      // denote an ending lane
    
      int ixold=ix-dn;
      vector<bool> peaksPaired(heatPeaks.size(), false); // matched new peaks
      vector<int> ipeakPaired; // test for duplicates

      for(unsigned ipeakold=0; ipeakold<lanes_cuty[icut-1].size();
	  ipeakold++){
	double iymaxold=lanes_cuty[icut-1][ipeakold];
	int iyold=int(iymaxold+0.5);
	double headold=headingMatrix[ixold][iyold];
	double iy_expected=iymaxold+tan(headold)*dn;  // quadratic grid!

	double max_iy_mismatch=dn*abs(tan(headold+param.maxHeadMismatch)
					-tan(headold-param.maxHeadMismatch));

	// find ipeak of minimum deviation for given ipeakold

	double diymin=1e10;
	int ipeakMatch=-1;
        for(unsigned ipeak=0; ipeak<heatPeaks.size(); ipeak++){
	    double iymax=heatPeaks[ipeak]; // double; parabolic approx!
	    double iy_diff=fabs(iy_expected-iymax);
	    if(iy_diff<diymin){
	      diymin=iy_diff;
	      ipeakMatch=ipeak;
	    }
	}

	// get the final matches (also filter out too large headings)

	bool headingOK=
	  (fabs(headold)<param.maxHeadMismatch)
	  ||(fabs(headold)- PI<param.maxHeadMismatch); // opposite dir
       
	if(headingOK&&(ipeakMatch>=0)&&(diymin<max_iy_mismatch)){

	  // check for duplicates
	  
	  bool duplicate=false;
	  for(unsigned ip=0; ip<ipeakPaired.size(); ip++){ //size initially 0
	    if(ipeakMatch==ipeakPaired[ip]){duplicate=true;}
	  }
	  
	  // add at appropriate pos to heatPeaksConnected and push_back
	  // duplicate check vector ipeakPaired
	  
	  if(!duplicate){
	    heatPeaksConnected[ipeakold]=heatPeaks[ipeakMatch];
	    peaksPaired[ipeakMatch]=true; // initialized to false
	    ipeakPaired.push_back(ipeakMatch);
	  }
	}
      
	//debug
	  
	if(debug){
	  double iymax_best=heatPeaks[ipeakMatch];
	  //double tanPlus=tan(headold+param.maxHeadMismatch);
	  //double tanMinus=tan(headold-param.maxHeadMismatch);
	  //double dy_maxMismatch=max_iy_mismatch*gridsize;
	  double yold=yminCompl+(iymaxold+0.5)*gridsize;
	  double y_expected=yminCompl+(iy_expected+0.5)*gridsize;
	  double y        =yminCompl+(iymax_best+0.5)*gridsize;
	  double y_mismatch=diymin*gridsize;
	  if(ipeakold==0){cout<<endl;}
	  cout
	    <<"after pairing: icut="<<icut
	    <<" ipeakold="<<ipeakold
	    <<" ipeakMatch="<<ipeakMatch
	    <<" xc="<<xc
	    //<<" xold="<<xminCompl+(ix-0.5)*gridsize
	    <<" yold="<<yold
	    <<" tanHeadOld="<<tan(headold)
	    <<" y_expected="<<y_expected
	    <<" y="<<y
		//  <<" tanPlus="<<tanPlus
		//  <<" tanMinus="<<tanMinus
	    <<" y_mismatch="<<y_mismatch
		//<<" diymin="<<diymin
		//<<" max_iy_mismatch="<<max_iy_mismatch
	    <<( (diymin<max_iy_mismatch)
		? "  matched!" : "  not matched")
	    <<endl;
	}
      }
    
 	    
      // ============================================================
      // identifyWriteLanes (4): Identify contiguous lanes by
      // the new peaks that found no partner will initiate new lanes
      // => new elements in the vector heatPeaksConnected
      // resulting in the first elements of heatPeaksConnected possibly -1
      // and an increased sie by the number of new peaks that found no
      // old partner

      if(debug){
        cout<<" heatPeaksConnected.size()="<<heatPeaksConnected.size()
	    <<" peaksPaired.size()="<<peaksPaired.size()<<endl;
	for(int ip=0; ip<int(peaksPaired.size()); ip++){
	  cout<<"ip="<<ip<<" peaksPaired[ip]="<<peaksPaired[ip]<<endl;
	}
      }
      
      for(unsigned ipeak=0; ipeak<heatPeaks.size(); ipeak++){
	if(!peaksPaired[ipeak]){
	  heatPeaksConnected.push_back(heatPeaks[ipeak]);
	}
      }

      // debug
      if(debug){

	cout<<"\nheatPeaksConnected after adding new lanes:"<<endl;

	for(unsigned ilane=0; ilane<heatPeaksConnected.size(); ilane++){
	  double y=yminCompl+(heatPeaksConnected[ilane]+0.5)*gridsize;
	  cout<<" xc="<<xminCompl+(ix+0.5)*gridsize
	      <<" y="<<y
	      <<endl;
	  }
      }
      
    } // else branch for icut>0

    lanes_cuty.push_back(heatPeaksConnected);
    //lanes_cuty.push_back(heatPeaks);

  } // icut loop

  // NOTE: filtered for a minimum of param.ncutmin icut points only at the end
  // of block 5 ("loop over ilane")
  
  // ####################################################################
  // identifyWriteLanes (5)
  // bring contiguous lines together by making iline the master index
  // ####################################################################

  // use interpolated iyd from lanes_cuty[icut][ilane]

  if(false){
    cout<<"\nanalyzeWriteDistBased (3c): "
        <<"\nbring contiguous lines together by making iline the master index"
        <<endl<<endl;
  }

  // lanesFinal[ilane][4*icut+0] = x (with grid interpolation)
  // lanesFinal[ilane][4*icut+1] = y 
  // lanesFinal[ilane][4*icut+2] = heat (propto traffic flow)
  // lanesFinal[ilane][4*icut+3] = direction ]-pi,pi]

  
  int ncut=lanes_cuty.size();
  int nLanesMax=lanes_cuty[ncut-1].size(); // because new lanes push_backed
  vector<vector<double>> lanesFinal;
  for(int ilane=0; ilane<nLanesMax; ilane++){
    vector<double> oneLane;
    for(int icut=0; icut<ncut; icut++){
      int ix=int(dn*(icut+0.5));
      double x=xminCompl+(ix+0.5)*gridsize;

      bool insideRange=(x<param.xmax); 
      
      if(ilane<int(lanes_cuty[icut].size())){
	double iyd=lanes_cuty[icut][ilane];
	int iy=int(iyd+0.5);
	if((heatmap[ix][iy]>0.5)&&insideRange){ //!!
	  oneLane.push_back(x);
	  oneLane.push_back(yminCompl+(iyd+0.5)*gridsize);
	  oneLane.push_back(heatmap[ix][iy]);
	  oneLane.push_back(headingMatrix[ix][iy]);
	  if(false){
	    //if(ix==290){
	    cout<<"ilane="<<ilane
		<<" x="<<xminCompl+(ix+0.5)*gridsize
		<<" y="<<yminCompl+(iyd+0.5)*gridsize
		<<" heat="<<heatmap[ix][iy]
		<<" direction="<<headingMatrix[ix][iy]
		<<endl;
	  }
	}
      }
    }

    
    // smooth the resulting ycoords vectors of each lane

    int dnCutSmooth=int(param.dsSmooth_lane/param.distBetweenCuts+0.5);

    int ncut=oneLane.size()/4;

    vector<double>lane_y(ncut);

    for(int icut=0; icut<ncut; icut++){
      lane_y[icut]=oneLane[4*icut+1];
      //cout<<"icut="<<icut<<" lane_y[icut]="<<lane_y[icut]<<endl;
    }
    vector<double>lane_y_smooth=stat.smoothTriang(lane_y,dnCutSmooth);
    if(true){
      cout<<"\nTest final lane smoothing: dnCutSmooth="<<dnCutSmooth<<endl;
      for(int icut=0; icut<ncut; icut++){
	
	cout<<"ilane="<<ilane
	    <<" lane_y[icut]="<<lane_y[icut]
	    <<" lane_y_smooth[icut]="<<lane_y_smooth[icut]
	    <<endl;
      }
    }


    // re-insert smoothed y values
    
    for(int icut=0; icut<ncut; icut++){
      oneLane[4*icut+1]=lane_y_smooth[icut];
    }

    // set "no data" entries back to error sign y=yminCompl
    
    // only add to lanesFinal if more than
    // param.ncutmin data points (Notice: 4 entries per data point!)

    if(int(oneLane.size())>=4*param.ncutmin){lanesFinal.push_back(oneLane);}

  } // loop over ilane


  // correction 1: Sort by increasing y coordinate of first lane datapoint
  // (second entry [1] of the subvector of lanesFinal)
  
  sort(lanesFinal.begin(), lanesFinal.end(), compare_by_secondEntry);
  
  // correction 2: Merge lanes with only one or two points missing
  // that are therefore identified as separate lanes

  double ydist_max=1.;   //[m]
  double maxMissingPoints=2;
  double dxCut=dn*gridsize;
     // xdist_max=dxCut*xdist_maxMissingPoints; gridsize=param.gridsize

  //merge[i]=0: no merging, everything is OK
  // merge[i]=1: lanesFinal[i] merged with lanesFinal[i+1]
  // merge[i]=-1: lanesFinal[i+1] merged with lanesFinal[i]

  int nLanes_beforeMerge=int(lanesFinal.size());
if(nLanes_beforeMerge==0){
  cerr<<"identifyWriteLanes: Warning: too little data;"
      <<" could not identify any heatmap peaks"<<endl
      <<" so no lanes"<<endl;
 }

  vector<int> merge(max(1,nLanes_beforeMerge-1));

  for(int il=0; il<nLanes_beforeMerge-1; il++){
    int nPoints0=int(lanesFinal[il].size())/4;
    int nPoints1=int(lanesFinal[il+1].size())/4;
    double xdist_incr=lanesFinal[il+1][0]-lanesFinal[il][4*(nPoints0-1)];
    double xdist_decr=lanesFinal[il][0]-lanesFinal[il+1][4*(nPoints1-1)];

    double ydist_incr=lanesFinal[il+1][1]-lanesFinal[il][4*(nPoints0-1)+1];
    double ydist_decr=lanesFinal[il][1]-lanesFinal[il+1][4*(nPoints1-1)+1];

    bool possiblyMergeUp=
      ((xdist_incr>0)&&(round(xdist_incr/dxCut)<=maxMissingPoints+1));
    bool possiblyMergeDown=
      ((xdist_decr>0)&&(round(xdist_decr/dxCut)<=maxMissingPoints+1));

    merge[il]=0;
    if(possiblyMergeUp&&(fabs(ydist_incr)<ydist_max)){ merge[il]=1;}
    if(possiblyMergeDown&&(fabs(ydist_decr)<ydist_max)){ merge[il]=-1;}

    if(false){
      cout<<"possibly merging lanes: il="<<il
	<<"\n  lanesFinal[il+1][0]="<<lanesFinal[il+1][0]
	<<" lanesFinal[il][4*(nPoints0-1)]="<<lanesFinal[il][4*(nPoints0-1)]

	<<"\n  xdist_incr="<<xdist_incr
	<<" xdist_decr="<<xdist_decr
	<<"\n  possiblyMergeUp="<<possiblyMergeUp
	<<" possiblyMergeDown="<<possiblyMergeDown
	<<"\n  (xdist_decr>0)="<<(xdist_decr>0)
	<<" round(xdist_decr/dxCut)="<<round(xdist_decr/dxCut)
	<<"\n  merge[il]="<<merge[il]
	<<endl;
    }
  }

  // do the merging
  int nMerge=0;
  for(int il0=0; il0<nLanes_beforeMerge-1; il0++){
    int il=il0-nMerge;
    int nPoints0=int(lanesFinal[il].size())/4;
    int nPoints1=int(lanesFinal[il+1].size())/4;
    double xdist_incr=lanesFinal[il+1][0]-lanesFinal[il][4*(nPoints0-1)];
    double xdist_decr=lanesFinal[il][0]-lanesFinal[il+1][4*(nPoints1-1)];
    if(merge[il0]!=0){
      int dngap=(merge[il0]==1) ? round(xdist_incr/dxCut)
	: round(xdist_decr/dxCut);
      int ilDown=il;
      int ilUp=il+1;
      if(merge[il0]==1){
	cout<<"merge original lane "<<il0
	    <<" at xc="<<lanesFinal[ilUp][0]-dxCut
	    <<" upwards to "<<il0+1
	    <<" actual lane "<<il<<" to "<<il+1
	    <<endl;
	
	for(int igap=1; igap<dngap; igap++){
	  // special case xCut (=first element of a point)
	  lanesFinal[ilDown].push_back(lanesFinal[ilUp][0]
				       -dxCut*(dngap-igap));
	  for(int ipoint=1; ipoint<4; ipoint++){
	    lanesFinal[ilDown].push_back(lanesFinal[ilUp][ipoint]);
	  }
	}
	
	for(int i=0; i<int(lanesFinal[ilUp].size()); i++){
	  lanesFinal[ilDown].push_back(lanesFinal[ilUp][i]);
	}	
      }
      
      else if(merge[il0]==-1){
	cout<<"merge original lane "<<il0+1
	    <<" at xc="<<lanesFinal[ilDown][0]-dxCut
	    <<" downwards to "<<il0
	    <<" actual lane "<<il0-nMerge+1<<" to "<<il0-nMerge
	    <<endl;
	for(int igap=1; igap<dngap; igap++){
	  // special case xCut (=first element of a point)
	  lanesFinal[ilUp].push_back(lanesFinal[ilDown][0]
				       -dxCut*(dngap-igap));
	  for(int ipoint=1; ipoint<4; ipoint++){
	    lanesFinal[ilUp].push_back(lanesFinal[ilDown][ipoint]);
	  }
	}
	for(int i=0; i<int(lanesFinal[ilDown].size()); i++){
	  lanesFinal[ilUp].push_back(lanesFinal[ilDown][i]);
	};
	lanesFinal[ilDown]=lanesFinal[ilUp];  // elementwise copy up to down
      }

      
      for(int il=ilUp; il<int(lanesFinal.size())-1; il++){
	lanesFinal[il]=lanesFinal[il+1];
      }
      nMerge++;
      lanesFinal.resize(nLanes_beforeMerge-nMerge);
    }
  }
  

 
  //##########################################################
  // identifyWriteLanes(6): write proj.peaks
  //##########################################################

  char fname_out[2048];
  sprintf(fname_out,"%s.peaks", projName.c_str());
  cout <<"\nanalyzeWriteDistBased: writing "<<fname_out<<endl;
  ofstream outfilePeaks(fname_out);

  outfilePeaks<<"#This file was produced by calling identifyWritePeaks"
	      <<"\n#from the program extract_pNEUMA "
	      <<projName<<" <WhatToDo=3>"
	      <<"\n#Pparameters:"
	      <<"\n#param.ds="<<param.ds
	      <<"\n#param.gridsize="<<gridsize
              <<"\n#param.dsSmooth_heading="<<param.dsSmooth_heading
              <<"\n#param.dsSmooth_grid="<<param.dsSmooth_grid
              <<"\n#param.distBetweenCuts="<<param.distBetweenCuts
	      <<endl;

  outfilePeaks<<"\n"<<param.rotAngle<<"\t #rotAngle (rotation angle)"<<endl;
  
  outfilePeaks<<"\n#xc[m]\ty[m]\t\tavgHead\t\theat\tix(debug)\n";
  
  
  for(unsigned icut=0; icut<lanes_cuty.size(); icut++){
    int ix=int(dn*(icut+0.5));
    double xc=xminCompl+(ix+0.5)*gridsize;
    
    for(unsigned ilane=0; ilane<lanes_cuty[icut].size(); ilane++){
      double iyd=lanes_cuty[icut][ilane];
      int iy=int(iyd+0.5);
      double y=yminCompl+(iyd+0.5)*gridsize;



      if(heatmap[ix][iy]>0){
        outfilePeaks<<xc<<"\t"
		    <<y<<"    \t"
		    <<headingMatrix[ix][iy]<<"\t"
	            <<heatmap[ix][iy]<<"\t"
	            <<ix
		    <<endl;
      }
    }
    outfilePeaks<<endl;
  }


  
  //##########################################################
  // identifyWriteLanes(7): write proj.lanes
  //##########################################################

  sprintf(fname_out,"%s.lanes", projName.c_str());
  cout <<"\nanalyzeWriteDistBased: writing "<<fname_out<<endl;

  ofstream outfileLanes(fname_out);

  outfileLanes<<"#This file was produced by calling identifyWriteLanes"
	      <<"\n#from the program extract_pNEUMA "
	      <<projName<<" <WhatToDo=3>"
	      <<"\n#Pparameters:"
	      <<"\n#param.ds="<<param.ds
	      <<"\n#param.gridsize="<<gridsize
              <<"\n#param.dsSmooth_heading="<<param.dsSmooth_heading
              <<"\n#param.dsSmooth_grid="<<param.dsSmooth_grid
              <<"\n#param.distBetweenCuts="<<param.distBetweenCuts
	      <<endl;

  outfileLanes<<"\n"<<param.rotAngle<<"\t #rotAngle (rotation angle)"<<endl;
  
  outfileLanes<<"\n#index\txc[m]\ty[m]\t\tavgHead\t\theat\n";
  
  for(unsigned ilane=0; ilane<lanesFinal.size(); ilane++){
    for(unsigned icut=0; icut<(lanesFinal[ilane].size())/4; icut++){
      outfileLanes<<ilane<<"\t"
		  <<lanesFinal[ilane][4*icut+0]<<"\t"
		  <<lanesFinal[ilane][4*icut+1]<<"   \t"// len approx tab
 		  <<lanesFinal[ilane][4*icut+3]<<"   \t"
 		  <<lanesFinal[ilane][4*icut+2]
		  <<endl;
    }
    outfileLanes<<endl;
  }
cout<<"end identifyWriteLanes(.)"<<endl;
}//  identifyWriteLanes


 

  



//#####################################################
int main(int argc, char* argv[]) {
//#####################################################

  //###############################
  // Cmdline parsing
  // ##############################

 
  char fname_in[MAXSTR];
  
  //cout <<"argc="<<argc<<endl;

  if ( !((argc==3)||(argc==4)||(argc==5))){
    cerr <<"INFO: extracts/manipulates the pNEUMA data" <<endl
	 <<"      from a csv file in the pNEUMA format" <<endl;
    cerr <<"\nCalling sequence:"<<endl<<endl
	 <<" extractTraj_pNEUMA <pNEUMA_csv> <whatDToDo> [further params]"
	 <<endl<<endl
	 <<"Cmd-line parameters:"<<endl<<endl
	 <<" pNEUMA_csv:   Input file as downloaded from pNEUMA site"<<endl
	 <<" WhatToDo:     0: reorganizes data in row per veh and time,"<<endl
	 <<"               1: Transform to a distance-based (instead of time),"
	 <<endl
	 <<"               2: calculates/writes heatmap,"<<endl
	 <<"               3: Identify and write heatpeaks and road axes,"<<endl
	 <<"               4: logical trajectories w/respect to a road,"<<endl
	 <<"               5: extract local neighborhood of a subject"<<endl
	 <<"               6: mass-produce .FCdata files for calibration\n"
	 <<"               7: create veh-type dependent lateral densities"<<endl
	 <<" roadAxisLane: [WhatToDo=4-7] lane of the road to take for log coords\n"
	 <<" reverse:      [WhatToDo=4] 0: not; 1: reverse direction"<<endl
	 <<" ID_subj:      [WhatToDo=5] veh ID to get neighbords and leaderCF"<<endl
	 <<" log_x:        [WhatToDo=7] Logical x coord for cross section"<<endl
	 <<endl
	 <<"Examples:"<<endl<<endl
	 <<" extractTraj_pNEUMA 20181024_d8_0900_0930.csv 0"<<endl
	 <<" extractTraj_pNEUMA 20181024_d8_0900_0930.csv 1"<<endl
	 <<" extractTraj_pNEUMA 20181024_d8_0900_0930.csv 2"<<endl
	 <<" extractTraj_pNEUMA 20181024_d8_0900_0930.csv 3"<<endl
	 <<" extractTraj_pNEUMA 20181024_d8_0900_0930.csv 4 2 0"<<endl
	 <<" extractTraj_pNEUMA 20181024_d8_0900_0930.csv 4 4 1"<<endl
	 <<" extractTraj_pNEUMA 20181024_d8_0900_0930.csv 5 2 1004"<<endl
	 <<" extractTraj_pNEUMA 20181024_d8_0900_0930.csv 5 2 1004"<<endl
	 <<" extractTraj_pNEUMA 20181024_d8_0900_0930.csv 5 4 2224"<<endl
	 <<" extractTraj_pNEUMA 20181024_d1_0900_0930.csv 6 2"<<endl
	 <<" extractTraj_pNEUMA 20181024_d8_0900_0930.csv 7 1 300"<<endl
	 <<endl
	 <<endl;
    exit(-1);
  }

  string str_infile=argv[1];
  size_t lastindex = str_infile.find_last_of("."); 
  string projName = str_infile.substr(0, lastindex);

  sprintf(fname_in,"%s",argv[1]); // input file name
  int WhatToDo=atoi(argv[2]);
  if((WhatToDo==4)&&(argc!=5)){
    cerr<<"error: You chose WhatToDo=4 but did not give the"
	  <<"\nroad axis lane to base the logical trajectories on"
	  <<"\nand/or not whether the reverse direction is true"
	  <<"\n(typically for the higher lane indices)"
	  <<"\nExamples:"
	  <<"\nextractTraj_pNEUMA 20181024_d8_0900_0930.csv 4 2 0"
	  <<"\nor"
	  <<"\nextractTraj_pNEUMA 20181024_d8_0900_0930.csv 4 4 1"
	  <<endl;
    exit(-1);
  }

  if((WhatToDo==5)&&(argc!=5)){
    cerr<<"error: You chose WhatToDo=5 but did not give the road axis lane"
	  <<"\nand/or the subject ID"
	  <<"\nExample for neighbors of vehicle ID 1003:"
	  <<"\nextractTraj_pNEUMA 20181024_d1_0900_0930.csv 5 2 1003"
	  <<"\nor"
	  <<"\nextractTraj_pNEUMA 20181024_d8_0900_0930.csv 5 2 1004"
	<<endl;
    exit(-1);
  }

   if((WhatToDo==6)&&(argc!=4)){
    cerr<<"error: You chose WhatToDo=6 but did not give the road axis lane"
	  <<"\nExample for the third lane from the right as axis:"
	  <<"\nextractTraj_pNEUMA 20181024_d1_0900_0930.csv 6 2"
	<<endl;
    exit(-1);
  }

  if((WhatToDo==7)&&(argc!=5)){
    cerr<<"error: You chose WhatToDo=7 but did not give the road axis lane"
	<<"\nand/or the longitudinal coordinate for the cross section"
	  <<"\nExample for the second lane from the right as axis:"
	  <<"\nextractTraj_pNEUMA 20181024_d8_0900_0930.csv 7 1 300"
	<<endl;
    exit(-1);
  }

  if(true){
    cout <<"WhatToDo="<<WhatToDo
	 <<endl;
  }


  

  // ###############################################################
  // data structures
  // ###############################################################


  // whatToDo=0 data structures for original data

  vector<vector<double>> data;
  // data[iveh]: tuples {lat=y, lon=x, v, accLon, accLat, time}
 
  vector<string> str_typeVec;
  vector<int> IDvec;  
  vector<int> typeVec;

  
  // whatToDo=1 data structures for distance-based and possibly filtered data

  vector< vector<double>> data_s; 

  vector< int> IDvecFiltered;   // refer to original IDvec
  vector< int> typeVecFiltered; // size may be < IDvec, typeVec

  // whatToDo=2 data structures for grid data
  
  vector< vector<int> > heatmap;
  vector< vector<double> >  headingMatrix;


  // whatToDo=3 data structures for peaks and lanes
  
  vector< vector<double> > lanesFinal; // 

  // whatToDo=4 data structures for a single road axis (one lane)
  
  RoadData road;

  
  // ###############################################################
  // input
  // ###############################################################


  // input (1): load parameters file

  char param_fname[1024];
  char vehProp_fname[1024];
  // sprintf(param_fname, "%s.%s%i", projName.c_str(), "param", WhatToDo);
  sprintf(param_fname, "%s.param", projName.c_str());
  sprintf(vehProp_fname, "%s", "vehicleProperties.param");

  TuningParameters param;
  vector<VehicleProperties> vehPropVec;
  init(WhatToDo, param, param_fname, vehPropVec, vehProp_fname);


  
  
   // intput (2): read and parse input csv

  if(WhatToDo<5){ // WhatToDo==5-7 only neads logical traj data
    
    cout <<"main: data input from "<< fname_in<<endl;
  
    ifstream  infile (fname_in);
    if(!infile){
       cerr << "Error opening file " << fname_in << " for reading" << endl;
       exit(-1);
    }

    string line;
    getline(infile, line); // header/title will not be parsed

    int i=0;
    
    while (getline(infile, line)){

      int oneID;
      string oneType;
      vector<double>oneTrajectory;
      parseOneLine(line, oneID, oneType, oneTrajectory);
      IDvec.push_back(oneID);
      str_typeVec.push_back(oneType);
      data.push_back(oneTrajectory);


      if(false){
        cout<<"after parseOneLine: i="<<i
  	    <<"  IDvec[i]="<<IDvec[i]
  	    <<"  str_typeVec[i]="<<str_typeVec[i]
   	    <<"  nData_line="<<data[i].size()
  	    <<"  tmin="<<data[i][5]
  	    <<setprecision(9)
  	    <<" lat0="<<data[i][0]
  	    <<" lon0="<<data[i][1]
            <<" y0="<<get_North(data[i][0])
  	    <<" x0"<<get_East(data[i][1])<<endl
  	    <<" lat1="<<data[i][6]
  	    <<" lon1="<<data[i][7]
            <<" y1="<<get_North(data[i][6])
  	    <<" x1"<<get_East(data[i][7])
  	    <<endl;
        //exit(0);
      }
      
      i++;
    }
  
    int nVeh=i;
    
  
    for(int i=0; i<nVeh; i++){
      typeVec.push_back(str_type2i(str_typeVec[i]));
    }
  
    cout <<"main: data input completed: nVeh="<<nVeh
         <<" data.size()="<<data.size()
         <<endl;
  }// WhatToDo<5


  
  // #########################################################
  // WhatToDo=0: Alternative file output separate rows for each veh and time
  // no other data manipulation as thinning out by saving only every
  // dit'th dtime step
  // #########################################################

  if(WhatToDo==0){
    writeRowPerVehTime(projName, param, IDvec, typeVec, data);
    return(0);
  }


  // #########################################################
  // 1: Transform to a distance-based (instead of time-based representation
  // (basis for heatmap and lane/road identification) and write
  // result to file 
  // #########################################################

  if(WhatToDo==1){
    transform2DistBased(param, typeVec,
			data, IDvecFiltered, typeVecFiltered,
			data_s);
    
    writeDistBased(projName, param, IDvecFiltered, typeVecFiltered, data_s);
    return(0);
  }

  // #########################################################
  // 2: Calculate and write heatmap
  // #########################################################

  if(WhatToDo==2){
    double xminCompl, yminCompl;
    
    transform2DistBased(param, typeVec,
			data, IDvecFiltered, typeVecFiltered,
			data_s);

    bool excludeVerticalOD=false; // here, general orientations sensible
    calculateHeatmap(param, data_s, excludeVerticalOD,
		     heatmap, headingMatrix, xminCompl, yminCompl);
    
    writeHeatmap(projName, param, excludeVerticalOD, // excl only chges fname
		 heatmap, headingMatrix, xminCompl, yminCompl);
    return(0);
  }
    

  // #########################################################
  // 3: Identify and write heatpeaks and axes of roads/lanes
  // #########################################################

  if(WhatToDo==3){
    double xminCompl, yminCompl;

    transform2DistBased(param, typeVec,
			data, IDvecFiltered, typeVecFiltered,
			data_s);
    
    bool excludeVerticalOD=true; // only horizontal orientations sensible
    calculateHeatmap(param, data_s, excludeVerticalOD,
		     heatmap, headingMatrix, xminCompl, yminCompl);

    writeHeatmap(projName, param, excludeVerticalOD, // excl only chges fname
		 heatmap, headingMatrix, xminCompl, yminCompl);

    identifyWriteLanes(projName, param, heatmap, headingMatrix,
		       xminCompl, yminCompl);
    
 }

  // #########################################################
  // 4: Identify and write logical trajectories related to a road
  // other filters as when writing heatmap or .lanes file:
  // spatial restricted, but all vehicles
  // #########################################################

  if(WhatToDo==4){


    int roadAxisLane=atoi(argv[3]);
    bool reverseDirection=(atoi(argv[4])==1) ? true : false;
  
    // get road structure from file
    
    char infile[2048];
    sprintf(infile,"%s.%s", projName.c_str(), "lanes");
    RoadData roadAxis=getRoadAxisFromFile(infile, roadAxisLane);


     // get traffic-light info
    
    char fnameTL[2048];
    sprintf(fnameTL,"%s.%s%i", projName.c_str(),"trafficLights",roadAxisLane);
   
    
    // do all the following with overwritten roadAngle

    transform2DistBased(param, typeVec,
			data, IDvecFiltered, typeVecFiltered,
			data_s);

    
    vector<Trajectory> logTrajs=
      calcLogicalTrajs(param, data, data_s,
		       IDvecFiltered, typeVec,
		       roadAxis, roadAxisLane, reverseDirection, fnameTL);

    writeLogicalTrajs(projName, param, logTrajs,
		      roadAxisLane, reverseDirection);

    cout<<"Notice: used latest .lanes file; if needed, run"
	<<" extract_pNEUMA with option WhatToDo=3 to create a new one\n";
    return(0);

 }



  
  // #########################################################
  // 5: get the neighbors for a certain vehicle index for the whole time
  // this veh exists
  // #########################################################

  if(WhatToDo==5){

    int roadAxisLane=atoi(argv[3]);
    int ID_subj  =atoi(argv[4]);


    // get logical trajs from file
    
    char infile[2048];
    sprintf(infile,"%s.road%i.%s", projName.c_str(),
	    roadAxisLane, "traj");
    vector<Trajectory> logTrajs=getLogicalTrajsFromFile(infile);

    //int diffLane=0; // corresponds to, e.g., plot of _road2_lane2
    //int ID=getIDnearestTo(115, 950, diffLane, param, logTrajs); //!!!
    //cout<<"ID="<<ID<<endl;
    
    // calculate and write the neighborhood
    
    bool onlyCF=false;
    bool noMotoCFtarget=false;
    vector<vehicle> leaderCF_tseries;
    bool writeCF;

    // test reverse (lane 4) with e.g., 20181024_d8_0830_0900_veh2244
    
    extractWriteNeighborhood(param, projName, onlyCF, noMotoCFtarget,
			     ID_subj, roadAxisLane, vehPropVec,
			     logTrajs, leaderCF_tseries, writeCF);

    return(0);
 }

  /*#########################################################
   6: mass-produce .FCdata files for calibration
   serial application of extractWriteNeighborhood
   (is economic since the logical trajectory file needs to be
   read only once)
   filters: 
     1. subject has minDistance=300 m
     2. not any of the leaders must be a motorcycle
     3. subject (follower) no motorcycle     "&&(traj.type!=0)"  
     4. subject (follower) no taxi           "&&(traj.type!=4)"  
     5. All the restrictions of the logical trajs (isReverse, 
        possibly boxConstraints)

  #########################################################*/

  if(WhatToDo==6){

    bool onlyCF=true;
    double minDistance=250;
    bool noMotoCFleader=true; // follower filter: no motos or taxis
    // long filter param.dxmax, lat filter param.dymax
      
    int roadAxisLane=atoi(argv[3]);


    // get logical trajs from file
    
    char infile[2048];
    sprintf(infile,"%s.road%i.%s", projName.c_str(),
	    roadAxisLane, "traj");
    vector<Trajectory> logTrajs=getLogicalTrajsFromFile(infile);
    
    // calculate and write the neighborhood

    for (int iveh=0; iveh<int(logTrajs.size()); iveh++){
      Trajectory traj=logTrajs[iveh];
      int ID_subj=traj.ID;
      int nt=int(traj.time.size());
      double distance=traj.x[nt-1]-traj.x[0];
      bool subjFilterPassed=
	(distance>minDistance)  // suffic. long traj (also exludes TL)
	&&(traj.type!=0)        // no motorcycles as subjects
	&&(traj.type!=4)        // no taxis as subjects
	;

      // leaderCF_tseries can an only be evaluated inside procedure
      // => passed as parameter and populated inside; access from outside
      // for informative purposes via reference parameter
      
      if(subjFilterPassed){
	vector<vehicle> leaderCF_tseries;
	bool writeCF;
	extractWriteNeighborhood(param, projName, onlyCF, noMotoCFleader,
				 ID_subj, roadAxisLane, vehPropVec,
				 logTrajs, leaderCF_tseries, writeCF);

	// analyze leaderCF_tseries
	bool debug=false;
	
	if(debug){
  	  vector<int> leaderIDs;
  	  vector<int> leaderTypes;
	  vector<double> leaderChangeTimes;
	  vector<double> subjLongpos;
	  int IDold=-9999; // =ID of leaderCFnone
	  bool isTLold=false;	
	  bool isNoneold=false;	
	  for(int it=0; it<int(leaderCF_tseries.size()); it++){
	    int IDlead=leaderCF_tseries[it].ID;
	    bool isTL=(leaderCF_tseries[it].type==6);
	    bool isNone=(leaderCF_tseries[it].ID==-9999);
	    //if(ID_subj==952){
	    if(false){
	      cout<<"it="<<it<<" IDlead="<<IDlead
	  	  <<" IDlead="<<IDlead
	  	  <<" type="<<leaderCF_tseries[it].type
	  	  <<" isTL="<<isTL
	  	  <<" isTLold="<<isTLold
	  	  <<" isNone="<<isNone
	  	  <<" isNoneold="<<isNoneold
	  	  <<endl;
	    }

	    
	    // only register nontrivial target changes
	    if( (it==0)
	        || ((IDlead!=IDold)&&(IDlead>0)&&(IDold>0)) // regular leaders
	        || (IDlead*IDold<0) // switch between regular and TL/none
	        || (isNone!=isNoneold) // switch between non and not none
	        || (fabs(leaderCF_tseries[it].x-leaderCF_tseries[it-1].x)>10)
	        ){
	      leaderIDs.push_back(IDlead);
	      leaderTypes.push_back(leaderCF_tseries[it].type);
	      leaderChangeTimes.push_back(leaderCF_tseries[it].time);
	      subjLongpos.push_back(traj.x[it]);
	      IDold=IDlead;
	      isTLold=isTL;
	      isNoneold=isNone;
	    }
	  }


	  if(writeCF){ // !!!write metainfo of leaders
	    cout<<"writing CF data of veh "<<ID_subj<<endl;
	    if(false){
	      cout<<"  veh ID="<<ID_subj<<":  leaderIDs={";
	      for(int ic=0; ic<int(leaderIDs.size()); ic++){
	        cout<<leaderIDs[ic]<<" ";
	      }
	      cout<<"}  times={";
	      for(int ic=0; ic<int(leaderChangeTimes.size()); ic++){
	        cout<<setprecision(0)<<fixed<<leaderChangeTimes[ic]<<" ";
	      }
	      cout<<"}  x={";
	      for(int ic=0; ic<int(leaderChangeTimes.size()); ic++){
	        cout<<setprecision(0)<<fixed<<subjLongpos[ic]<<" ";
	      }
	      cout<<"}\n";
	    }
	  }
	} // if(debug)

	
      }
    }

    return(0);
 }

  // #########################################################
  // 7: get vehicle-dependent density profiles by analyzing the
  // corresponding .road<lane>.traj file for a given
  // logical coordinate x_log
  // #########################################################

  if(WhatToDo==7){

    int roadAxisLane=atoi(argv[3]);
    double x_center  =atof(argv[4]); // (will use 5 cross sect around
                                     // to get more data)

    
    // take trajectories at x_log +/- x_range; 2*x_range=vmax*dt
    
    double dt=DT*param.dit;
    double x_range=30*dt/2;  

    
    // define the attributes of the cross section
    // (y is relative to lane axis) and the histogram data container

    const int nx=5;
    const int dx=20; // take nx crosss section x_center +/- (nx/2)*dx
    
    const double y_halfWidth=6;  // lane axis +/- 2 lanes left/right
    const int ny=100;            // number of histogram points
    const double dy=2*y_halfWidth/ny;
    
    const int nTypes=6; // {0=moto,1=car,2=medveh,3=truck,4=taxi,5=bus}
    vector<vector<int>> histograms(ny,vector<int>(nTypes,0)); // last arg init to 0
    

    // get the logical traj data structure from file
    
    char infile[2048];
    sprintf(infile,"%s.road%i.%s", projName.c_str(),
	    roadAxisLane, "traj");
    vector<Trajectory> logTrajs=getLogicalTrajsFromFile(infile);


    // fill the histograms
    
    for(int iveh=0; iveh<int(logTrajs.size()); iveh++){
      Trajectory traj=logTrajs[iveh];
      for(int ix=0; ix<nx; ix++){
	double x_cs=x_center-(nx/2+ix)*dx;
        bool success=false;
	
        for(int it=0; (!success)&&(it<int(traj.time.size())); it++){
  	  success=( (traj.x[it]>x_cs-x_range)&&(traj.x[it]<=x_cs+x_range));
 	  if(success){
	    int iy=int(ny/2 + round(traj.y[it]/dy));  // round return double
	    
	    if((iy>=0)&&(iy<ny)&&(traj.type>=0)&&(traj.type<nTypes)){
	      histograms[iy][traj.type]++;
	    }
	  }
	}
      }
    }
    
    
   
    // get output in file

    char fname_out[2048];
    sprintf(fname_out,"%s.road%i_x%i.hist", projName.c_str(), roadAxisLane,
	    int(round(x_center)));
	    
    cout<<"\nWriting "<<fname_out<<endl;
	    
    ofstream outfile(fname_out);

    outfile
      <<"#This file was produced by calling\n#extractTraj_pNEUMA "
      <<projName
      <<"  [WhatToDo=]7 [roadAxisLane=]"<<roadAxisLane
      <<"  [x_center=]"<<x_center
      <<"\n#Road axis data file: "<<projName<<".lanes";
    
    outfile
      <<"\n#Histograms for x_center="<<x_center
      <<"\n#y\tmoto\tcar\tmedveh\ttruck\ttaxi\tbus"<<endl;
    
    for(int iy=0; iy<ny; iy++){
      outfile
	  <<(iy-ny/2)*dy
	  <<"\t"<<histograms[iy][0]
	  <<"\t"<<histograms[iy][1]
	  <<"\t"<<histograms[iy][2]
	  <<"\t"<<histograms[iy][3]
	  <<"\t"<<histograms[iy][4]
	  <<"\t"<<histograms[iy][5]
	  <<endl;
    }
 
    
    return(0);
  }

  
  return(0);
}




// #########################################################
// #########################################################




  /*
  //###############################
  // Example random numbers (#include "RandomUtils.h" written by Arne
  // ##############################

  cout <<"\nrandom numbers: myRand() starts with fixed seed by default"<<endl;
  double rand= myRand()-0.5; // rand sim G(-1/2,1/2)
  cout <<"\nrandom number rand sim G(-1/2,1/2): rand="<<rand<<endl;
  rand= myRand()-0.5; // rand sim G(-1/2,1/2)
  cout <<"random number rand sim G(-1/2,1/2): rand="<<rand<<endl;

  // Random seed:
  // srand(seed); e.g., srand(42);
  // Arne's function setRandomSeed(); works only OUTSIDE of scripts;
  // otherwise, obviously the starting time of script used!

  //#####################################################
  */


