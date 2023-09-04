
/*
  (mar18)
  Template fuer main() Programm; Ort:
  
  ~/versionedProjects/lib/templates/templateMain.cpp

  makefile dazu:
  
  ~/versionedProjects/lib/templates/makefile
 
  Achtung! Auch ohne .h File muss man bei $OBJECTS immer auch das 
  File mit der Main-Methode dazunehmen!
  (sonst "ld: undefined reference to main"
*/

using namespace std;

// c
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

//alternatively there is <cstdio> which declares
//everything in namespace std
//but the explicit "using namespace std;" puts
//everything in global namespace


// c++  (the first is nearly always needed)
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip> //!! <feb19>
#include <vector>

// own (delete if not needed)
#include "general.h"
#include "Statistics.h"
#include "RandomUtils.h" // contains, e.g.,  myRand()
#include "InOut.h"
#include "Math.h"


// constants

static const int NDATA_MAX=5000;// max. number of data points
static const int MAXSTR=500;// max. string length


//#####################################################
//#####################################################


double type2typeID(string type){
  return (type=="bike") ? 0
    : (type=="car") ? 1
    : (type=="truck") ? 2 : 3;
}


void showInfo(){
  
    cerr <<"\nUsage: analyzeDistributions <project><extension>"<<endl
	 <<" extension may be one of FCdata or traj"<<endl
	 <<" Input is <project>.<extension> which may be concatenated"
	 <<" from several .FCdata or .traj files"<<endl
	 <<" output are files such as <project>_<extension.accDistr,"
	 <<" or <project>_<extension>.TTCdistr"<<endl
	 <<endl<<"Examples:"<<endl
	 <<"  analyzeDistributions d8 FCdata"<<endl
	 <<"  analyzeDistributions 20181029_d8_0930_1000.road2 traj"<<endl
	 <<endl <<endl;
}



int main(int argc, char* argv[]) {

//#####################################################
  // cmdline parsing
//#####################################################
  
  if (argc!=3){
    showInfo();
    exit(-1);
  }

  char   projectName[MAXSTR];
  char   extension[10];
  sprintf(projectName,"%s",argv[1]);
  sprintf(extension,"%s",argv[2]);
  string str_extension(extension);
  
  //#####################################################
  // input
  //#####################################################


  InOut inout;
  char   fnameIn[MAXSTR+10];
  sprintf(fnameIn,"%s.%s",projectName, extension);
  int nData=inout.getNumberOfDataLines(fnameIn);
  cout<<"nData="<<nData<<endl;

  vector<double> time;
  vector<double> vx;
  vector<double> lead_vx;
  vector<double> sx;
  vector<int> leadID;

  //t[s] flwType x[m] y[m] hdg vx[m/s] vy[m/s] leadID  leadTyp dx[m] dy[m] gap[m] lead_vx lead_vy

  if(str_extension == "FCdata"){
    inout.get_col(fnameIn, 1, nData, time);
    inout.get_col(fnameIn, 6, nData, vx);
    inout.get_col(fnameIn, 8, nData, leadID);
    inout.get_col(fnameIn, 12, nData, sx);
    inout.get_col(fnameIn, 13, nData, lead_vx);
  }
  else if(str_extension == "traj"){
    inout.get_col(fnameIn, 3, nData, time);
    inout.get_col(fnameIn, 7, nData, vx);
    inout.get_col(fnameIn, 1, nData, leadID); // own ID here but OK for reset
    // inout.get_col(fnameIn, 12, nData, sx);
    sx.resize(nData);
    lead_vx.resize(nData);
    for(int i=0; i<nData; i++){
      sx[i]=-999; // err value since val not given
      lead_vx[i]=-999;
    } 
  }
  else{
    cerr<<"Extension is neither FCdata nor traj";
    showInfo();
    exit(-1);
  }

 
   
  //#####################################################
  // calculations
  //#####################################################

  // accelerations and inverse TTC
  
  double dt=time[1]-time[0];
  double SMALL_VAL=1e-9;
  if((dt<0.02)||(dt>1.0)){
    cerr<<"analyzeDistributions: error:"
	<<" calculated sampling interval not in [0.02,1]"<<endl;
    exit(-1);
  }
  vector<double> accx;
  vector<double> TTC_slow;
  vector<double> TTC_fast;
  vector<double> bkin_slow;
  vector<double> bkin_fast;
  vector<double> invbkin_slow;
  vector<double> invbkin_fast;
  vector<double> critIDM_slow;
  vector<double> critIDM_fast;

  // accelerations
  
  int iValid=0;
  for(int i=1; i<nData-1; i++){
    bool newLeader=((leadID[i+1]!=leadID[i])||(leadID[i]!=leadID[i-1]));

    // needed because of concatenation (no followerID because in orig fnames)
    bool newFollower=((fabs(time[i+1]-time[i]-dt)>SMALL_VAL) ||
		      (fabs(time[i]-time[i-1]-dt)>SMALL_VAL));
    
    bool valid=(vx[i+1]>0) || (vx[i-1]>0); // ignore stopped episodes
    if(valid){
      accx.push_back((newLeader||newFollower) ? 0 : 0.5*(vx[i+1]-vx[i-1])/dt);
      iValid++;
    }
  }
  int nValid_acc=iValid;
  cout<<"accel distributions: nValid_acc="<<nValid_acc<<endl;

  // TTC, bkin and others

  double vc=6.; // threshold between slow and fast
  iValid=0;
  int iValid_slow=0;
  int iValid_fast=0;
  if(str_extension == "FCdata"){
    for(int i=1; i<nData-1; i++){
      bool valid=((vx[i]-lead_vx[i])>0); // ignore stopped episodes
      if(valid){
	iValid++;
	if(vx[i]<=vc){
          TTC_slow.push_back(sx[i]/(vx[i]-lead_vx[i]));
          bkin_slow.push_back(0.5*pow(vx[i]-lead_vx[i],2)/sx[i]);
          invbkin_slow.push_back(2*sx[i]/pow(vx[i]-lead_vx[i],2));
          critIDM_slow.push_back(sx[i]/(vx[i]*(vx[i]-lead_vx[i])));
	  iValid_slow++;
	}
	else{
          TTC_fast.push_back(sx[i]/(vx[i]-lead_vx[i]));
          bkin_fast.push_back(0.5*pow(vx[i]-lead_vx[i],2)/sx[i]);
          invbkin_fast.push_back(2*sx[i]/pow(vx[i]-lead_vx[i],2));
          critIDM_fast.push_back(sx[i]/(vx[i]*(vx[i]-lead_vx[i])));
	  iValid_fast++;
	}
	
      }
    }
  }
  int nValid_TTC_slow=iValid_slow;
  int nValid_TTC_fast=iValid_fast;
  int nValid_TTC=iValid;
  cout<<"TTC distributions: nValid_TTC_slow="<<nValid_TTC_slow
      <<" nValid_TTC_fast="<<nValid_TTC_fast
      <<" nValid_TTC="<<nValid_TTC
      <<endl;

  

  // histograms 
  
  Statistics stat;
  int nClass_acc=2*int(0.5*0.2*sqrt(nValid_acc))+1;
  int nClass_TTC=2*int(0.5*0.1*sqrt(nValid_TTC))+1;

  // histogram acceleration ax
  
  double hist_acc[nClass_acc];  // double: compatibility write_array
  double accmin=-5;
  double accmax=5;
  double classwidth_acc=(accmax-accmin)/nClass_acc;
  stat.calculate_histogram(accx, accmin, classwidth_acc,
			   nClass_acc, hist_acc);
  // restrict count of center column (corresp. to acc=0)
  // to the max of the neighboring columns

  int countSidesMax=max(hist_acc[nClass_acc/2-1], hist_acc[nClass_acc/2+1]);
  hist_acc[nClass_acc/2]=countSidesMax;


  // histograms TTC, bkin and others (in x direction)

  double hist_TTC_slow[nClass_TTC];  // nClass_TTC may be =0;
  double hist_TTC_fast[nClass_TTC];  
  double hist_bkin_slow[nClass_TTC]; 
  double hist_bkin_fast[nClass_TTC];  
  double hist_invbkin_slow[nClass_TTC]; 
  double hist_invbkin_fast[nClass_TTC];  
  double hist_critIDM_slow[nClass_TTC]; 
  double hist_critIDM_fast[nClass_TTC];  
  double TTCmin=0.;
  double TTCmax=10;
  double bkinmin=0.;
  double bkinmax=3;
  double invbkinmin=0.;
  double invbkinmax=5;
  double critIDMmin=0.;
  double critIDMmax=5;
  double classWidth_TTC=(TTCmax-TTCmin)/nClass_TTC;
  double classWidth_bkin=(bkinmax-bkinmin)/nClass_TTC;
  double classWidth_invbkin=(invbkinmax-invbkinmin)/nClass_TTC;
  double classWidth_critIDM=(critIDMmax-critIDMmin)/nClass_TTC;
  
  if(str_extension == "FCdata"){
    stat.calculate_histogram(TTC_slow, TTCmin, classWidth_TTC,
			     nClass_TTC, hist_TTC_slow);
    stat.calculate_histogram(TTC_fast, TTCmin, classWidth_TTC,
			     nClass_TTC, hist_TTC_fast);
    stat.calculate_histogram(bkin_slow, bkinmin, classWidth_bkin,
			     nClass_TTC, hist_bkin_slow);
    stat.calculate_histogram(bkin_fast, bkinmin, classWidth_bkin,
			     nClass_TTC, hist_bkin_fast);
    stat.calculate_histogram(invbkin_slow, invbkinmin, classWidth_invbkin,
			     nClass_TTC, hist_invbkin_slow);
    stat.calculate_histogram(invbkin_fast, invbkinmin, classWidth_invbkin,
			     nClass_TTC, hist_invbkin_fast);
    stat.calculate_histogram(critIDM_slow, critIDMmin, classWidth_critIDM,
			     nClass_TTC, hist_critIDM_slow);
    stat.calculate_histogram(critIDM_fast, critIDMmin, classWidth_critIDM,
			     nClass_TTC, hist_critIDM_fast);
  }


 


  //#####################################################
  // output
  //#####################################################


  // acc distributions
  
  double classCenter[nClass_acc]; // class center array for output
  for(int i=0; i<nClass_acc; i++){
    classCenter[i]=accmin+(i+0.5)*classwidth_acc;
  }


  char   fnameAccDistr[MAXSTR+9];
  sprintf(fnameAccDistr,"%s.accDistr",projectName);
  char   titleStr[MAXSTR];

  sprintf(titleStr,"%s","#generated by analyzeDistributions\n#value\tnHist");

  cout <<" writing to "<<fnameAccDistr<<" ..."<<endl;
  inout.write_array(fnameAccDistr, nClass_acc, classCenter, hist_acc,
		    titleStr);

  // TTC distributions
  
  if(str_extension == "FCdata"){

    double classCenter_TTC[nClass_TTC]; // class center array for output
    double classCenter_bkin[nClass_TTC]; // class center array for output
    double classCenter_invbkin[nClass_TTC]; // class center array for output
    double classCenter_critIDM[nClass_TTC]; // class center array for output

    for(int i=0; i<nClass_TTC; i++){
      classCenter_TTC[i]=TTCmin+(i+0.5)*classWidth_TTC;
      classCenter_bkin[i]=bkinmin+(i+0.5)*classWidth_bkin;
      classCenter_invbkin[i]=invbkinmin+(i+0.5)*classWidth_invbkin;
      classCenter_critIDM[i]=critIDMmin+(i+0.5)*classWidth_critIDM;
    }


    
    char   fnameTTCdistr[MAXSTR+9];
    sprintf(fnameTTCdistr,"%s.TTCdistr",projectName);

    sprintf(titleStr,"#generated by analyzeDistributions, speed threshold vc=%.1f\n#valTTC\tn_TTC_slow\tn_TTC_fast\tvalbkin\tn_bkin_slow\tn_bkin_fast\tvalinvbkin\tn_invbkin_slow\tn_invbkin_fast\tvalcritIDM\tn_critIDM_slow\tn_critIDM_fast",vc);

    cout <<" writing to "<<fnameTTCdistr<<" ..."<<endl;
    inout.write_array(fnameTTCdistr, nClass_TTC,
		      classCenter_TTC, hist_TTC_slow, hist_TTC_fast,
		      classCenter_bkin, hist_bkin_slow, hist_bkin_fast,
		      classCenter_invbkin, hist_invbkin_slow, hist_invbkin_fast,
		      classCenter_critIDM, hist_critIDM_slow, hist_critIDM_fast,
		      titleStr);
  }
 
  return(0);
}

// #########################################################







