
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


int main(int argc, char* argv[]) {

  char   projectName[MAXSTR];
  if ( (argc!=5)&&(argc!=6)){
    cerr <<"\nUsage: analyzeCalibr <file> <hist_min> <classWidth> <nClass> [logScaling]"
	 <<endl<<"logScaling=0: linear; =1: logarithmic in the param"<<endl
	 <<"  with <file> containing the calibrated values in the 3rd col"
	 <<"\nExample:"<<endl
	 <<" analyzeCalibr d1_0900_0930_results_IDM_SSE_s.T -0.1 0.2 20\n"
	 <<" analyzeCalibr d1_0900_0930_results_IDM_SSE_s.GOF 0 0.2 30\n"
	 <<" analyzeCalibr d1_0900_0930_results_IDM_SSE_s.GOF 0 0.2 30 0\n"
	 <<" analyzeCalibr d1_0900_0930_results_IDM_SSE_s.GOF -5.1 0.2 30 1\n"
	 <<endl;
    exit(-1);
  }

  sprintf(projectName,"%s",argv[1]);
  double hist_min=atof(argv[2]);
  double classWidth=atof(argv[3]);
  int nClass=atoi(argv[4]);
  bool logScale=false;
  if (argc==6){logScale=((atoi(argv[5])==1)) ? true : false;}


 

  


//#####################################################
// input
//#####################################################

  char   fnameIn[MAXSTR];
  sprintf(fnameIn,"%s",projectName);

  InOut inout;
  int nData=inout.getNumberOfDataLines(fnameIn);
  double calibrValues[nData];

    // outfile   param    calibrVal   error   errorVal

  inout.get_col(fnameIn, 3, nData, calibrValues);

  if(logScale){
    for(int i=0; i<nData; i++){
      calibrValues[i]=(calibrValues[i]>0) ? log(calibrValues[i]) : -9999;
    }
  }

//#####################################################
// Do some calculations
//#####################################################

  Statistics stat;

  
  double histogram[nClass];  // double: compatibility write_array
  stat.calculate_histogram(calibrValues, nData, hist_min, classWidth,
			   nClass, histogram);
  double classCenter[nClass];
  for(int i=0; i<nClass; i++){
    classCenter[i]=hist_min+(i+0.5)*classWidth;
  }
  
//#####################################################
// histogram file Output
//#####################################################

  char   fnameOut[MAXSTR+10];
  char   titleStr[MAXSTR];

  sprintf(fnameOut,"%s.hist",projectName);
  sprintf(titleStr,"%s","#generated by analyzeCalibr\n#value\tnHist");
  if(logScale){
    sprintf(fnameOut,"%s_log.hist",projectName);
    sprintf(titleStr,"%s","#generated by analyzeCalibr\n#value\tnHist 1");
  }

  cout <<" writing to "<<fnameOut<<" ..."<<endl;
  inout.write_array(fnameOut, nClass, classCenter, histogram, titleStr);

  return(0);
}

// #########################################################







