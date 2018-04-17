#include "vectorfileread.hxx"

#include<vector>
#include<string>
#include<fstream>
#include<iostream>
#include<sstream>

using namespace std;

#include <cstdlib>

namespace FRC {
  VOL_T vectorfileread(const string &filename) {
    ifstream infile(filename.c_str());
    VOL_T data;
    while(infile) {
      FIBRE_T fibre;
      string s;
      if (!getline(infile,s)) break;
      istringstream ss_s(s);
      while(ss_s) {
	string x;
	if (!getline(ss_s,x,' ')) break;
	istringstream ss_x(x);
	PNT_T pnt;
	while(ss_x) {
	  string element;
	  if (!getline(ss_x,element,',')) break;
	  pnt.push_back(atof(element.c_str()));
	}
	fibre.push_back(pnt);
      }
      if (!fibre.empty())
	data.push_back(fibre);
    }
    return data;
  }

  void vectorfilewrite(const string &filename, const VOL_T &fibres) {
    ofstream fileout(filename.c_str());
    int i,j,k;
    for (i=0;i<fibres.size();i++) {
      //      std::cout <<"#allfibres:"<<fibres.size()<<std::endl;
      for (j=0;j<fibres[i].size();j++) {
	//	std::cout <<"#pnts:"<<fibres[i].size()<<std::endl;
	for (k=0;k<fibres[i][j].size()-1;k++) {
	  //	  std::cout <<fibres[i][j][k]<<",";
	  fileout <<fibres[i][j][k]<<",";	  
	}
	//	std::cout <<fibres[i][j][k]<<" ";
	fileout <<fibres[i][j][k]<<" ";
      }
      fileout <<std::endl;
      //      std::cout <<std::endl;
    }
    fileout.close();
    //    std::cout <<"Completed vectorfilewrite"<<std::endl;
  }

  void vectorfilewrite(const string &filename, const FIBRE_T &pnts) {
    ofstream fileout(filename.c_str());
    int i,j;
    for (i=0;i<pnts.size();i++) {
      for (j=0;j<pnts[i].size()-1;j++) 
	fileout <<pnts[i][j]<<",";
      fileout <<pnts[i][j]<<std::endl;
    }
    fileout.close();
  }
}

