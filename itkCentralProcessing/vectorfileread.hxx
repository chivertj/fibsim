#ifndef VECTORFILEREAD_HEADER
#define VECTORFILEREAD_HEADER

#include<vector>
#include<string>
#include<fstream>
#include<iostream>
#include<sstream>

using namespace std;



namespace FRC {
  typedef float ELEM_T;
  typedef vector<ELEM_T> PNT_T;
  typedef vector<PNT_T> FIBRE_T;
  typedef vector<FIBRE_T> VOL_T;
  VOL_T vectorfileread(const string &filename);
  void vectorfilewrite(const string &filename, const VOL_T &fibres);
  void vectorfilewrite(const string &filename, const FIBRE_T &pnts);
}

#endif
