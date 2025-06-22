#ifndef EXP_DAT
#define EXP_DAT

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <map>
#include <numeric> 
#include "ktlMoleculeRandom.h"

struct ScatteringCenters {
    std::vector<point> coordinates;
    std::vector<char> types;
    std::vector<std::vector<double>> distances;
};

struct ExperimentalData {
   std::vector<double> q;
   std::vector<double> I;
   std::vector<double> I_err;
};

class experimentalData{
 public:
  // best fitting hydration parameter
  double C2;
  experimentalData(const char* scatterFile);
  void calculate_spline_coefficients(const std::vector<double>& x, const std::vector<double>& y,std::vector<double>& A, std::vector<double>& B,std::vector<double>& C, std::vector<double>& D);
  double evaluate_spline(double x_val, const std::vector<double>& x,const std::vector<double>& A, const std::vector<double>& B,const std::vector<double>& C, const std::vector<double>& D);
  bool binDataCheck(double &dMax,double &qmin,double &qmax);
  void subsetScatteringData(std::vector<double>& A,std::vector<double>& B,std::vector<double>& C,double kmin, double kmax,std::vector<double>& A_selected,std::vector<double>& B_selected,std::vector<double>& C_selected);
  int setPhases(double &dMax,double &kmin,double &kmax);
  std::vector<point> calculate_geometric_normals(std::vector<point>& ca_coords);
  std::vector<point> place_side_chains( std::vector<point>& ca_coords,std::vector<point>& geometric_vectors,std::vector<char>& residue_names);
  std::vector<std::vector<double> > calculate_distances( std::vector<point>& coordinates,int& molIndex); 
  ScatteringCenters process_structure(std::vector<point>& ca_coords,std::vector<char>& residue_names,int& molIndex);
  std::vector<double> calculate_saxs_implicit( ScatteringCenters& centers );
  std::vector<point> flatten_coords( std::vector<std::vector<point> >& coords);
  std::vector<char> flatten_residueNames( std::vector<std::vector<std::string> >& aminoList);
  double calculateChiSquared(std::vector<ktlMolecule> &mol,double &qmin,double &qmax,std::vector<std::vector<double> > &mixtureVals);
  double calculateChiSquaredTest(std::vector<ktlMolecule> &mol,double &qmin,double &qmax,std::vector<std::vector<double> > &mixtureVals);
  std::vector<double> calculate_intensity_at_experimental_q(std::vector<double>& I_mod);
  std::vector<double> calculate_intensity_at_experimental_q(std::vector<double>& q_mod,std::vector<double>& I_mod,std::vector<double>& expQRange);
  double calculateChiSquaredUpdate(ktlMolecule& molNew,int& k,double &qmin,double &qmax,std::vector<std::vector<double> > &mixtureVals);
  double calculateChiSquared_Weighted(std::vector<ktlMolecule> &mol,double &qmin,double &qmax,std::vector<std::vector<double> > &mixtureVals);
  double calculateChiSquaredUpdate_Weighted(ktlMolecule& molNew,int& k,double &qmin,double &qmax,std::vector<std::vector<double> > &mixtureVals);
  void writeScatteringToFile(std::vector<std::vector<double> > &mixtureVals,const char* filename);
  void writeScatteringToFile_ChiSq(std::vector<std::vector<double> > &mixtureVals,const char* filename);
private:
  std::map< char, std::vector<double> > form_factors;
  std::vector<std::pair<double,double> > exprDat;
  std::vector<double> qvals;
  double absKmin;
  double absKmax;
  double kMin;
  double kMax;
  int noDistBins;
  std::vector<double> maxDist;
  std::vector<std::pair<double,double> > distBins;
  std::vector<double> experimentalIntensity;
  std::vector<std::vector<double> > scatVec;
  std::vector<ScatteringCenters> scs;
  std::vector<std::vector<double> > Ivec;
  std::vector<std::vector<double> > IvecOnData;
  std::vector<double> exprIInterp;
  std::vector<double> q;
  std::vector<double> exprQset;
  std::vector<double> exprIset;
  std::vector<double> exprEset;
  std::vector<double> exprQSubset;
  std::vector<double> exprISubset;
  std::vector<double> exprESubset;
};

#endif
