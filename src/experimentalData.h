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
#include <limits>

#include "ktlMoleculeRandom.h" // brings in point + ktlMolecule

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
  double C2{};
  explicit experimentalData(const char* scatterFile);

  // interpolation helpers
  void calculate_spline_coefficients(const std::vector<double>& x, const std::vector<double>& y,
                                     std::vector<double>& A, std::vector<double>& B,
                                     std::vector<double>& C, std::vector<double>& D);
  double evaluate_spline(double x_val, const std::vector<double>& x,
                         const std::vector<double>& A, const std::vector<double>& B,
                         const std::vector<double>& C, const std::vector<double>& D);

  // data preparation
  bool binDataCheck(double &dMax,double &qmin,double &qmax);
  void subsetScatteringData(std::vector<double>& A,std::vector<double>& B,std::vector<double>& C,
                            double kmin, double kmax,
                            std::vector<double>& A_selected,std::vector<double>& B_selected,std::vector<double>& C_selected);
  int  setPhases(double &dMax,double &kmin,double &kmax);

  // geometry + placement
  std::vector<point> calculate_geometric_normals(const std::vector<point>& ca_coords) const;
  std::vector<point> place_side_chains(const std::vector<point>& ca_coords,
                                       const std::vector<point>& geometric_vectors,
                                       const std::vector<char>& residue_names) const;
  std::vector<std::vector<double> > calculate_distances(const std::vector<point>& coordinates,int& molIndex);
  ScatteringCenters process_structure(const std::vector<point>& ca_coords,
                                      const std::vector<char>& residue_names,int& molIndex);
  std::vector<double> calculate_saxs_implicit(ScatteringCenters& centers);

  // flatteners
  std::vector<point> flatten_coords(const std::vector<std::vector<point>>& coords) const;
  std::vector<char>  flatten_residueNames(const std::vector<std::vector<std::string>>& aminoList) const;

  // chi^2
  double calculateChiSquared(std::vector<ktlMolecule> &mol,double &qmin,double &qmax,std::vector<std::vector<double> > &mixtureVals);
  double calculateChiSquaredUpdate(ktlMolecule& molNew,int& k,double &qmin,double &qmax,std::vector<std::vector<double> > &mixtureVals);
  double calculateChiSquared_Weighted(std::vector<ktlMolecule> &mol,double &qmin,double &qmax,std::vector<std::vector<double> > &mixtureVals);
  double calculateChiSquaredUpdate_Weighted(ktlMolecule& molNew,int& k,double &qmin,double &qmax,std::vector<std::vector<double> > &mixtureVals);

  // IO
  void writeScatteringToFile(std::vector<std::vector<double> > &mixtureVals,const char* filename);
  void writeScatteringToFile_ChiSq(std::vector<std::vector<double> > &mixtureVals,const char* filename);

  // convenience (kept for compatibility)
  std::vector<double> calculate_intensity_at_experimental_q(std::vector<double>& I_mod);
  std::vector<double> calculate_intensity_at_experimental_q(std::vector<double>& q_mod,std::vector<double>& I_mod,std::vector<double>& expQRange);

  // legacy wrapper, if referenced elsewhere
  double calculateChiSquaredTest(std::vector<ktlMolecule> &mol,double &qmin,double &qmax,std::vector<std::vector<double> > &mixtureVals){
      return calculateChiSquared(mol,qmin,qmax,mixtureVals);
  }

private:
  std::map<char, std::vector<double> > form_factors;
  std::vector<std::pair<double,double> > exprDat;
  std::vector<double> qvals;
  double absKmin{};
  double absKmax{};
  double kMin{};
  double kMax{};
  int noDistBins{10000};

  // per-run containers (make sure to clear/resize appropriately)
  std::vector<double> maxDist;
  std::vector<std::pair<double,double> > distBins;
  std::vector<double> experimentalIntensity;
  std::vector<std::vector<double> > scatVec;

  std::vector<ScatteringCenters> scs;
  std::vector<std::vector<double> > Ivec;
  std::vector<std::vector<double> > IvecOnData;
  std::vector<double> exprIInterp;

  // experimental data
  std::vector<double> exprQset;
  std::vector<double> exprIset;
  std::vector<double> exprEset;
  std::vector<double> exprQSubset;
  std::vector<double> exprISubset;
  std::vector<double> exprESubset;

  // model q grid (0.00 .. 0.20 by 0.01)
  std::vector<double> q;
};

#endif // EXP_DAT
