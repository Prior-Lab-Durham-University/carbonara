#include "experimentalData.h"

static inline double safe_log(double x){
    const double eps = 1e-300;
    return std::log(x < eps ? eps : x);
}

experimentalData::experimentalData(const char* scatterFile){
  noDistBins = 10000;

  // read in the scattering file
  std::ifstream scatterdat(scatterFile);
  std::string sctline;
  double kval{}, Ival{}, Eval{};

  if(scatterdat.is_open()){
     while(std::getline(scatterdat,sctline)){
        if (sctline.empty()) continue;
        std::stringstream ss(sctline);
        if (!(ss >> kval >> Ival >> Eval)) continue;

        scatVec.push_back({kval, Ival, Eval});
        exprQset.push_back(kval);
        exprIset.push_back(Ival);
        exprEset.push_back(Eval);
     }
  }else{
    std::cout<<"no scattering file found\n";
  }

  if (!scatVec.empty()){
      absKmin = scatVec.front()[0];
      absKmax = scatVec.back()[0];
  } else {
      absKmin = 0.0;
      absKmax = 0.0;
  }

  // pre-tabulated form factors
  form_factors['B'] = {5.325233, 5.321616, 5.321073, 5.330379, 5.367118, 5.418221, 5.548391, 5.694481, 5.802268, 5.880531, 5.900372, 5.865756, 5.778460, 5.611933, 5.353311, 5.065751, 4.844303, 4.655859, 4.507514, 4.382296, 4.200807};
  form_factors['A'] = {6.860783, 6.884742, 6.936093, 6.987066, 7.002729, 7.001643, 7.011408, 7.059012, 7.136653, 7.204735, 7.278987, 7.358991, 7.468390, 7.587982, 7.698403, 7.792850, 7.835658, 7.747937, 7.795444, 8.008073, 8.268595};
  form_factors['R'] = {15.586594, 15.603896, 15.644747, 15.706051, 15.731217, 15.786658, 15.751397, 15.723773, 15.667085, 15.568371, 15.501597, 15.510941, 15.660628, 15.938443, 16.391161, 17.138668, 18.150400, 19.201813, 19.957359, 20.556793, 21.194653};
  form_factors['N'] = {14.042562, 14.039913, 14.030098, 14.046182, 14.072546, 14.079157, 13.947478, 13.831409, 13.740083, 13.667616, 13.673161, 13.796733, 14.041808, 14.428295, 14.873251, 15.052513, 14.851926, 14.920277, 15.730271, 17.125431, 18.337770};
  form_factors['D'] = {18.276451, 18.266953, 18.238106, 18.188072, 18.080282, 17.894873, 17.509510, 17.006626, 16.529716, 16.041245, 15.567556, 15.144760, 14.858373, 14.703475, 14.630036, 14.854190, 15.309014, 15.771427, 16.001890, 16.368433, 16.704248};
  form_factors['C'] = {10.351711, 10.346272, 10.323308, 10.293743, 10.243892, 10.212900, 10.157944, 10.163668, 10.267132, 10.432532, 10.664062, 10.914792, 11.118233, 11.204067, 11.197291, 11.052269, 10.771274, 10.325391, 9.898160, 9.588455, 9.458176};
  form_factors['Q'] = {13.376416, 13.384198, 13.387307, 13.362279, 13.276151, 13.152504, 12.912460, 12.665339, 12.532801, 12.473034, 12.546851, 12.777804, 13.243794, 13.995014, 15.000033, 15.920972, 16.500341, 17.029819, 17.050764, 16.716652, 16.530117};
  form_factors['E'] = {21.686304, 21.668835, 21.592291, 21.449541, 21.200043, 20.835613, 20.150393, 19.448851, 18.805557, 18.179039, 17.611521, 17.140610, 16.856714, 16.806103, 16.947565, 17.285227, 17.668753, 18.353695, 19.242361, 19.913992, 20.734098};
  form_factors['G'] = {10.030974, 10.004344, 9.943201, 9.875965, 9.802672, 9.715462, 9.611974, 9.495539, 9.422803, 9.374266, 9.338644, 9.343287, 9.433417, 9.612279, 9.818403, 10.084248, 10.328732, 10.578753, 10.550336, 10.467358, 10.294571};
  form_factors['H'] = {13.493807, 13.499746, 13.497142, 13.453699, 13.357458, 13.200852, 12.953443, 12.703934, 12.458268, 12.168073, 11.901849, 11.685180, 11.584588, 11.724971, 12.761358, 13.564467, 14.011192, 14.246178, 14.511696, 15.254260};
  form_factors['I'] = {0.027768, 0.028121, 0.028573, 0.023142, 0.007413, 0.003274, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.010318, 0.029916, 0.021361, 0.011234, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  form_factors['L'] = {0.001379, 0.001181, 0.000961, 0.000118, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.012047, 0.007686, 0.000000, 0.000000, 0.001314, 0.006526, 0.009532, 0.008636};
  form_factors['K'] = {9.349114, 9.359500, 9.370693, 9.372478, 9.334329, 9.252807, 9.074130, 8.906091, 8.755677, 8.595029, 8.477226, 8.410988, 8.422702, 8.465599, 8.387974, 8.202576, 7.873840, 7.777987, 8.014545, 8.801325, 10.069035};
  form_factors['M'] = {5.460539, 5.465856, 5.488052, 5.560043, 5.633809, 5.795099, 6.002742, 6.323314, 6.623469, 6.956409, 7.325475, 7.701069, 8.105552, 8.563451, 9.163991, 9.708839, 9.891646, 9.833041, 9.734925, 9.476597, 9.286131};
  form_factors['F'] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.001274, 0.020930, 0.086051, 0.251851, 0.409124, 0.562665, 0.686578, 0.832251, 1.048074, 1.359721, 1.720046, 2.006726, 2.087998, 2.053319, 2.235860, 2.837866};
  form_factors['P'] = {4.716053, 4.711367, 4.685131, 4.642177, 4.559615, 4.457505, 4.098892, 3.705757, 3.372101, 3.053588, 2.797038, 2.656890, 2.760468, 3.205317, 3.912106, 4.736746, 5.429208, 6.331475, 7.254241, 7.861680, 7.935802};
  form_factors['S'] = {8.985313, 9.003293, 9.028045, 9.040872, 8.994022, 8.931241, 8.708052, 8.462237, 8.256139, 8.054568, 7.901446, 7.829412, 7.892257, 8.112808, 8.520312, 9.034179, 9.444173, 9.954925, 10.321150, 9.990824, 9.688856};
  form_factors['T'] = {6.628838, 6.637546, 6.641963, 6.629081, 6.596369, 6.589334, 6.517858, 6.426005, 6.370910, 6.309315, 6.288209, 6.310956, 6.417308, 6.626843, 6.961301, 7.252404, 7.440239, 7.624983, 7.722063, 7.757763, 8.270276};
  form_factors['W'] = {2.888272, 2.949734, 3.104849, 3.325679, 3.557907, 3.812995, 4.198963, 4.642990, 5.063122, 5.487359, 5.964144, 6.477489, 7.032223, 7.543501, 7.751148, 7.804211, 7.761533, 8.023214, 8.654520, 9.687876, 10.764302};
  form_factors['Y'] = {5.100512, 5.130500, 5.228212, 5.390429, 5.591930, 5.801286, 5.978783, 6.084779, 6.201869, 6.252052, 6.266182, 6.257378, 6.282480, 6.352914, 6.603123, 6.989089, 7.407878, 7.660326, 7.929643, 8.116093, 8.308664};
  form_factors['V'] = {0.028584, 0.028193, 0.026019, 0.023774, 0.010495, 0.002455, 0.000000, 0.000000, 0.006006, 0.012980, 0.016876, 0.014590, 0.028809, 0.029842, 0.008057, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};

  // model q grid (0.00 .. 0.20 step 0.01)
  q.clear();
  for (int i = 0; i <= 19; i++) { q.push_back(0.2/19.0 * i); }
}

void experimentalData::calculate_spline_coefficients(const std::vector<double>& x, const std::vector<double>& y,
                                 std::vector<double>& A, std::vector<double>& B,
                                 std::vector<double>& C, std::vector<double>& D) {
    const size_t n = x.size();
    std::vector<double> h(n - 1), alpha(n - 1), l(n), mu(n), z(n);

    for (size_t i = 0; i < n - 1; i++) h[i] = x[i + 1] - x[i];
    for (size_t i = 1; i < n - 1; i++) alpha[i] = (3.0/h[i])*(y[i+1] - y[i]) - (3.0/h[i-1])*(y[i] - y[i-1]);

    l[0] = 1.0; mu[0] = 0.0; z[0] = 0.0;
    for (size_t i = 1; i < n - 1; i++) {
        l[i] = 2.0 * (x[i+1] - x[i-1]) - h[i-1]*mu[i-1];
        mu[i] = h[i] / l[i];
        z[i]  = (alpha[i] - h[i-1]*z[i-1]) / l[i];
    }
    l[n-1] = 1.0; z[n-1] = 0.0;

    std::vector<double> c(n, 0.0);
    for (int j = static_cast<int>(n) - 2; j >= 0; --j) {
        c[static_cast<size_t>(j)] = z[static_cast<size_t>(j)] - mu[static_cast<size_t>(j)]*c[static_cast<size_t>(j+1)];
    }

    A.resize(n-1); B.resize(n-1); C.resize(n-1); D.resize(n-1);
    for (size_t i = 0; i < n - 1; i++) {
        A[i] = (c[i+1] - c[i])/(3.0*h[i]);
        B[i] = c[i];
        C[i] = (y[i+1] - y[i])/h[i] - h[i]*(c[i+1] + 2.0*c[i])/3.0;
        D[i] = y[i];
    }
}

double experimentalData::evaluate_spline(double x_val, const std::vector<double>& x,
                      const std::vector<double>& A, const std::vector<double>& B,
                      const std::vector<double>& C, const std::vector<double>& D) {
    size_t i = 0;
    while (i < x.size() - 2 && x[i + 1] <= x_val) i++;
    const double dx = x_val - x[i];
    return ((A[i]*dx + B[i])*dx + C[i])*dx + D[i];
}

bool experimentalData::binDataCheck(double &dMax,double &qmin,double &qmax){
  const double dq = (qmax-qmin)/double(noDistBins);
  bool goodsplit = true;
  size_t k = 0;
  for(int j=1;j<=noDistBins;j++){
    const double qminBin = qmin + (j-1)*dq;
    const double qmaxBin = qmin + j*dq;
    std::vector<double> intensities;

    while(k < scatVec.size() && scatVec[k][0] < qminBin) { k++; }
    while(k < scatVec.size() && scatVec[k][0] >= qminBin && scatVec[k][0] <= qmaxBin){
      intensities.push_back(scatVec[k][1]);
      k++;
    }
    if(intensities.size() < 2){ goodsplit = false; }
  }
  return goodsplit;
}

void experimentalData::subsetScatteringData(std::vector<double>& A,std::vector<double>& B,std::vector<double>& C,
                      double kmin, double kmax,
                      std::vector<double>& A_selected,
                      std::vector<double>& B_selected,
                      std::vector<double>& C_selected)
{
    A_selected.clear(); B_selected.clear(); C_selected.clear();
    for (size_t i = 0; i < A.size(); ++i) {
        if (A[i] >= kmin && A[i] <= kmax) {
            A_selected.push_back(A[i]);
            B_selected.push_back(B[i]);
            C_selected.push_back(C[i]);
        }
    }
}

int experimentalData::setPhases(double &dMax,double &qmin,double &qmax){
  int noDistBinsTemp = int(2.1*std::ceil((qmax-qmin)*dMax/3.14159265359));

  if(noDistBins != noDistBinsTemp){
    noDistBins = noDistBinsTemp;

    bool binsOk=false;
    while(!binsOk && noDistBins>2){
      binsOk = binDataCheck(dMax,qmin,qmax);
      if(!binsOk){ noDistBins--; }
    }

    const double dr = (dMax) / double(noDistBins);
    const double dq = (qmax - qmin) / double(noDistBins);
    experimentalIntensity.clear(); distBins.clear(); qvals.clear();

    size_t k=0;
    for(int i=1;i<=noDistBins;i++){
      std::pair<double,double> rDistSets{(i-1)*dr, i*dr};
      distBins.push_back(rDistSets);
      for(int j=1;j<=noDistBins;j++){
        const double q = qmin + (j-0.5)*dq;
        if (i==1){
          qvals.push_back(q);
          const double qminBin = qmin + (j-1)*dq;
          const double qmaxBin = qmin + j*dq;
          std::vector<double> intensities;

          while(k < scatVec.size() && scatVec[k][0] < qminBin) { k++; }
          while(k < scatVec.size() && scatVec[k][0] >= qminBin && scatVec[k][0] <= qmaxBin){
            intensities.push_back(scatVec[k][1]);
            k++;
          }
          std::sort(intensities.begin(), intensities.end());
          if (!intensities.empty()){
            int size = static_cast<int>(intensities.size());
            int medPt = static_cast<int>(std::round(double(size)/2.0));
            medPt = std::max(1, std::min(medPt, size));
            experimentalIntensity.push_back(intensities[medPt-1]);
          } else {
            experimentalIntensity.push_back(0.0);
          }
        }
      }
    }
  }

  subsetScatteringData(exprQset,exprIset,exprEset,qmin,qmax,exprQSubset,exprISubset,exprESubset);
  return noDistBins;
}

std::vector<point> experimentalData::calculate_geometric_normals(const std::vector<point>& ca_coords) const {
    std::vector<point> ca_vectors;
    std::vector<point> normals;

    for (size_t i = 1; i < ca_coords.size(); i++) {
        ca_vectors.push_back(ca_coords[i] - ca_coords[i-1]);
    }
    for (size_t i = 1; i < ca_vectors.size(); i++) {
        point norm = ca_vectors[i] - ca_vectors[i-1];
        norm = norm * (-1.0) / norm.length();
        normals.push_back(norm);
    }
    if (!ca_vectors.empty()){
        point first_norm = ca_vectors[0] * (-1.0) / ca_vectors[0].length();
        point final_norm = ca_vectors.back() / ca_vectors.back().length();
        normals.insert( normals.begin(), first_norm );
        normals.push_back( final_norm );
    }
    return normals;
}

std::vector<point> experimentalData::place_side_chains(const std::vector<point>& ca_coords,
                                      const std::vector<point>& geometric_vectors,
                                      const std::vector<char>& residue_names) const {
    std::map<char, double> residue_distances = {
        {'R', 4.2662}, {'N', 2.5349}, {'D', 2.5558}, {'C', 2.3839},
        {'Q', 3.1861}, {'E', 3.2541}, {'H', 3.1861}, {'I', 2.3115},
        {'L', 2.6183}, {'K', 3.6349}, {'M', 3.1912}, {'F', 3.4033},
        {'P', 1.8773}, {'U', 1.5419}, {'S', 1.9661}, {'T', 1.9533},
        {'W', 3.8916}, {'Y', 3.8807}, {'V', 1.9555}, {'G', 0.0}, {'A', 0.0}
    };

    std::vector<point> side_chain_positions;
    side_chain_positions.reserve(ca_coords.size());

    for (size_t i = 0; i < ca_coords.size(); i++) {
        double distance = 0.0;
        auto it = residue_distances.find(residue_names[i]);
        if (it != residue_distances.end()) distance = it->second;
        side_chain_positions.push_back(ca_coords[i] + geometric_vectors[i] * distance);
    }
    return side_chain_positions;
}

std::vector<std::vector<double> > experimentalData::calculate_distances(const std::vector<point>& coordinates, int& molIndex) {
    const size_t n = coordinates.size();
    std::vector< std::vector<double> > distances(n, std::vector<double>(n, 0.0));
    double min_dist = std::numeric_limits<double>::max();

    for (size_t i = 0; i < n; i++){
        for (size_t j = i; j < n; j++){
            if (i==j) {
                distances[i][j] = 0.0;
            } else {
                double dist = coordinates[i].eDist(coordinates[j]);
                if (molIndex >= 0 && static_cast<size_t>(molIndex) < maxDist.size()) {
		  if (dist > maxDist[molIndex]) maxDist[molIndex] = dist;
		}
                distances[i][j] = dist;
                distances[j][i] = dist;
            }
        }
    }
    (void)min_dist; // not used, but keep structure
    return distances;
}

ScatteringCenters experimentalData::process_structure(const std::vector<point>& ca_coords,
						      const std::vector<char>& residue_names,int& molIndex){
    std::vector<point> geometric_normals = calculate_geometric_normals(ca_coords);
    std::vector<point> side_chain_positions = place_side_chains(ca_coords, geometric_normals, residue_names);

    ScatteringCenters centers;

    for (size_t i=0; i < residue_names.size(); i++) {
        if (residue_names[i] != 'A' && residue_names[i] != 'G') {
            centers.coordinates.push_back(ca_coords[i]);
            centers.types.push_back('B');
            centers.coordinates.push_back(side_chain_positions[i]);
            centers.types.push_back(residue_names[i]);
        }
    }
    for (size_t i=0; i < residue_names.size(); i++) {
        if (residue_names[i] == 'A' || residue_names[i] == 'G') {
            centers.coordinates.push_back(ca_coords[i]);
            centers.types.push_back(residue_names[i]);
        }
    }
    centers.distances = calculate_distances(centers.coordinates,molIndex);
    return centers;
}

std::vector<double> experimentalData::calculate_saxs_implicit(ScatteringCenters& centers) {
    std::vector<double> I(q.size(), 0.0);
    std::vector< std::vector<double> > ff_matrix;
    ff_matrix.reserve(centers.types.size());
    for (char type : centers.types) {
        ff_matrix.push_back( form_factors.at(type) );
    }
    for (size_t q_idx = 0; q_idx < q.size(); q_idx++) {
        double qv = q[q_idx];
        for (size_t i = 0; i < centers.coordinates.size(); i++) {
            for (size_t j = 0; j < centers.coordinates.size(); j++) {
                double r_ij = centers.distances[i][j];
                double ff_product = ff_matrix[i][q_idx] * ff_matrix[j][q_idx];
                double sinc = (qv < 1e-12) ? 1.0 : (r_ij == 0.0 ? 1.0 : std::sin(qv * r_ij) / (qv * r_ij));
                I[q_idx] += ff_product * sinc;
            }
        }
    }
    return I;
}

std::vector<point> experimentalData::flatten_coords(const std::vector<std::vector<point>>& coords) const {
   std::vector<point> flat;
   size_t total = 0;
   for (const auto& section : coords) total += section.size();
   flat.reserve(total);
   for (const auto& section : coords) {
       flat.insert(flat.end(), section.begin(), section.end());
   }
   return flat;
}

static std::vector<char> convertStringsToChars(const std::vector<std::string>& strings) {
    std::vector<char> out;
    size_t total = 0;
    for (const auto& s : strings) total += s.size();
    out.reserve(total);
    for (const auto& str : strings) {
        out.insert(out.end(), str.begin(), str.end());
    }
    return out;
}

std::vector<char> experimentalData::flatten_residueNames(const std::vector<std::vector<std::string>>& aminoList) const {
   std::vector<char> flat;
   for (const auto& section : aminoList) {
     std::vector<char> charSec = convertStringsToChars(section);
     flat.insert(flat.end(), charSec.begin(), charSec.end());
   }
   return flat;
}

std::vector<double> experimentalData::calculate_intensity_at_experimental_q(std::vector<double>& I_mod) {
    size_t N = q.size() - 1;
    std::vector<double> A(N), B(N), C(N), D(N);
    calculate_spline_coefficients(q, I_mod, A, B, C, D);

    std::vector<double> I_interp;
    I_interp.reserve(qvals.size());
    for(size_t i=0;i<qvals.size();i++){
      I_interp.push_back(evaluate_spline(qvals[i], q, A, B, C, D));
    }
    return I_interp;
}

std::vector<double> experimentalData::calculate_intensity_at_experimental_q(std::vector<double>& q_mod,std::vector<double>& I_mod,std::vector<double>& expQRange) {
    size_t N = q_mod.size() - 1;
    std::vector<double> A(N), B(N), C(N), D(N);
    calculate_spline_coefficients(q_mod, I_mod, A, B, C, D);

    std::vector<double> I_interp;
    I_interp.reserve(expQRange.size());
    for(size_t i = 0; i<expQRange.size(); i++) {
      I_interp.push_back(evaluate_spline(expQRange[i], q_mod, A, B, C, D));
    }
    return I_interp;
}

double experimentalData::calculateChiSquared(std::vector<ktlMolecule> &mol,double &qmin,double &qmax,std::vector<std::vector<double> > &mixtureVals){
   kMin = qmin; kMax = qmax;

   // reset per-run containers
   maxDist.assign(mol.size(), 0.0);
   scs.clear(); scs.resize(mol.size());
   Ivec.clear(); Ivec.resize(mol.size());
   IvecOnData.clear(); IvecOnData.resize(mol.size());

   // build centers per molecule
   for(size_t i=0;i<mol.size();i++){
     std::vector<std::vector<point> > coords = mol[i].getCoordinates();
     std::vector<point> chainCoordinates = flatten_coords(coords);
     std::vector<std::vector<std::string> > aminoList = mol[i].getAminoList();
     std::vector<char> chainResidues = flatten_residueNames(aminoList);
     int idx = static_cast<int>(i);
     scs[i] = process_structure(chainCoordinates,chainResidues,idx);
   }

   double dMax = 0.0;
   if (!maxDist.empty()) dMax = *std::max_element(begin(maxDist), end(maxDist));
   (void)setPhases(dMax,qmin,qmax);

   for(size_t i=0;i<mol.size();i++){
     Ivec[i] = calculate_saxs_implicit(scs[i]);
   }

   for(size_t i=0;i<mol.size();i++){
     IvecOnData[i] = calculate_intensity_at_experimental_q(q,Ivec[i],exprQSubset);
   }

   exprIInterp = calculate_intensity_at_experimental_q(qvals,experimentalIntensity,exprQSubset);

   double best = std::numeric_limits<double>::infinity();

   for(size_t i = 0; i<mixtureVals.size(); i++){
     std::vector<double> Icomb(IvecOnData[0].size(),0.0);
     for(size_t j = 0; j<mixtureVals[i].size(); j++){
       for(size_t k = 0; k<IvecOnData[j].size(); k++){
         Icomb[k] += IvecOnData[j][k]*mixtureVals[i][j]; // FIX: accumulate
       }
     }

     double logDifMean=0.0;
     int noMean=0;
     std::vector<double> logdifs;
     logdifs.reserve(Icomb.size());

     for(size_t l=0;l<Icomb.size();l++){
       double logScatDif= safe_log(Icomb[l]) - safe_log(exprIInterp[l]);
       logdifs.push_back(logScatDif);
       if(exprQSubset[l] < kMin + 0.01){
         logDifMean += logScatDif;
         noMean++;
       }
     }
     if (noMean>0) logDifMean /= double(noMean);

     double predTemp =0.0;
     for(size_t l=0;l<Icomb.size();l++){
       double scatDif = logdifs[l] - logDifMean;
       predTemp += scatDif*scatDif;
     }
     best = std::min(best, predTemp);
   }
   return best / std::max<size_t>(1, exprQSubset.size()-1);
}

double experimentalData::calculateChiSquaredUpdate(ktlMolecule& molNew,int& k,double &qmin,double &qmax,std::vector<std::vector<double> > &mixtureVals){
   kMin = qmin; kMax = qmax;

   if (static_cast<size_t>(k) >= scs.size()){
       scs.resize(k+1);
       Ivec.resize(k+1);
       IvecOnData.resize(k+1);
       maxDist.resize(k+1, 0.0);
   }
   maxDist[k] = 0.0;

   std::vector<std::vector<point> > coords = molNew.getCoordinates();
   std::vector<point> chainCoordinates = flatten_coords(coords);
   std::vector<std::vector<std::string> > aminoList = molNew.getAminoList();
   std::vector<char> chainResidues = flatten_residueNames(aminoList);
   scs[k] = process_structure(chainCoordinates,chainResidues,k);

   double dMax = 0.0;
   if (!maxDist.empty()) dMax = *std::max_element(begin(maxDist), end(maxDist));
   (void)setPhases(dMax,qmin,qmax);

   Ivec[k] = calculate_saxs_implicit(scs[k]);
   IvecOnData[k] = calculate_intensity_at_experimental_q(q,Ivec[k],exprQSubset);

   exprIInterp = calculate_intensity_at_experimental_q(qvals,experimentalIntensity,exprQSubset);

   double best = std::numeric_limits<double>::infinity();

   for(size_t i = 0; i<mixtureVals.size(); i++){
     std::vector<double> Icomb(IvecOnData[0].size(),0.0);
     for(size_t j = 0; j<mixtureVals[i].size(); j++){
       for(size_t kk = 0; kk<IvecOnData[j].size(); kk++){
         Icomb[kk] += IvecOnData[j][kk]*mixtureVals[i][j];
       }
     }

     double logDifMean=0.0;
     int noMean=0;
     std::vector<double> logdifs;
     logdifs.reserve(Icomb.size());
     for(size_t l=0;l<Icomb.size();l++){
       double logScatDif= safe_log(Icomb[l]) - safe_log(exprIInterp[l]);
       logdifs.push_back(logScatDif);
       if(exprQSubset[l] < kMin + 0.01){
         logDifMean += logScatDif;
         noMean++;
       }
     }
     if (noMean>0) logDifMean /= double(noMean);

     double predTemp =0.0;
     for(size_t l=0;l<Icomb.size();l++){
       double scatDif = logdifs[l] - logDifMean;
       predTemp += scatDif*scatDif;
     }
     best = std::min(best, predTemp);
   }
   return best / std::max<size_t>(1, exprQSubset.size()-1);
}

double experimentalData::calculateChiSquared_Weighted(std::vector<ktlMolecule> &mol,double &qmin,double &qmax,std::vector<std::vector<double> > &mixtureVals){
   kMin = qmin; kMax = qmax;

   maxDist.assign(mol.size(), 0.0);
   scs.clear(); scs.resize(mol.size());
   Ivec.clear(); Ivec.resize(mol.size());
   IvecOnData.clear(); IvecOnData.resize(mol.size());

   for(size_t i=0;i<mol.size();i++){
     std::vector<std::vector<point> > coords = mol[i].getCoordinates();
     std::vector<point> chainCoordinates = flatten_coords(coords);
     std::vector<std::vector<std::string> > aminoList = mol[i].getAminoList();
     std::vector<char> chainResidues = flatten_residueNames(aminoList);
     int idx = static_cast<int>(i);
     scs[i] = process_structure(chainCoordinates,chainResidues,idx);
   }

   double dMax = 0.0;
   if (!maxDist.empty()) dMax = *std::max_element(begin(maxDist), end(maxDist));
   (void)setPhases(dMax,qmin,qmax);

   for(size_t i=0;i<mol.size();i++) Ivec[i] = calculate_saxs_implicit(scs[i]);
   for(size_t i=0;i<mol.size();i++) IvecOnData[i] = calculate_intensity_at_experimental_q(q,Ivec[i],exprQSubset);

   double best = std::numeric_limits<double>::infinity();

   for(size_t i = 0; i<mixtureVals.size(); i++){
     std::vector<double> Icomb(IvecOnData[0].size(),0.0);
     for(size_t j = 0; j<mixtureVals[i].size(); j++){
       for(size_t k = 0; k<IvecOnData[j].size(); k++){
         Icomb[k] += IvecOnData[j][k]*mixtureVals[i][j]; // FIX: accumulate
       }
     }

     double logDifMean=0.0;
     int noMean=0;
     std::vector<double> logdifs;
     logdifs.reserve(Icomb.size());

     for(size_t l=0;l<Icomb.size();l++){
       double logScatDif= safe_log(Icomb[l]) - safe_log(exprISubset[l]);
       logdifs.push_back(logScatDif);
       if(exprQSubset[l] < kMin + 0.075){
         logDifMean += logScatDif;
         noMean++;
       }
     }
     if (noMean>0) logDifMean /= double(noMean);

     // rescale
     for(size_t l=0;l<Icomb.size();l++){
       double logScatScaled = safe_log(Icomb[l]) - logDifMean;
       Icomb[l] = std::exp(logScatScaled);
     }

     // chi^2
     std::vector<double> Imodel = Icomb;
     double predTemp = 0.0;
     for(size_t l = 0; l<exprQSubset.size(); l++){
       double dif = Imodel[l] - exprISubset[l];
       predTemp += dif*dif/(exprESubset[l]*exprESubset[l] + 1e-300);
     }
     predTemp = std::abs(predTemp/std::max<size_t>(1, exprQSubset.size()-1) - 1.0);
     best = std::min(best, predTemp);
   }
   return best;
}

double experimentalData::calculateChiSquaredUpdate_Weighted(ktlMolecule& molNew,int& k,double &qmin,double &qmax,std::vector<std::vector<double> > &mixtureVals){
   kMin = qmin; kMax = qmax;

   if (static_cast<size_t>(k) >= scs.size()){
       scs.resize(k+1);
       Ivec.resize(k+1);
       IvecOnData.resize(k+1);
       maxDist.resize(k+1, 0.0);
   }
   maxDist[k] = 0.0;

   std::vector<std::vector<point> > coords=molNew.getCoordinates();
   std::vector<point> chainCoordinates = flatten_coords(coords);
   std::vector<std::vector<std::string> > aminoList = molNew.getAminoList();
   std::vector<char> chainResidues = flatten_residueNames(aminoList);
   scs[k] = process_structure(chainCoordinates,chainResidues,k);

   double dMax = 0.0;
   if (!maxDist.empty()) dMax = *std::max_element(begin(maxDist), end(maxDist));
   (void)setPhases(dMax,qmin,qmax);

   Ivec[k] = calculate_saxs_implicit(scs[k]);
   IvecOnData[k] = calculate_intensity_at_experimental_q(q,Ivec[k],exprQSubset);

   double best = std::numeric_limits<double>::infinity();

   for(size_t i = 0; i<mixtureVals.size(); i++){
     std::vector<double> Icomb(IvecOnData[0].size(),0.0);
     for(size_t j = 0; j<mixtureVals[i].size(); j++){
       for(size_t kk = 0; kk<IvecOnData[j].size(); kk++){
         Icomb[kk] += IvecOnData[j][kk]*mixtureVals[i][j];
       }
     }

     double logDifMean=0.0;
     int noMean=0;
     std::vector<double> logdifs;
     logdifs.reserve(Icomb.size());
     for(size_t l=0;l<Icomb.size();l++){
       double logScatDif= safe_log(Icomb[l]) - safe_log(exprISubset[l]);
       logdifs.push_back(logScatDif);
       if(exprQSubset[l] < kMin + 0.075){
         logDifMean += logScatDif;
         noMean++;
       }
     }
     if (noMean>0) logDifMean /= double(noMean);

     for(size_t l=0;l<Icomb.size();l++){
       double logScatScaled = safe_log(Icomb[l]) - logDifMean;
       Icomb[l] = std::exp(logScatScaled);
     }

     std::vector<double> Imodel = Icomb;
     double predTemp = 0.0;
     for(size_t l = 0; l<exprQSubset.size(); l++){
       double dif = Imodel[l] - exprISubset[l];
       predTemp += dif*dif/(exprESubset[l]*exprESubset[l] + 1e-300);
     }
     predTemp = std::abs(predTemp/std::max<size_t>(1, exprQSubset.size()-1) - 1.0);
     best = std::min(best, predTemp);
   }
   return best;
}

void experimentalData::writeScatteringToFile(std::vector<std::vector<double> > &mixtureVals,const char* filename){
  size_t best = 0;
  double logDifMeanFinal = 0.0;
  double bestScore = std::numeric_limits<double>::infinity();
  std::vector<double> IcombBest;

  exprIInterp = calculate_intensity_at_experimental_q(qvals,experimentalIntensity,exprQSubset);

  for(size_t i = 0; i<mixtureVals.size(); i++){
    std::vector<double> Icomb(IvecOnData[0].size(),0.0);
    for(size_t j = 0; j<mixtureVals[i].size(); j++){
      for(size_t k = 0; k<IvecOnData[j].size(); k++){
        Icomb[k] += IvecOnData[j][k]*mixtureVals[i][j]; // FIX accumulate
      }
    }

    double logDifMean=0.0;
    int noMean=0;
    std::vector<double> logdifs;
    logdifs.reserve(Icomb.size());
    for(size_t l=0;l<Icomb.size();l++){
      double logScatDif= safe_log(Icomb[l]) - safe_log(exprIInterp[l]);
      logdifs.push_back(logScatDif);
      if(exprQSubset[l] < kMin + 0.01){
        logDifMean += logScatDif;
        noMean++;
      }
    }
    if (noMean>0) logDifMean /= double(noMean);

    double predTemp = 0.0;
    for(size_t l=0;l<Icomb.size();l++){
      double scatDif = logdifs[l] - logDifMean;
      predTemp += scatDif*scatDif;
    }
    if(predTemp < bestScore){
        bestScore = predTemp;
        best = i;
        logDifMeanFinal = logDifMean;
        IcombBest = Icomb;
    }
  }

  std::ofstream myfile(filename);
  for(size_t i=0;i<IcombBest.size();i++){
    myfile << exprQSubset[i] << " "
           << IcombBest[i] << " "
           << (safe_log(IcombBest[i]) - logDifMeanFinal) << " "
           << safe_log(exprISubset[i]) << "\n";
  }
  for(size_t i=0;i<mixtureVals[best].size();i++){
    myfile << mixtureVals[best][i] << (i+1==mixtureVals[best].size() ? "\n" : " ");
  }
  myfile.close();
}

void experimentalData::writeScatteringToFile_ChiSq(std::vector<std::vector<double> > &mixtureVals,const char* filename){
  size_t best = 0;
  double bestScore = std::numeric_limits<double>::infinity();
  std::vector<double> IcombBest;

  for(size_t i = 0; i<mixtureVals.size(); i++){
    std::vector<double> Icomb(IvecOnData[0].size(),0.0);
    for(size_t j = 0; j<mixtureVals[i].size(); j++){
      for(size_t k = 0; k<IvecOnData[j].size(); k++){
        Icomb[k] += IvecOnData[j][k]*mixtureVals[i][j];
      }
    }

    double logDifMean=0.0;
    int noMean=0;
    std::vector<double> logdifs;
    logdifs.reserve(Icomb.size());
    for(size_t l=0;l<Icomb.size();l++){
      double logScatDif= safe_log(Icomb[l]) - safe_log(exprISubset[l]);
      logdifs.push_back(logScatDif);
      if(exprQSubset[l] < kMin + 0.075){
        logDifMean += logScatDif;
        noMean++;
      }
    }
    if (noMean>0) logDifMean /= double(noMean);

    for(size_t l=0;l<Icomb.size();l++){
      double logScatScaled = safe_log(Icomb[l]) - logDifMean;
      Icomb[l] = std::exp(logScatScaled);
    }

    std::vector<double> Imodel = Icomb;
    double predTemp = 0.0;
    for(size_t l = 0; l<exprQSubset.size(); l++){
      double dif = Imodel[l] - exprISubset[l];
      predTemp += dif*dif/(exprESubset[l]*exprESubset[l] + 1e-300);
    }
    predTemp = std::abs(predTemp/std::max<size_t>(1, exprQSubset.size()-1) - 1.0);
    if(predTemp < bestScore){
      bestScore = predTemp;
      IcombBest = Imodel;
      best = i;
    }
  }

  std::ofstream myfile(filename);
  for(size_t l = 0; l<exprQSubset.size(); l++){
    myfile << exprQSubset[l] << " " << IcombBest[l] << " "
           << safe_log(IcombBest[l]) << " " << safe_log(exprISubset[l]) << "\n";
  }
  for(size_t i=0;i<mixtureVals[best].size();i++){
    myfile << mixtureVals[best][i] << (i+1==mixtureVals[best].size() ? "\n" : " ");
  }
  myfile.close();
}


double experimentalData::calculateChiSquaredUpdate_FrozenGrid(const ktlMolecule& molNew, int k,
                                                              double &qmin, double &qmax,
                                                              const std::vector<std::vector<double>> &mixtureVals)
{
    // Keep these (used for the low-q centering threshold)
    kMin = qmin; kMax = qmax;

    // ---- Build local scattering centers WITHOUT touching ed state ----
    // Avoid updating maxDist[...] by passing a sentinel index that fails the check.
    int scratchMolIndex = -1;

    // Flatten inputs
    const auto coordsSections = molNew.getCoordinates();
    const auto chainCoords    = flatten_coords(coordsSections);
    const auto aminoList      = molNew.getAminoList();
    const auto residues       = flatten_residueNames(aminoList);

    // Local centers & intensity for the updated molecule only
    ScatteringCenters sc_local = process_structure(chainCoords, residues, scratchMolIndex);
    std::vector<double> I_local = calculate_saxs_implicit(sc_local);

    // Interpolate this molecule onto the CURRENT, COMMITTED experimental grid.
    // IMPORTANT: do NOT call setPhases here; just use exprQSubset as-is.
    std::vector<double> I_onData_k = calculate_intensity_at_experimental_q(q, I_local, exprQSubset);

    // We rely on the previously-committed experimental binning & centering.
    // If exprIInterp is empty (e.g. first call didn’t build it), build it once without rebinning.
    if (exprIInterp.empty()) {
        exprIInterp = calculate_intensity_at_experimental_q(qvals, experimentalIntensity, exprQSubset);
    }

    // Size the combination vector to the current grid (not IvecOnData[0], which may be stale-sized)
    const size_t N = exprQSubset.size();
    std::vector<double> Icomb(N, 0.0);

    // Combine mixtures using existing IvecOnData for j != k and local I_onData_k for j == k.
    for (size_t row = 0; row < mixtureVals.size(); ++row) {
        std::fill(Icomb.begin(), Icomb.end(), 0.0);

        for (size_t j = 0; j < mixtureVals[row].size(); ++j) {
            const double w = mixtureVals[row][j];

            const std::vector<double> *src = nullptr;
            if (static_cast<int>(j) == k) {
                src = &I_onData_k;
            } else if (j < IvecOnData.size()) {
                src = &IvecOnData[j]; // committed, on the same grid after last full rebuild
            }

            if (src && !src->empty()) {
                const size_t M = std::min(N, src->size()); // paranoia guard
                for (size_t t = 0; t < M; ++t) Icomb[t] += (*src)[t] * w;
            }
        }

        // Log-recentering vs experimental (same as your unweighted path)
        double logDifMean = 0.0; int noMean = 0;
        std::vector<double> logdifs; logdifs.reserve(N);
        for (size_t t = 0; t < N; ++t) {
            const double l = safe_log(Icomb[t]) - safe_log(exprIInterp[t]);
            logdifs.push_back(l);
            if (exprQSubset[t] < kMin + 0.01) { logDifMean += l; ++noMean; }
        }
        if (noMean > 0) logDifMean /= double(noMean);

        double predTemp = 0.0;
        for (size_t t = 0; t < N; ++t) {
            const double d = logdifs[t] - logDifMean;
            predTemp += d * d;
        }

        // Track best over mixture rows
        // We return normalized chi-like value, as your full path does.
        static thread_local double best = std::numeric_limits<double>::infinity();
        if (row == 0) best = std::numeric_limits<double>::infinity();
        if (predTemp < best) best = predTemp;
        if (row + 1 == mixtureVals.size()) {
            return best / std::max<size_t>(1, N - 1);
        }
    }
    // Shouldn’t reach here
    return std::numeric_limits<double>::infinity();
}

double experimentalData::calculateChiSquaredUpdate_Weighted_FrozenGrid(const ktlMolecule& molNew, int k,
                                                                       double &qmin, double &qmax,
                                                                       const std::vector<std::vector<double>> &mixtureVals)
{
    kMin = qmin; kMax = qmax;

    int scratchMolIndex = -1;

    const auto coordsSections = molNew.getCoordinates();
    const auto chainCoords    = flatten_coords(coordsSections);
    const auto aminoList      = molNew.getAminoList();
    const auto residues       = flatten_residueNames(aminoList);

    ScatteringCenters sc_local = process_structure(chainCoords, residues, scratchMolIndex);
    std::vector<double> I_local = calculate_saxs_implicit(sc_local);
    std::vector<double> I_onData_k = calculate_intensity_at_experimental_q(q, I_local, exprQSubset);

    const size_t N = exprQSubset.size();
    std::vector<double> Icomb(N, 0.0);

    // Weighted path uses exprISubset/exprESubset (already subset to current grid by last full rebuild)
    if (exprISubset.size() != N || exprESubset.size() != N) {
        // If this ever trips, your program previously changed the grid without a full rebuild.
        // Best effort: clamp to min size.
    }

    double best = std::numeric_limits<double>::infinity();

    for (size_t row = 0; row < mixtureVals.size(); ++row) {
        std::fill(Icomb.begin(), Icomb.end(), 0.0);

        for (size_t j = 0; j < mixtureVals[row].size(); ++j) {
            const double w = mixtureVals[row][j];
            const std::vector<double> *src = nullptr;
            if (static_cast<int>(j) == k) {
                src = &I_onData_k;
            } else if (j < IvecOnData.size()) {
                src = &IvecOnData[j];
            }
            if (src && !src->empty()) {
                const size_t M = std::min(N, src->size());
                for (size_t t = 0; t < M; ++t) Icomb[t] += (*src)[t] * w;
            }
        }

        // Low-q log centering → rescale model back to linear (as in your weighted code)
        double logDifMean = 0.0; int noMean = 0;
        for (size_t t = 0; t < N; ++t) {
            const double l = safe_log(Icomb[t]) - safe_log(exprISubset[t]);
            if (exprQSubset[t] < kMin + 0.075) { logDifMean += l; ++noMean; }
        }
        if (noMean > 0) logDifMean /= double(noMean);

        for (size_t t = 0; t < N; ++t) {
            Icomb[t] = std::exp(safe_log(Icomb[t]) - logDifMean);
        }

        // Classic chi^2 with errors
        double chi = 0.0;
        for (size_t t = 0; t < N; ++t) {
            const double dif = Icomb[t] - exprISubset[t];
            chi += dif * dif / (exprESubset[t] * exprESubset[t] + 1e-300);
        }
        chi = std::abs(chi / std::max<size_t>(1, N - 1) - 1.0);
        best = std::min(best, chi);
    }
    return best;
}

double experimentalData::scoreWeightedSandbox(const std::vector<ktlMolecule>& mol,
                                              double &kmin, double &kmax,
                                              const std::vector<std::vector<double>>& mixtureVals)
{
  EDStateSnapshot snap = snapshot();
  double chi = calculateChiSquared_Weighted(const_cast<std::vector<ktlMolecule>&>(mol),
                                            kmin, kmax,
                                            const_cast<std::vector<std::vector<double>>&>(mixtureVals));
  restore(snap);
  return chi;
}

double experimentalData::scoreWeightedCommit(const std::vector<ktlMolecule>& mol,
                                             double &kmin, double &kmax,
                                             const std::vector<std::vector<double>>& mixtureVals)
{
  // This is just the normal stateful path you already have:
  return calculateChiSquared_Weighted(const_cast<std::vector<ktlMolecule>&>(mol),
                                      kmin, kmax,
                                      const_cast<std::vector<std::vector<double>>&>(mixtureVals));
}

double experimentalData::scoreWeightedSandboxSingleReplace(const std::vector<ktlMolecule>& current,
                                                           int replaceIndex,
                                                           const ktlMolecule& molNew,
                                                           double &kmin, double &kmax,
                                                           const std::vector<std::vector<double>>& mixtureVals)
{
  std::vector<ktlMolecule> molLike = current;   // shallow copy of handles; deep-enough for your API
  if (replaceIndex >= 0 && replaceIndex < (int)molLike.size()) {
    molLike[replaceIndex] = molNew;
  }
  return scoreWeightedSandbox(molLike, kmin, kmax, mixtureVals);
}
