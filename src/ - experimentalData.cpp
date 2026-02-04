#include "experimentalData.h"


experimentalData::experimentalData(const char* scatterFile){
  noDistBins = 10000;
  /************************************************
    read in the scattering
   ***********************************************/
  std::ifstream scatterdat;
  scatterdat.open(scatterFile);
  std::string sctline;
  double kval; double Ival;double Eval;

  if(scatterdat.is_open()){
     while(!scatterdat.eof()){
      std::getline(scatterdat,sctline);
      std::stringstream ss(sctline);
      ss>>kval;
      ss.ignore();
      ss>>Ival;
      ss.ignore();
      ss>>Eval;

      // std::cout << "\n kval: " << kval << ", Ival: " << Ival  << "\n";
      std::vector<double> ln;
      ln.push_back(kval);
      ln.push_back(Ival);
      ln.push_back(Eval);
      scatVec.push_back(ln);
      exprQset.push_back(kval);
      exprIset.push_back(Ival);
      exprEset.push_back(Eval);
     }
  }else{
    std::cout<<"no scattering file found\n";
  }
  // set the minimum and maximum possibe k values, given the data
  absKmin = scatVec[0][0];
  absKmax = scatVec[scatVec.size()-1][0];

    form_factors['B']= {6.07295, 6.13716, 6.20118, 6.26095, 6.31677, 6.3603, 6.39493, 6.42365, 6.42927, 6.40175, 6.34744, 6.26177, 6.14727, 6.01218, 5.85868, 5.69907, 5.55587, 5.42385, 5.28944, 5.17487, 5.08316};
    form_factors['A']= {8.46551, 8.46612, 8.46554, 8.46411, 8.45628, 8.45991, 8.48668, 8.54026, 8.59042, 8.6273, 8.63936, 8.63948, 8.63224, 8.63738, 8.69052, 8.79047, 8.95796, 9.19116, 9.4747, 9.77982, 10.09104};
    form_factors['R']= {15.97302, 15.95439, 15.93992, 15.93458, 15.96004, 16.00658, 16.08708, 16.21123, 16.35174, 16.50562, 16.67811, 16.87833, 17.10721, 17.37167, 17.66991, 18.01768, 18.41591, 18.85884, 19.3417, 19.84922, 20.36093};
    form_factors['N']= {13.99426, 13.91744, 13.84798, 13.78002, 13.72236, 13.67659, 13.65847, 13.70364, 13.80423, 13.96873, 14.21251, 14.53057, 14.90458, 15.31608, 15.74217, 16.16634, 16.57062, 16.97316, 17.3803, 17.78838, 18.18951};
    form_factors['D']= {18.85878, 18.7376, 18.62254, 18.50494, 18.38204, 18.22496, 18.04046, 17.83065, 17.60334, 17.37469, 17.1935, 17.07013, 17.01717, 17.03528, 17.12922, 17.30479, 17.54592, 17.84864, 18.18575, 18.54153, 18.90354};
    form_factors['C']= {10.88972, 10.86846, 10.85705, 10.84811, 10.85066, 10.88551, 10.95716, 11.05368, 11.15169, 11.24243, 11.34001, 11.4258, 11.50652, 11.58231, 11.62641, 11.63897, 11.62479, 11.60013, 11.57779, 11.56811, 11.56486}; 
    form_factors['Q']= {12.54351, 12.502, 12.46029, 12.41088, 12.36713, 12.33481, 12.33067, 12.39231, 12.52887, 12.74887, 13.05905, 13.44563, 13.89615, 14.39604, 14.92818, 15.48414, 16.04007, 16.59412, 17.13793, 17.68049, 18.21785};
    form_factors['E']= {21.38957, 21.22922, 21.05806, 20.8519, 20.59365, 20.27548, 19.90859, 19.52861, 19.16213, 18.8474, 18.62646, 18.5265, 18.56544, 18.7526, 19.08349, 19.54812, 20.14113, 20.82622, 21.56534, 22.33433, 23.11584};
    form_factors['G']= {11.23402, 11.21512, 11.19615, 11.1726, 11.13785, 11.08978, 11.02476, 10.9486, 10.8257, 10.66288, 10.48178, 10.30495, 10.15137, 10.03115, 9.95346, 9.91543, 9.93208, 10.00603, 10.12878, 10.28136, 10.45941};
    form_factors['H']= {15.26566, 15.14504, 15.02802, 14.8994, 14.75838, 14.60824, 14.45493, 14.31679, 14.18449, 14.06403, 13.9742, 13.92349, 13.92036, 13.96552, 14.06403, 14.2039, 14.38556, 14.60109, 14.83937, 15.09678, 15.35291};
    form_factors['I']= {0.06703, 0.04749, 0.03185, 0.01767, 0.01464, 0.00784, 0.00968, 0.02409, 0.0546, 0.08487, 0.12603, 0.16871, 0.19876, 0.21761, 0.23043, 0.24545, 0.2674, 0.28691, 0.30337, 0.31387, 0.34518};
    form_factors['L']= {0.17687, 0.13367, 0.09214, 0.05337, 0.02589, 0.00433, 0.0, 0.0, 0.00351, 0.00229, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00338, 0.0119, 0.01411, 0.02238, 0.04727, 0.0893};
    form_factors['K']= {10.6941, 10.62362, 10.55538, 10.48627, 10.41843, 10.33769, 10.24058, 10.16454, 10.10982, 10.08123, 10.09476, 10.13953, 10.20381, 10.28853, 10.39257, 10.50933, 10.62698, 10.74104, 10.84485, 10.9335, 11.0132};
    form_factors['M']= {8.14354, 8.15392, 8.17153, 8.19493, 8.23621, 8.28633, 8.34606, 8.42172, 8.49896, 8.56944, 8.64216, 8.69808, 8.73262, 8.73494, 8.69461, 8.61863, 8.50121, 8.35641, 8.19502, 8.02893, 7.85267};
    form_factors['F']= {0.01564, 0.01366, 0.01471, 0.00907, 0.0137, 0.02166, 0.05748, 0.13515, 0.2191, 0.31089, 0.40335, 0.49465, 0.57905, 0.65076, 0.71762, 0.78588, 0.8564, 0.94182, 1.05262, 1.18754, 1.3387};
    form_factors['P']= {4.43189, 4.29078, 4.15107, 4.00949, 3.87962, 3.75139, 3.63015, 3.5358, 3.46906, 3.44222, 3.47424, 3.57316, 3.73358, 3.94682, 4.19912, 4.48073, 4.78045, 5.10558, 5.45824, 5.82542, 6.19645};
    form_factors['S']= {9.90277, 9.8204, 9.7353, 9.64759, 9.54248, 9.41374, 9.29078, 9.20646, 9.17263, 9.20734, 9.32111, 9.51194, 9.76299, 10.05826, 10.38308, 10.72841, 11.07408, 11.40136, 11.70828, 11.99646, 12.27696};
    form_factors['T']={6.9251, 6.88782, 6.84852, 6.79931, 6.75147, 6.70858, 6.67011, 6.64255, 6.62106, 6.61326, 6.63036, 6.67231, 6.7329, 6.8272, 6.96443, 7.13604, 7.33751, 7.56716, 7.81268, 8.07158, 8.34234};
    form_factors['W']= {5.61551, 5.71504, 5.82552, 5.94816, 6.10656, 6.30016, 6.52076, 6.77241, 7.02658, 7.28262, 7.54053, 7.78763, 8.01592, 8.22504, 8.4026, 8.5429, 8.65547, 8.74182, 8.81065, 8.87705, 8.94489};
    form_factors['Y']= {6.3531, 6.36018, 6.38203, 6.41262, 6.4588, 6.52538, 6.62405, 6.76144, 6.91401, 7.05843, 7.19746, 7.32661, 7.43155, 7.50282, 7.54997, 7.58132, 7.59015, 7.59654, 7.59816, 7.60097, 7.60695};
    form_factors['V'] ={0.26858, 0.19696, 0.1278, 0.06453, 0.02555, 0.00256, 0.00026, 0.00495, 0.00745, 0.01011, 0.01635, 0.03125, 0.05043, 0.06951, 0.09509, 0.11501, 0.11912, 0.11049, 0.10911, 0.10318, 0.10142};


  
   
    
    // set up "idealised q values for Josh's coefficients (to change maybe)"

    for (int i = 0; i <= 20; i++) {
        q.push_back(0.01 * i);
    }
    
}


/***********************************************

For interpolating the fitting law to experimental data

 **********************************************/

void experimentalData::calculate_spline_coefficients(const std::vector<double>& x, const std::vector<double>& y,
                                 std::vector<double>& A, std::vector<double>& B,
                                 std::vector<double>& C, std::vector<double>& D) {
    size_t n = x.size();
    std::vector<double> h(n - 1), alpha(n - 1), l(n), mu(n), z(n);
    
    // Step 1: Compute h_i = x_{i+1} - x_i
    for (size_t i = 0; i < n - 1; i++) {
        h[i] = x[i + 1] - x[i];
    }
    
    // Step 2: Compute alpha_i
    for (size_t i = 1; i < n - 1; i++) {
        alpha[i] = (3.0/h[i])*(y[i+1] - y[i]) - (3.0/h[i-1])*(y[i] - y[i-1]);
    }
    
    // Step 3: Solve tridiagonal system
    l[0] = 1.0;
    mu[0] = 0.0;
    z[0] = 0.0;
    
    for (size_t i = 1; i < n - 1; i++) {
        l[i] = 2.0 * (x[i+1] - x[i-1]) - h[i-1]*mu[i-1];
        mu[i] = h[i]/l[i];
        z[i] = (alpha[i] - h[i-1]*z[i-1])/l[i];
    }
    
    l[n-1] = 1.0;
    z[n-1] = 0.0;
    
    // Step 4: Compute coefficients
    std::vector<double> c(n);
    c[n-1] = 0.0;
    
    for (size_t j = n - 2; j < n; j--) {  // Note: using size_t, so we break when j underflows
        c[j] = z[j] - mu[j]*c[j+1];
    }
    
    // Step 5: Compute spline coefficients
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
    // Find the appropriate interval
    size_t i = 0;
    while (i < x.size() - 2 && x[i + 1] <= x_val) i++;
    
    // Compute relative x
    double dx = x_val - x[i];
    
    // Evaluate cubic polynomial
    return ((A[i]*dx + B[i])*dx + C[i])*dx + D[i];
}





bool experimentalData::binDataCheck(double &dMax,double &qmin,double &qmax){
  double dr = dMax/double(noDistBins);
  double dq = (qmax-qmin)/double(noDistBins);

  bool goodsplit=true;
  int k=0;
  for(int j=1;j<=noDistBins;j++){
    double q = qmin + (j-0.5)*dq;
    double qminBin = qmin + (j-1)*dq;
    double qmaxBin = qmin + j*dq;
    std::vector<double> intensities;
    // if we have selected a higher q than the lowest experimental data, seacrch for the minimum point
    while(scatVec[k][0]<qminBin){
      k++;
    }
    // std::cout << "\n j: " << j << ", qminBin: " << qminBin << ", qmaxBin: " << qmaxBin <<", scatVec[k] first: " << scatVec[k].first << ", scatVec[k] second: " << scatVec[k].second << "\n";

    while(scatVec[k][0]>=qminBin &&scatVec[k][0]<=qmaxBin){

          intensities.push_back(scatVec[k][1]);
          k++;
    }
    if(intensities.size()<2){
      goodsplit=false;
    }

    // std::cout << "\n j: " << j << ", intensity size: " << intensities.size() << "\n";

  }
  return  goodsplit;
}

void experimentalData::subsetScatteringData(std::vector<double>& A,
                      std::vector<double>& B,
                      std::vector<double>& C,
                      double kmin, double kmax,
                      std::vector<double>& A_selected,
                      std::vector<double>& B_selected,
                      std::vector<double>& C_selected)
{
    A_selected.clear();
    B_selected.clear();
    C_selected.clear();

    for (size_t i = 0; i < A.size(); ++i) {
        if (A[i] >= kmin && A[i] <= kmax) {
            A_selected.push_back(A[i]);
            B_selected.push_back(B[i]);
            C_selected.push_back(C[i]);
        }
    }
}


/****************************

Calculate the range of q values using Shannon bins*2 nad calculate the experimental intensities for this

 **************************/

int experimentalData::setPhases(double &dMax,double &qmin,double &qmax){
  int noDistBinsTemp = int(2.1*std::ceil((qmax-qmin)*dMax/3.14159265359));
  // std::cout<<"\n no dist bins tmp " << noDistBinsTemp << "\n";
  // std::cout<<"\n no dist bins glb " << noDistBins << "\n";

  if(noDistBins!=noDistBinsTemp){
    noDistBins=noDistBinsTemp;
    //std::cout<<"\n no dist bins A " << noDistBins << "\n";

    //check if we can bin
    bool binsOk=false;
    // check we can bin the data with this number of bins
    while(binsOk==false && noDistBins>2){
      binsOk = binDataCheck(dMax,qmin,qmax);
      if(binsOk==false){
	noDistBins=noDistBins-1;
	// std::cout<<"\n no dist bins check: " << noDistBins << "\n";

      }
    }
    //std::cout<<"\n no dist bins " << dMax <<" " << qmin <<" " << qmax <<" "<<noDistBins<<"\n";
    double dr = dMax/double(noDistBins);
    double dq = (qmax-qmin)/double(noDistBins);
    experimentalIntensity.clear();distBins.clear();
    int k=0;
    qvals.clear();
    for(int i=1;i<=noDistBins;i++){
      std::vector<double> rset;
      double r = (i-0.5)*dr;
      std::pair<double,double> rDistSets;
      rDistSets.first=(i-1)*dr;
      rDistSets.second=i*dr;
      distBins.push_back(rDistSets);
      double q;
      for(int j=1;j<=noDistBins;j++){
	q =qmin+(j-0.5)*dq;
	if(i==1){
	  qvals.push_back(q);
	  double qminBin = qmin+(j-1)*dq;
	  double qmaxBin = qmin+j*dq;
	  std::vector<double> intensities;
	  std::vector<double> errors; 
	  // for the logged standard deviation.
	  double intensitySum = 0.0;
	  int noIntensities=0;
	  std::vector<double> loggedIntensities;
	  // if we have selected a higher q than the lowest experimental data, seacrch for the minimum point
	  while(scatVec[k][0]<qminBin){
	    k++;
	  }
	  while(scatVec[k][0]>=qminBin &&scatVec[k][0]<=qmaxBin){
	    intensities.push_back(scatVec[k][1]);
	    if(scatVec[k][1]>0.0){
	      double loggedInten = std::log(scatVec[k][1]);
	      intensitySum = intensitySum + loggedInten;
	      noIntensities++;
	      loggedIntensities.push_back(loggedInten);
	    }
	    k++;
	  }
	  //get the median
	  std::sort(intensities.begin(),intensities.end());
	  int size = intensities.size();
	  int medPt =int(std::round(float(size)/2.0));
	  experimentalIntensity.push_back(intensities[medPt-1]);
	  // calulate the mean
	  double mean = intensitySum/double(noIntensities);
	}
      }
    }
  }
  // now take the median of the scattering data
  // std::cout<<"\n fin setPhases * \n";

  //finally  chop the experimental data to this range
  subsetScatteringData(exprQset,exprIset,exprEset,qmin,qmax,exprQSubset,exprISubset,exprESubset);
  return noDistBins;
}

/**************************************************************

Calculate the predcited normal directions

 ************************************************************/


std::vector<point> experimentalData::calculate_geometric_normals(std::vector<point>& ca_coords) {
    
    std::vector<point> ca_vectors;
    std::vector<point> normals;
    
    // Calculate vectors between consecutive CA atoms
    for (size_t i = 1; i < ca_coords.size(); i++) {
        ca_vectors.push_back(ca_coords[i] - ca_coords[i-1]);
    }
    
    // calculate the geometric *outward* facing normals
    for (size_t i = 1; i < ca_vectors.size(); i++) {
        
        point norm = ca_vectors[i] - ca_vectors[i-1];
	
        norm = norm * (-1.0) / norm.length();   // outward
        normals.push_back(norm);
    }
    
    point first_norm = ca_vectors[0] * (-1.0) / ca_vectors[0].length();
    point final_norm = ca_vectors.back() / ca_vectors.back().length();
    
    
    normals.insert( normals.begin(), first_norm );
    normals.push_back( final_norm );

    
    
    return normals;
}


/****************************************************************

    Place the sidechains 

*****************************************************************/


std::vector<point> experimentalData::place_side_chains( std::vector<point>& ca_coords,
                                      std::vector<point>& geometric_vectors,
                                      std::vector<char>& residue_names) {
    
    std::map<char, double> residue_distances = {
        {'R', 4.2662}, {'N', 2.5349}, {'D', 2.5558}, {'C', 2.3839},
        {'Q', 3.1861}, {'E', 3.2541}, {'H', 3.1861}, {'I', 2.3115},
        {'L', 2.6183}, {'K', 3.6349}, {'M', 3.1912}, {'F', 3.4033},
        {'P', 1.8773}, {'U', 1.5419}, {'S', 1.9661}, {'T', 1.9533},
        {'W', 3.8916}, {'Y', 3.8807}, {'V', 1.9555}, {'G', 0.0}, {'A', 0.0}
    };
    
    std::vector<point> side_chain_positions;
    
    for (size_t i = 0; i < ca_coords.size(); i++) {
        
        double distance = residue_distances[residue_names[i]];
        side_chain_positions.push_back(ca_coords[i] + geometric_vectors[i] * distance);
        
    }
    
    return side_chain_positions;
}

std::vector<std::vector<double> > experimentalData::calculate_distances( std::vector<point>& coordinates, int& molIndex) {
    
    size_t n = coordinates.size();
    std::vector< std::vector<double> > distances(n, std::vector<double>(n));
    double min_dist = std::numeric_limits<double>::max();
    
    for (size_t i = 0; i < n; i++){
        for (size_t j = i; j < n; j++){
            
            if (i==j) { distances[i][j] = 0.0;
                
            } else {
                double dist = coordinates[i].eDist(coordinates[j]);
		if(dist>maxDist[molIndex]){
		  //std::cout<<"new biggest dist "<<dist<<"\n";
		  //coordinates[i].printPoint();
		  //coordinates[j].printPoint();
		  maxDist[molIndex] = dist;
		}
                distances[i][j] = dist;
                distances[j][i] = dist;
            }
        }
    }
    
    return distances;
}

 
/*********************************************************************


calculate all the scattering centres

 ********************************************************************/

ScatteringCenters experimentalData::process_structure(std::vector<point>& ca_coords,
						      std::vector<char>& residue_names,int& molIndex){
    
    std::vector<point> geometric_normals = calculate_geometric_normals(ca_coords);
    std::vector<point> side_chain_positions = place_side_chains(ca_coords, geometric_normals, residue_names);
        
    ScatteringCenters centers;
   
  
    // add backbone + sidechains centres for non-GLY/ALA residues
    for (size_t i=0; i < residue_names.size(); i++) {
        
        if (residue_names[i] != 'A' && residue_names[i] != 'G') {
            
            // backbone
            centers.coordinates.push_back(ca_coords[i]);
            centers.types.push_back('B');

	    
            // side chain
            centers.coordinates.push_back(side_chain_positions[i]);
	   
            centers.types.push_back(residue_names[i]);
        }
    }
    
    // add GLY + ALA CA postions
    for (size_t i=0; i < residue_names.size(); i++) {
        
        if (residue_names[i] == 'A' || residue_names[i] == 'G') {
            
            // backbone
            centers.coordinates.push_back(ca_coords[i]);
            centers.types.push_back(residue_names[i]);
        }
    }
    
    centers.distances = calculate_distances(centers.coordinates,molIndex);
    
    return centers;
}

/****************************************************

Routine to calculatge predicted I

 ***************************************************/
 

std::vector<double> experimentalData::calculate_saxs_implicit( ScatteringCenters& centers ) {
    
    std::vector<double> I(q.size(), 0.0);
    
    std::vector< std::vector<double> > ff_matrix;
    for (char type : centers.types) {
        ff_matrix.push_back( form_factors.at(type) );
    }
    // calculate scattering at each q
    for (int q_idx = 0; q_idx < q.size(); q_idx++) {
        
        double qv =q[q_idx];
	
        for (size_t i = 0; i < centers.coordinates.size(); i++) {
            for (size_t j = 0; j < centers.coordinates.size(); j++) {
                double r_ij = centers.distances[i][j];
                double ff_product = ff_matrix[i][q_idx] * ff_matrix[j][q_idx];
                
                double sinc;
                if (qv<0.0000000001) {
                    sinc = 1.0;
                } else {
                    sinc = (r_ij == 0) ? 1.0 : std::sin(qv * r_ij) / (qv * r_ij);
                }
                
                I[q_idx] += ff_product * sinc;
            }
        }
    }
    
    return I;
}


/********************************************************************

flatten_coords put a molecule into one vector of coordinates for the scattering calc

 *******************************************************************/

std::vector<point> experimentalData::flatten_coords( std::vector<std::vector<point> >& coords) {
   std::vector<point> flat;
   for (const auto& section : coords) {
       flat.insert(flat.end(), section.begin(), section.end());
   }
   return flat;
}


std::vector<char> convertStringsToChars(const std::vector<std::string>& strings) {
    return std::accumulate(strings.begin(), strings.end(), std::vector<char>(),
        [](std::vector<char>& result, const std::string& str) {
            result.insert(result.end(), str.begin(), str.end());
            return result;
        });
}

std::vector<char> experimentalData::flatten_residueNames( std::vector<std::vector<std::string> >& aminoList) {
   std::vector<char> flat;
   for (const auto& section : aminoList) {
     std::vector<char> charSec =convertStringsToChars(section);
     flat.insert(flat.end(), charSec.begin(),charSec.end());
   }
   return flat;
}

/*********************************************

For interpolation to experimental data

 *********************************************/

std::vector<double> experimentalData::calculate_intensity_at_experimental_q(std::vector<double>& I_mod) {
    
    size_t N = q.size() - 1;
    std::vector<double> A(N), B(N), C(N), D(N);
    
    calculate_spline_coefficients(q, I_mod, A, B, C, D);
    
    std::vector<double> I_interp;
    for(int i=0;i<qvals.size();i++){   
      I_interp.push_back(evaluate_spline(qvals[i], q, A, B, C, D));  
    }
    return I_interp;
}


 double experimentalData::calculateChiSquared(std::vector<ktlMolecule> &mol,double &qmin,double &qmax,std::vector<std::vector<double> > &mixtureVals){
    kMin = qmin;
   kMax = qmax;
   for(int i=0;i<mol.size();i++){
     maxDist.push_back(0.0);
   }
   // loop over all molecules considered and calculate their side chains and distances
   for(int i=0;i<mol.size();i++){
     //get calpha's
     std::vector<std::vector<point> > coords=mol[i].getCoordinates();
     std::vector<point> chainCoordinates = flatten_coords(coords);
     // get side chains
     std::vector<std::vector<std::string> > aminoList = mol[i].getAminoList();
     std::vector<char> chainResidues = flatten_residueNames(aminoList);
     ScatteringCenters calphsAndSideChains  = process_structure(chainCoordinates,chainResidues,i);
     scs.push_back(calphsAndSideChains);
   }

   // use the maximum distance size to set the number of q values we fit to:
   double dMax = *std::max_element(begin(maxDist), end(maxDist));
   int noBins = setPhases(dMax,qmin,qmax);
   
   // now loop over all molecules and calculate their scattering

   Ivec.clear();
   for(int i=0;i<mol.size();i++){
     Ivec.push_back(calculate_saxs_implicit(scs[i]));
   }

   IvecOnData.clear();
   /*for(int i=0;i<mol.size();i++){
     IvecOnData.push_back(calculate_intensity_at_experimental_q(Ivec[i]));
     }*/


   for(int i=0;i<mol.size();i++){
     IvecOnData.push_back( calculate_intensity_at_experimental_q(q,Ivec[i],exprQSubset));
   }
   

   // interpolate the experimental data at the shannon

   exprIInterp.clear();
   
   exprIInterp =calculate_intensity_at_experimental_q(qvals,experimentalIntensity,exprQSubset);
   
   double pred=10000.0;

   for(int i =0;i<mixtureVals.size();i++){
     std::vector<double> Icomb(IvecOnData[0].size(),0.0);
     for(int j =0;j<mixtureVals[i].size();j++){
       for(int k =0;k<IvecOnData[j].size();k++){
	 Icomb[k] = Icomb[k] + IvecOnData[j][k]*mixtureVals[i][j];
       }
     }
     // calculate the "chi squared fit" first calculae all the distances and work out the scale factor
     double logDifMean=0;
     int noMean=0;
     std::vector<double> logdifs;
     
     for(int l=0;l<Icomb.size();l++){
    //std::cout<<std::log(scatVals[i])+(std::log(experimentalIntensity[0])-std::log(scatVals[0]))<<" "<<std::log(experimentalIntensity[i])<<"\n";
       double logScatDif= std::log(Icomb[l]) - std::log(exprIInterp[l]);
       logdifs.push_back(logScatDif);
       if(exprQSubset[l]<kMin+0.01){
	 logDifMean =  logDifMean + logScatDif;
	 noMean++;
       }
     }
     logDifMean = logDifMean/double(noMean);

     // finally calculate the "chi squared value "
     double predTemp =0.0;
     for(int l=0;l<Icomb.size();l++){
       double scatDif = logdifs[l] - logDifMean;
       predTemp = predTemp + scatDif*scatDif;
     }
     if(predTemp<pred){
       pred = predTemp;
     }
   }
   //std::cout<<"in ed norm "<<pred/(IvecOnData[0].size()-1)<<"\n";
   return pred/(exprQSubset.size()-1);
 }





double experimentalData::calculateChiSquaredUpdate(ktlMolecule& molNew,int& k,double &qmin,double &qmax,std::vector<std::vector<double> > &mixtureVals){
   kMin = qmin;
   kMax = qmax;
   // loop over all molecules considered and calculate their side chains and distances
   maxDist[k]=0.0;
   std::vector<std::vector<point> > coords=molNew.getCoordinates();
   std::vector<point> chainCoordinates = flatten_coords(coords);
     // get side chains
   std::vector<std::vector<std::string> > aminoList = molNew.getAminoList();
   std::vector<char> chainResidues = flatten_residueNames(aminoList);
   ScatteringCenters calphsAndSideChains  = process_structure(chainCoordinates,chainResidues,k);
   scs[k] = calphsAndSideChains;
   

   // use the maximum distance size to set the number of q values we fit to:
   double dMax = *std::max_element(begin(maxDist), end(maxDist));
   int noBins = setPhases(dMax,qmin,qmax);
   // now loop over all molecules and calculate their scattering

   Ivec[k]=calculate_saxs_implicit(scs[k]);
   IvecOnData[k] = calculate_intensity_at_experimental_q(q,Ivec[k],exprQSubset);

   exprIInterp.clear();
   exprIInterp = calculate_intensity_at_experimental_q(qvals,experimentalIntensity,exprQSubset);
   double pred=10000.0;

   for(int i =0;i<mixtureVals.size();i++){
     std::vector<double> Icomb(IvecOnData[0].size(),0.0);
     for(int j =0;j<mixtureVals[i].size();j++){
       for(int k =0;k<IvecOnData[j].size();k++){
	 Icomb[k] =  Icomb[k] +  IvecOnData[j][k]*mixtureVals[i][j];
       }
     }
     // calculate the "chi squared fit" first calculae all the distances and work out the scale factor
     double logDifMean=0;
     int noMean =0;
     std::vector<double> logdifs;
     for(int l=0;l<Icomb.size();l++){
    //std::cout<<std::log(scatVals[i])+(std::log(experimentalIntensity[0])-std::log(scatVals[0]))<<" "<<std::log(experimentalIntensity[i])<<"\n";
       double logScatDif= std::log(Icomb[l]) - std::log(exprIInterp[l]);
       logdifs.push_back(logScatDif);
       if(exprQSubset[l]<kMin+0.01){
	 logDifMean =  logDifMean + logScatDif;
	 noMean++;
       }
     }
     logDifMean = logDifMean/double(noMean);
     // finally calculate the "chi squared value "
     double predTemp =0.0;
     for(int l=0;l<Icomb.size();l++){
       double scatDif = logdifs[l] - logDifMean;
       predTemp = predTemp + scatDif*scatDif;
     }
     if(predTemp<pred){
       pred = predTemp;
     }
   }
   return pred/(exprQSubset.size()-1);
 }



void experimentalData::writeScatteringToFile(std::vector<std::vector<double> > &mixtureVals,const char* filename){
  int best =0;
  double logDifMeanFinal = 10000.0;
  double pred =10000.0;
  std::vector<double> IcombBest;
  exprIInterp = calculate_intensity_at_experimental_q(qvals,experimentalIntensity,exprQSubset);
  for(int i =0;i<mixtureVals.size();i++){
    std::vector<double> Icomb(IvecOnData[0].size(),0.0);
    for(int j =0;j<mixtureVals[i].size();j++){
      for(int k =0;k<IvecOnData[j].size();k++){
	Icomb[k] =  Icomb[k] + IvecOnData[j][k]*mixtureVals[i][j];
      }
    }
     // calculate the "chi squared fit" first calculae all the distances and work out the scale factor
    double logDifMean=0;
    int noMean =0;
    std::vector<double> logdifs;
    for(int l=0;l<Icomb.size();l++){
      //std::cout<<std::log(scatVals[i])+(std::log(experimentalIntensity[0])-std::log(scatVals[0]))<<" "<<std::log(experimentalIntensity[i])<<"\n";
      double logScatDif= std::log(Icomb[l]) - std::log(exprIInterp[l]);
      logdifs.push_back(logScatDif);
      if(exprQSubset[l]<kMin+0.01){
	logDifMean =  logDifMean + logScatDif;
	noMean++;
       }
    }
    logDifMean = logDifMean/double(noMean);
    // finally calculate the "chi squared value "
    double predTemp =10000.0;
    for(int l=0;l<Icomb.size();l++){
      double scatDif = logdifs[l] - logDifMean;
      predTemp = predTemp + scatDif*scatDif;
    }
    if(predTemp<pred || i==0){
	pred = predTemp;
	best = i;
	logDifMeanFinal = logDifMean;
	IcombBest= Icomb;
    }
  }
  
  // write to file
  std::ofstream myfile;
  myfile.open(filename);
  for(int i=0;i<IcombBest.size();i++){
    myfile<<exprQSubset[i]<<" "<<IcombBest[i]<<" "<<std::log(IcombBest[i])- logDifMeanFinal<<" "<<std::log(exprISubset[i])<<"\n";
  }
  // finally add the percentage combination values
  for(int i=0;i<mixtureVals[best].size();i++){
    if(i==(mixtureVals[best].size()-1)){
      myfile<<mixtureVals[best][i]<<"\n";
    }else{
      myfile<<mixtureVals[best][i]<<" ";
    }
  }
  myfile.close();
 }



std::vector<double> experimentalData::calculate_intensity_at_experimental_q(std::vector<double>& q_mod,std::vector<double>& I_mod,std::vector<double>& expQRange) {
    
    size_t N = q_mod.size() - 1;
    std::vector<double> A(N), B(N), C(N), D(N);
    calculate_spline_coefficients(q_mod, I_mod, A, B, C, D);
    
    std::vector<double> I_interp;
    for(int i =0; i<expQRange.size();i++) {
      I_interp.push_back(evaluate_spline(expQRange[i], q_mod, A, B, C, D));   
    }
    return I_interp;
}


 double experimentalData::calculateChiSquared_Weighted(std::vector<ktlMolecule> &mol,double &qmin,double &qmax,std::vector<std::vector<double> > &mixtureVals){
   kMin = qmin;
   kMax = qmax;
   for(int i=0;i<mol.size();i++){
     maxDist.push_back(0.0);
   }
   // loop over all molecules considered and calculate their side chains and distances
   for(int i=0;i<mol.size();i++){
     //get calpha's
     std::vector<std::vector<point> > coords=mol[i].getCoordinates();
     std::vector<point> chainCoordinates = flatten_coords(coords);
     // get side chains
     std::vector<std::vector<std::string> > aminoList = mol[i].getAminoList();
     std::vector<char> chainResidues = flatten_residueNames(aminoList);
     ScatteringCenters calphsAndSideChains  = process_structure(chainCoordinates,chainResidues,i);
     scs.push_back(calphsAndSideChains);
   }

   // use the maximum distance size to set the number of q values we fit to:
   double dMax = *std::max_element(begin(maxDist), end(maxDist));
 
   int noBins = setPhases(dMax,qmin,qmax);

   // now loop over all molecules and calculate their scattering

   Ivec.clear();
   for(int i=0;i<mol.size();i++){
     Ivec.push_back(calculate_saxs_implicit(scs[i]));
   }

   IvecOnData.clear();
   
   for(int i=0;i<mol.size();i++){
     IvecOnData.push_back(calculate_intensity_at_experimental_q(q,Ivec[i],exprQSubset));
   }

   
   double pred=10000.0;

   for(int i =0;i<mixtureVals.size();i++){
     std::vector<double> Icomb(IvecOnData[0].size(),0.0);
     for(int j =0;j<mixtureVals[i].size();j++){
       for(int k =0;k<IvecOnData[j].size();k++){
	 Icomb[k] =  Icomb[k] + IvecOnData[j][k]*mixtureVals[i][j];
       }
     }
     // calculate the "chi squared fit" first calculae all the distances and work out the scale factor
     double logDifMean=0;
     int noMean=0;
     std::vector<double> logdifs;
     
     for(int l=0;l<Icomb.size();l++){
       // double logScatDif= std::log(Icomb[l]) - std::log(experimentalIntensity[l]);
       double logScatDif= std::log(Icomb[l]) - std::log(exprISubset[l]);
       logdifs.push_back(logScatDif);
       if(exprQSubset[l]<kMin+0.1){
	 logDifMean =  logDifMean + logScatDif;
	 noMean++;
       }
     }
     std::cout<<logDifMean<<"\n";
     logDifMean = logDifMean/double(noMean);

   
     for(int l=0;l<Icomb.size();l++){
       double logScatScaled = std::log(Icomb[l])-logDifMean;
       double scatScaled = std::exp(logScatScaled);
       Icomb[l] = scatScaled;
     }
     // now (spline interpolate)
     //std::vector<double> Imodel = calculate_intensity_at_experimental_q(qvals,Icomb,exprQSubset);
     std::vector<double> Imodel = Icomb;
     double predTemp =0.0;
    
     for(int l =0;l<exprQSubset.size();l++){
       double scatInterp = Imodel[l];
       double dif = scatInterp - exprISubset[l];
       predTemp = predTemp + dif*dif/(exprESubset[l]*exprESubset[l]);
       //predTemp = predTemp + dif*dif;
     }
     predTemp = std::abs(predTemp/(exprQSubset.size()-1) -1);
     std::cout<<predTemp<<"\n";
     if(predTemp<pred){
      pred = predTemp;
     }
     //std::cout<<"in ed norm (0)"<<pred<<" "<<predTemp<<"\n";
   }
   return pred;
 }

double experimentalData::calculateChiSquaredUpdate_Weighted(ktlMolecule& molNew,int& k,double &qmin,double &qmax,std::vector<std::vector<double> > &mixtureVals){
   // loop over all molecules considered and calculate their side chains and distances
   kMin = qmin;
   kMax = qmax;
   maxDist[k]=0.0;
   std::vector<std::vector<point> > coords=molNew.getCoordinates();
   std::vector<point> chainCoordinates = flatten_coords(coords);
     // get side chains
   std::vector<std::vector<std::string> > aminoList = molNew.getAminoList();
   std::vector<char> chainResidues = flatten_residueNames(aminoList);
   ScatteringCenters calphsAndSideChains  = process_structure(chainCoordinates,chainResidues,k);
   scs[k] = calphsAndSideChains;

   // use the maximum distance size to set the number of q values we fit to:
   double dMax = *std::max_element(begin(maxDist), end(maxDist));
   int noBins = setPhases(dMax,qmin,qmax);
   // now loop over all molecules and calculate their scattering

   Ivec[k]=calculate_saxs_implicit(scs[k]);
   IvecOnData[k] = calculate_intensity_at_experimental_q(q,Ivec[k],exprQSubset);
   
   double pred=10000.0;

   for(int i =0;i<mixtureVals.size();i++){
     std::vector<double> Icomb(IvecOnData[0].size(),0.0);
     for(int j =0;j<mixtureVals[i].size();j++){
       for(int k =0;k<IvecOnData[j].size();k++){
	 Icomb[k] =  Icomb[k] + IvecOnData[j][k]*mixtureVals[i][j];
       }
     }
     // calculate the "chi squared fit" first calculae all the distances and work out the scale factor
     double logDifMean=0;
     int noMean=0;
     std::vector<double> logdifs;
     
     for(int l=0;l<Icomb.size();l++){
       // double logScatDif= std::log(Icomb[l]) - std::log(experimentalIntensity[l]);
       double logScatDif= std::log(Icomb[l]) - std::log(exprISubset[l]);
       logdifs.push_back(logScatDif);
       if(exprQSubset[l]<kMin+0.1){
	 logDifMean =  logDifMean + logScatDif;
	 noMean++;
       }
     }
     logDifMean = logDifMean/double(noMean);
   
   
     for(int l=0;l<Icomb.size();l++){
       double logScatScaled = std::log(Icomb[l])-logDifMean;
       double scatScaled = std::exp(logScatScaled);
       Icomb[l] = scatScaled;
     }
     // now (spline interpolate)
     //std::vector<double> Imodel = calculate_intensity_at_experimental_q(qvals,Icomb,exprQSubset);
     std::vector<double> Imodel = Icomb;
     double predTemp =0.0;
    
     for(int l =0;l<exprQSubset.size();l++){
       double scatInterp = Imodel[l];
       double dif = scatInterp - exprISubset[l];
       //std::cout<<l<<" "<<dif<<" "<<exprESubset[l]<<"\n";
       predTemp = predTemp + dif*dif/(exprESubset[l]*exprESubset[l]);
       //predTemp = predTemp + dif*dif;
     }
     predTemp = std::abs(predTemp/(exprQSubset.size()-1) -1);
     if(predTemp<pred){
      pred = predTemp;
     }
     //std::cout<<"in ed update"<<pred<<" "<<predTemp<<"\n";
   }
   return pred;
 }


void experimentalData::writeScatteringToFile_ChiSq(std::vector<std::vector<double> > &mixtureVals,const char* filename){
  int best =0;
  double logDifMeanFinal = 10000.0;
  double pred =10000.0;
  std::vector<double> IcombBest;
  for(int i =0;i<mixtureVals.size();i++){
    std::vector<double> Icomb(IvecOnData[0].size(),0.0);
    for(int j =0;j<mixtureVals[i].size();j++){
      for(int k =0;k<IvecOnData[j].size();k++){
	Icomb[k] =  Icomb[k] + IvecOnData[j][k]*mixtureVals[i][j];
      }
    }
     double logDifMean=0;
     int noMean =0;
     std::vector<double> logdifs;
     for(int l=0;l<Icomb.size();l++){
       // double logScatDif= std::log(Icomb[l]) - std::log(experimentalIntensity[l]);
       double logScatDif= std::log(Icomb[l]) - std::log(exprISubset[l]);
       logdifs.push_back(logScatDif);
       if(exprQSubset[l]<kMin+0.1){
	 logDifMean =  logDifMean + logScatDif;
	 noMean++;
       }
     }
     logDifMean = logDifMean/double(noMean);
     for(int l=0;l<Icomb.size();l++){
       double logScatScaled = std::log(Icomb[l])-logDifMean;
       double scatScaled = std::exp(logScatScaled);
       Icomb[l] = scatScaled;
     }
     // now (spline interpolate)
     std::vector<double> Imodel = Icomb;
    // finally calculate the "chi squared value "
    pred=100000.0;
    double predTemp =0.0;
    for(int l =0;l<exprQSubset.size();l++){
      double scatInterp = Imodel[l];
      double dif = scatInterp - exprISubset[l];
      predTemp = predTemp + dif*dif/(exprESubset[l]*exprESubset[l]);
    }
    predTemp = std::abs(predTemp/(exprQSubset.size()-1) -1);
    if(predTemp<pred){
      pred = predTemp;
      IcombBest= Imodel;
	  best = i;	
    }
  }
  // write to file
  std::ofstream myfile;
  myfile.open(filename);
  for(int l =0;l<exprQSubset.size();l++){
    myfile<<exprQSubset[l]<<" "<<IcombBest[l]<<" "<<std::log(IcombBest[l])<<" "<<std::log(exprISubset[l])<<"\n";
  }
  // finally add the percentage combination values
  for(int i=0;i<mixtureVals[best].size();i++){
    if(i==(mixtureVals[best].size()-1)){
      myfile<<mixtureVals[best][i]<<"\n";
    }else{
      myfile<<mixtureVals[best][i]<<" ";
    }
  }
  myfile.close();
 }
