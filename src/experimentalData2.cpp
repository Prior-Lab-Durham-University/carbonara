#include "experimentalData.h"

experimentalData::experimentalData(const char* scatterFile){
  noDistBins =10000;
  /************************************************
    read in the scattering 
   ***********************************************/
  std::ifstream scatterdat;
  scatterdat.open(scatterFile);
  std::string sctline;
  double kval;double Ival;
  if(scatterdat.is_open()){
     while(!scatterdat.eof()){
      std::getline(scatterdat,sctline);
      std::stringstream ss(sctline);
      ss>>kval;
      ss.ignore();
      ss>>Ival;
      std::pair<double,double> pr;
      pr.first=kval;pr.second=Ival;
      scatVec.push_back(pr);
     }
  }else{
    std::cout<<"no scattering file found\n";
  }
  // set the minimum and maximum possibe k values, given the data
  absKmin = scatVec[0].first;
  absKmax = scatVec[scatVec.size()-1].first;
  /************************
     load in the scattering function paramers 
   ***********************/
  std::ifstream myfile;
  std::string output;
  myfile.open("averagedSolutions/3_4BeadParams.dat");
  double val;
  if (myfile.is_open()) { 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
      ss>>val;
      gaussianParams.push_back(val);
    }
  }
  a1=gaussianParams[0];
  a2=gaussianParams[1];
  a3=gaussianParams[2];
  a4=gaussianParams[3];
  a5=gaussianParams[4];
  b1=gaussianParams[5];
  b2=gaussianParams[6];
  b3=gaussianParams[7];
  b4=gaussianParams[8];
  b5=gaussianParams[9];
  c=gaussianParams[10];
  myfile.close();
  myfile.open("averagedSolutions/Hyd3_4BeadParams.dat");
  if (myfile.is_open()) { 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
      ss>>val;
      gaussianParamsHyd.push_back(val);
    }
  }
  myfile.close();
  a1H=gaussianParamsHyd[0];
  a2H=gaussianParamsHyd[1];
  a3H=gaussianParamsHyd[2];
  a4H=gaussianParamsHyd[3];
  a5H=gaussianParamsHyd[4];
  b1H=gaussianParamsHyd[5];
  b2H=gaussianParamsHyd[6];
  b3H=gaussianParamsHyd[7];
  b4H=gaussianParamsHyd[8];
  b5H=gaussianParamsHyd[9];
  cH=gaussianParamsHyd[10];
}

double experimentalData::gaussianFixed(double &s){
  double out =0.0;
  out = out + a1*std::exp(-b1*s*s);
  out = out + a2*std::exp(-b2*s*s);
  out = out + a3*std::exp(-b3*s*s);
  out = out + a4*std::exp(-b4*s*s); 
  out = out + a5*std::exp(-b5*s*s); 
  out = out +c;
  return out;
}

double experimentalData::gaussianHydFixed(double &s){
  double out =0.0;
  out = out + a1H*std::exp(-b1H*s*s);
  out = out + a2H*std::exp(-b2H*s*s);
  out = out + a3H*std::exp(-b3H*s*s);
  out = out + a4H*std::exp(-b4H*s*s); 
  out = out + a5H*std::exp(-b5H*s*s); 
  out = out +cH;
  return out;
}

double experimentalData::hydrationScatter(double &s){
  double fhydro = 0.702260*std::exp(-23.945605*s*s)+0.763666*std::exp(-74.897918*s*s)+0.248678*std::exp(-6.773289*s*s)+0.261323*std::exp(-233.583447*s*s) + 0.023017*std::exp(-1.337531*s*s)+0.000425;
  double foxy = 2.960427*std::exp(-14.182259*s*s)+2.508818*std::exp(-5.936858*s*s)+0.637853*std::exp(-0.112726*s*s)+0.722838*std::exp(-34.958481*s*s) + 1.142756*std::exp(-0.390240*s*s)+0.027014;
  return 2.0*fhydro + foxy;
}


bool experimentalData::binDataCheck(double &dMax,double &qmin,double &qmax){
  double dr = dMax/double(noDistBins);
  double dq = (qmax-qmin)/double(noDistBins);
  bool goodsplit=true;
  int k=0;
  for(int j=1;j<=noDistBins;j++){
    double q =qmin+(j-0.5)*dq;
    double qminBin = qmin+(j-1)*dq;
    double qmaxBin = qmin+j*dq;
    std::vector<double> intensities;
    // if we have selected a higher q than the lowest experimental data, seacrch for the minimum point
    while(scatVec[k].first<qminBin){
      k++;
    }
    while(scatVec[k].first>qminBin &&scatVec[k].first<=qmaxBin){
          intensities.push_back(scatVec[k].second);
          k++;
    }
    if(intensities.size()<2){
      goodsplit=false;
    }
  }
  return  goodsplit;
}


int experimentalData::setPhases(double &dMax,double &qmin,double &qmax){
  int noDistBinsTemp = int(2.1*std::ceil((qmax-qmin)*dMax/3.14159265359));
  if(noDistBins!=noDistBinsTemp){
    noDistBins=noDistBinsTemp;
    //check if we can bin
    bool binsOk=false;
    // check we can bin the data with this number of bins
    //std::cout<<"no dist bins "<<dMax<<" "<<noDistBins<<"\n";
    while(binsOk==false && noDistBins>2){
      binsOk= binDataCheck(dMax,qmin,qmax);
      if(binsOk==false){
	noDistBins=noDistBins-1;
      }
    }
    //std::cout<<"no dist bins "<<dMax<<" "<<noDistBins<<"\n";
    double dr = dMax/double(noDistBins);
    double dq = (qmax-qmin)/double(noDistBins);
    experimentalIntensity.clear();distBins.clear();scatPhases.clear();
    experimentalSD.clear();
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
	  // for the logged standard deviation. 
	  double intensitySum = 0.0;
	  int noIntensities=0;
	  std::vector<double> loggedIntensities;
	  // if we have selected a higher q than the lowest experimental data, seacrch for the minimum point
	  while(scatVec[k].first<qminBin){
	    k++;
	  }
	  while(scatVec[k].first>qminBin &&scatVec[k].first<=qmaxBin){
	    intensities.push_back(scatVec[k].second);
	    if(scatVec[k].second>0.0){
	      double loggedInten = std::log(scatVec[k].second);
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
	  // and the standard deviation
	  double SD=0.0;
	  for(int k=0;k<loggedIntensities.size();k++){
	    double dif = loggedIntensities[k]-mean;
	    SD=SD+dif*dif;
	  }					    
	  experimentalSD.push_back(SD/(double(loggedIntensities.size()-1)));
	}     
	double krprod = q*r;
	rset.push_back((std::sin(krprod)/krprod));
      }
      scatPhases.push_back(rset);
    }
  }
  // now take the median of the scattering data
  return noDistBins;
}

double experimentalData::fitToScattering(std::vector<double> &molDists,std::vector<double> &solDists,std::vector<double> &solMolDists,int &molSz,int &solSz){
  int k=0;
  /*****************************************

	   molecule distance binning
 
  ******************************************/
  std::vector<int> molDistNo;
  for(int i=0;i<distBins.size();i++){
    std::pair<double,double> binPr = distBins[i];
    int startK=k;
    //std::cout<<binPr.first<<" "<<binPr.second<<" "<<molDists[k]<<"\n";
    while(molDists[k]>binPr.first &&molDists[k]<=binPr.second){
      //std::cout<<molDists[k]<<"\n";
      k++;
    }
    int endK =k;
    if(i!=0){
      molDistNo.push_back(2.0*(endK-startK));
    }else{
      //if we are considering the first bin it must include the zero distances e,eg, the square of hte number of molecules (why suqe 
      molDistNo.push_back(2.0*(endK-startK)+double(molSz));
    }
    //std::cout<<"mol numbers "<<0.5*(binPr.first+binPr.second)<<" "<<2.0*(endK-startK)<<"\n";
  }
   /*****************************************

	   solvent distance binning
 
  ******************************************/
  std::vector<int> solDistNo;
  k=0;
  for(int i=0;i<distBins.size();i++){
    std::pair<double,double> binPr = distBins[i];
    int startK=k;
    while(solDists[k]>binPr.first && solDists[k]<=binPr.second){
      k++;
    }
    int endK =k;
    if(i!=0){
      solDistNo.push_back(2.0*(endK-startK));
    }else{
      solDistNo.push_back(2.0*(endK-startK)+double(solSz));
    }
    //std::cout<<"sol numbers "<<0.5*(binPr.first+binPr.second)<<" "<<2.0*(endK-startK)<<"\n";
  }
   /*****************************************

	  solvent-molecule distance binning
 
  ******************************************/
  std::vector<int> solMolDistNo;
  k=0;
  for(int i=0;i<distBins.size();i++){
    std::pair<double,double> binPr = distBins[i];
    int startK=k;
    while(solMolDists[k]>binPr.first && solMolDists[k]<=binPr.second){
      k++;
    }
    int endK =k;
    if(i!=0){
      solMolDistNo.push_back(2.0*(endK-startK));
      //std::cout<<binPr.first<<" "<<binPr.second<<" "<<2.0*(endK-startK)<<"\n";
    }else{
      solMolDistNo.push_back(2.0*(endK-startK));
      //std::cout<<binPr.first<<" "<<binPr.second<<" "<<(2.0*(endK-startK)+2.0*double(solSz*molSz))<<"\n";
    }
  }
  /***************************

    now construct the scattering formula
  
  ***********************************/
  std::vector<double> scatVals;
  for(int i=0;i<scatPhases[0].size();i++){
    double scatk=0.0;
    // here we get the scatteirng formuale (decay with q)
    double aminoScat = gaussianFixed(qvals[i]);
    double solScat = gaussianHydFixed(qvals[i]);
    for(int j=0;j<scatPhases.size();j++){
      scatk =  scatk + scatPhases[j][i]*(aminoScat*aminoScat*molDistNo[j] + solScat*solScat*solDistNo[j] + solScat*aminoScat*solMolDistNo[j]);
    }
    scatVals.push_back(scatk);
  }
    // now take the log
  std::vector<double> logdifs;
  double logDifMean=0.0;
  for(int i=0;i<scatVals.size();i++){
    //std::cout<<std::log(scatVals[i])+(std::log(experimentalIntensity[0])-std::log(scatVals[0]))<<" "<<std::log(experimentalIntensity[i])<<"\n";
    double logScatDif= std::log(scatVals[i]) - std::log(experimentalIntensity[i]);
    logdifs.push_back(logScatDif);
    logDifMean =  logDifMean + logScatDif;
  }
  logDifMean = logDifMean/scatVals.size();
  double pred = 0.0;
  for(int i=0;i<scatVals.size();i++){
    double scatDif = logdifs[i] - logDifMean;
    pred = pred + scatDif*scatDif;
  }
  return pred/(scatVals.size()-1);
}





double experimentalData::fitToScatteringMultiple(std::vector<std::vector<double> > &molDists,std::vector<std::vector<double> > &solDists,std::vector<std::vector<double> > &solMolDists,std::vector<int> &molSz,std::vector<int> &solSz,std::vector<double> &percentageCombinations){
  /*****************************************

	   molecule distance binning
 
  ******************************************/
  std::vector<double> molDistNo(distBins.size(),0);
  for(int j=0;j<molDists.size();j++){
    int k=0;
    for(int i=0;i<distBins.size();i++){
      std::pair<double,double> binPr = distBins[i];
      int startK=k;
      //std::cout<<binPr.first<<" "<<binPr.second<<" "<<molDists[k]<<"\n";
      while(molDists[j][k]>binPr.first &&molDists[j][k]<=binPr.second){
	//std::cout<<molDists[k]<<"\n";
	k++;
      }
      int endK =k;
      if(i!=0){
	molDistNo[i]=molDistNo[i] + percentageCombinations[j]*(2.0*double(endK-startK));
      }else{
	//if we are considering the first bin it must include the zero distances e,eg, the square of hte number of molecules (why suqe 
	molDistNo[i]=molDistNo[i]+ percentageCombinations[j]*(2.0*double(endK-startK)+double(molSz[j]));
      }
    }
    //std::cout<<"mol numbers "<<0.5*(binPr.first+binPr.second)<<" "<<2.0*(endK-startK)<<"\n";
  }
   /*****************************************

	   solvent distance binning
 
  ******************************************/
 
  std::vector<double> solDistNo(distBins.size(),0);
  for(int j=0;j<solDists.size();j++){
    int k=0;
    for(int i=0;i<distBins.size();i++){
      std::pair<double,double> binPr = distBins[i];
      int startK=k;
      while(solDists[j][k]>binPr.first && solDists[j][k]<=binPr.second){
	k++;
      }
      int endK =k;
      if(i!=0){
	solDistNo[i]=solDistNo[i] +  percentageCombinations[j]*2.0*(endK-startK);
      }else{
	solDistNo[i]=solDistNo[i] +  percentageCombinations[j]*(2.0*(endK-startK)+double(solSz[j]));
      }
    //std::cout<<"sol numbers "<<0.5*(binPr.first+binPr.second)<<" "<<2.0*(endK-startK)<<"\n";
    }
  }
   /*****************************************

	  solvent-molecule distance binning
 
   ******************************************/
  std::vector<double> solMolDistNo(distBins.size(),0.0);
  for(int j=0;j<solDists.size();j++){
    int k=0;
    for(int i=0;i<distBins.size();i++){
      std::pair<double,double> binPr = distBins[i];
      int startK=k;
      while(solMolDists[j][k]>binPr.first && solMolDists[j][k]<=binPr.second){
	k++;
      }
      int endK =k;
      if(i!=0){
	solMolDistNo[i] = solMolDistNo[i]+ percentageCombinations[j]*2.0*(endK-startK);
	//std::cout<<binPr.first<<" "<<binPr.second<<" "<<2.0*(endK-startK)<<"\n";
      }else{
	solMolDistNo[i] = solMolDistNo[i] +  percentageCombinations[j]*2.0*(endK-startK);
      }
    }
  }
  
  /***************************

    now construct the scattering formula
  
  ***********************************/
  double pred=1000.0;
  for(int hi=0;hi<5;hi++){
    double hydRat = 0.75 + 0.5*hi/4.0;
    std::vector<double> scatVals;  
    for(int i=0;i<scatPhases[0].size();i++){
      double scatk=0.0;
      // here we get the scatteirng formuale (decay with q)
      double aminoScat = gaussianFixed(qvals[i]);
      double solScat = gaussianHydFixed(qvals[i]);
      for(int j=0;j<scatPhases.size();j++){
	//std::cout<<"phases "<<scatPhases[j][i]<<" "<<molDistNo[j]<<" "<<solDistNo[j]<<" "<<solMolDistNo[j]<<"\n";
	scatk =  scatk + scatPhases[j][i]*(aminoScat*aminoScat*molDistNo[j] +hydRat*hydRat*solScat*solScat*solDistNo[j] + hydRat*solScat*aminoScat*solMolDistNo[j]);
      }
      //std::cout<<i<<" "<<scatk<<"\n";
      scatVals.push_back(scatk);
    }
    // now take the log
    std::vector<double> logdifs;
    double logDifMean=0.0;
    int noMean=0;
    for(int i=0;i<scatVals.size();i++){
      //std::cout<<std::log(scatVals[i])<<" "<<std::log(scatVals[i])+(std::log(experimentalIntensity[0])-std::log(scatVals[0]))<<" "<<std::log(experimentalIntensity[i])<<"\n";
      double logScatDif= std::log(scatVals[i]) - std::log(experimentalIntensity[i]);
      logdifs.push_back(logScatDif);
      if(qvals[i]<0.1){
	logDifMean =  logDifMean + logScatDif;
	noMean++;
      }
    }
    logDifMean = logDifMean/double(noMean);
    //logDifMean =  std::log(scatVals[0]) - std::log(experimentalIntensity[0]);
    double predTemp = 0.0;
    for(int i=0;i<scatVals.size();i++){
      double scatDif = logdifs[i] - logDifMean;
      //std::cout<<i<<" "<<experimentalSD[i]<<"\n";
      // predTemp = predTemp + scatDif*scatDif/experimentalSD[i];
      //std::cout<<qvals[i]<<"\n";
      predTemp = predTemp + std::exp(-50000.0*qvals[i]*qvals[i]*qvals[i]*qvals[i]*qvals[i]*qvals[i])*scatDif*scatDif;
      //predTemp = predTemp + scatDif*scatDif;  
    }
    if(predTemp<pred){
      pred = predTemp;
    }
  }
  return pred/(scatPhases[0].size()-1);
}








void experimentalData::writeScatteringToFile(std::vector<double> &molDists,std::vector<double> &solDists,std::vector<double> &solMolDists,int &molSz,int &solSz,const char* filename){
  int k=0;
  /*****************************************

	   molecule distance binning
 
  ******************************************/
  std::vector<int> molDistNo;
  for(int i=0;i<distBins.size();i++){
    std::pair<double,double> binPr = distBins[i];
    int startK=k;
    //std::cout<<binPr.first<<" "<<binPr.second<<"\n";
    while(molDists[k]>binPr.first &&molDists[k]<=binPr.second){
      k++;
    }
    int endK =k;
    if(i!=0){
      molDistNo.push_back(2.0*(endK-startK));
    }else{
      //if we are considering the first bin it must include the zero distances e,eg, the square of hte number of molecules (why suqe 
      molDistNo.push_back(2.0*(endK-startK)+double(molSz));
    }
  }
   /*****************************************

	   solvent distance binning
 
  ******************************************/
  std::vector<int> solDistNo;
  k=0;
  for(int i=0;i<distBins.size();i++){
    std::pair<double,double> binPr = distBins[i];
    int startK=k;
    while(solDists[k]>binPr.first && solDists[k]<=binPr.second){
      k++;
    }
    int endK =k;
    if(i!=0){
      solDistNo.push_back(2.0*(endK-startK));
    }else{
      solDistNo.push_back(2.0*(endK-startK)+double(solSz));
    }
  }
   /*****************************************

	  solvent-molecule distance binning
 
  ******************************************/
  std::vector<int> solMolDistNo;
  k=0;
  for(int i=0;i<distBins.size();i++){
    std::pair<double,double> binPr = distBins[i];
    int startK=k;
    while(solMolDists[k]>binPr.first && solMolDists[k]<=binPr.second){
      k++;
    }
    int endK =k;
    if(i!=0){
      solMolDistNo.push_back(2.0*(endK-startK));
      //std::cout<<binPr.first<<" "<<binPr.second<<" "<<2.0*(endK-startK)<<"\n";
    }else{
      solMolDistNo.push_back(2.0*(endK-startK));
      //std::cout<<binPr.first<<" "<<binPr.second<<" "<<(2.0*(endK-startK)+2.0*double(solSz*molSz))<<"\n";
    }
  }
  /***************************

    now construct the scattering formula
  
  ***********************************/
  std::vector<double> scatVals;
  for(int i=0;i<scatPhases[0].size();i++){
    double scatk=0.0;
    // here we get the scatteirng formuale (decay with q)
    double aminoScat = gaussianFixed(qvals[i]);
    double solScat = gaussianHydFixed(qvals[i]);
    for(int j=0;j<scatPhases.size();j++){
      scatk =  scatk + scatPhases[j][i]*(aminoScat*aminoScat*molDistNo[j] + solScat*solScat*solDistNo[j] + solScat*aminoScat*solMolDistNo[j]);
    }
    scatVals.push_back(scatk);
  }
    // now take the log
  std::vector<std::pair<double,double> > scatterData;
  std::ofstream myfile;
  myfile.open(filename);
  std::vector<double> logdifs;
  double logDifMean=0.0;
  for(int i=0;i<scatVals.size();i++){
    // std::cout<<std::log(scatVals[i])+(std::log(experimentalIntensity[0])-std::log(scatVals[0]))<<" "<<std::log(experimentalIntensity[i])<<"\n";
    double logScatDif= std::log(scatVals[i]) - std::log(experimentalIntensity[i]);
    logdifs.push_back(logScatDif);
    logDifMean =  logDifMean + logScatDif;
  }
  logDifMean = logDifMean/scatVals.size();
  for(int i=0;i<scatVals.size();i++){
    myfile<<qvals[i]<<" "<<scatVals[i]<<" "<<std::log(scatVals[i])- logDifMean<<" "<<std::log(experimentalIntensity[i])<<"\n";
  }
  myfile.close();
}


void experimentalData::writeScatteringToFileMultiple(std::vector<std::vector<double> > &molDists,std::vector<std::vector<double> > &solDists,std::vector<std::vector<double> > &solMolDists,std::vector<int> &molSz,std::vector<int> &solSz,std::vector<double> &percentageCombinations,const char* filename){
  /*****************************************

	   molecule distance binning
 
  ******************************************/
  std::vector<double> molDistNo(distBins.size(),0);
  for(int j=0;j<molDists.size();j++){
    int k=0;
    for(int i=0;i<distBins.size();i++){
      std::pair<double,double> binPr = distBins[i];
      int startK=k;
      //std::cout<<binPr.first<<" "<<binPr.second<<" "<<molDists[k]<<"\n";
      while(molDists[j][k]>binPr.first &&molDists[j][k]<=binPr.second){
	//std::cout<<molDists[k]<<"\n";
	k++;
      }
      int endK =k;
      if(i!=0){
	molDistNo[i]=molDistNo[i] + percentageCombinations[j]*(2.0*double(endK-startK));
      }else{
	//if we are considering the first bin it must include the zero distances e,eg, the square of hte number of molecules (why suqe 
	molDistNo[i]=molDistNo[i]+ percentageCombinations[j]*(2.0*double(endK-startK)+double(molSz[j]));
      }
    }
    //std::cout<<"mol numbers "<<0.5*(binPr.first+binPr.second)<<" "<<2.0*(endK-startK)<<"\n";
  }
   /*****************************************

	   solvent distance binning
 
  ******************************************/
 
  std::vector<double> solDistNo(distBins.size(),0);
  for(int j=0;j<solDists.size();j++){
    int k=0;
    for(int i=0;i<distBins.size();i++){
      std::pair<double,double> binPr = distBins[i];
      int startK=k;
      while(solDists[j][k]>binPr.first && solDists[j][k]<=binPr.second){
	k++;
      }
      int endK =k;
      if(i!=0){
	solDistNo[i]=solDistNo[i] +  percentageCombinations[j]*2.0*(endK-startK);
      }else{
	solDistNo[i]=solDistNo[i] +  percentageCombinations[j]*(2.0*(endK-startK)+double(solSz[j]));
      }
    //std::cout<<"sol numbers "<<0.5*(binPr.first+binPr.second)<<" "<<2.0*(endK-startK)<<"\n";
    }
  }
   /*****************************************

	  solvent-molecule distance binning
 
   ******************************************/
  std::vector<double> solMolDistNo(distBins.size(),0.0);
  for(int j=0;j<solDists.size();j++){
    int k=0;
    for(int i=0;i<distBins.size();i++){
      std::pair<double,double> binPr = distBins[i];
      int startK=k;
      while(solMolDists[j][k]>binPr.first && solMolDists[j][k]<=binPr.second){
	k++;
      }
      int endK =k;
      if(i!=0){
	solMolDistNo[i] = solMolDistNo[i]+ percentageCombinations[j]*2.0*(endK-startK);
	//std::cout<<binPr.first<<" "<<binPr.second<<" "<<2.0*(endK-startK)<<"\n";
      }else{
	solMolDistNo[i] = solMolDistNo[i] +  percentageCombinations[j]*2.0*(endK-startK)+double(solSz[j]);
	//std::cout<<binPr.first<<" "<<binPr.second<<" "<<(2.0*(endK-startK)+2.0*double(solSz*molSz))<<"\n";
      }
    }
  }
  
   /***************************

    now construct the scattering formula
  
  ***********************************/
  double pred=1000.0;
  double hydRatOut =1.0;
  for(int hi=0;hi<5;hi++){
    double hydRat = 0.75 + 0.5*hi/4.0;
    std::vector<double> scatVals;  
    for(int i=0;i<scatPhases[0].size();i++){
      double scatk=0.0;
      // here we get the scatteirng formuale (decay with q)
      double aminoScat = gaussianFixed(qvals[i]);
      double solScat = gaussianHydFixed(qvals[i]);
      for(int j=0;j<scatPhases.size();j++){
	//  std::cout<<"phases "<<scatPhases[j][i]<<" "<<molDistNo[j]<<" "<<solDistNo[j]<<" "<<solMolDistNo[j]<<"\n";
	scatk =  scatk + scatPhases[j][i]*(aminoScat*aminoScat*molDistNo[j] +hydRat*hydRat*solScat*solScat*solDistNo[j] + hydRat*solScat*aminoScat*solMolDistNo[j]);
      }
      //std::cout<<i<<" "<<scatk<<"\n";
      scatVals.push_back(scatk);
    }
    // now take the log
    std::vector<double> logdifs;
    double logDifMean=0.0;
    int noMean=0;
    for(int i=0;i<scatVals.size();i++){
      //std::cout<<std::log(scatVals[i])<<" "<<std::log(scatVals[i])+(std::log(experimentalIntensity[0])-std::log(scatVals[0]))<<" "<<std::log(experimentalIntensity[i])<<"\n";
      double logScatDif= std::log(scatVals[i]) - std::log(experimentalIntensity[i]);
      logdifs.push_back(logScatDif);
      if(qvals[i]<0.1){
	logDifMean =  logDifMean + logScatDif;
	noMean++;
      }
    }
    logDifMean = logDifMean/double(noMean);
    //logDifMean =  std::log(scatVals[0]) - std::log(experimentalIntensity[0]);
    double predTemp = 0.0;
    for(int i=0;i<scatVals.size();i++){
      double scatDif = logdifs[i] - logDifMean;
      predTemp = predTemp + scatDif*scatDif/experimentalSD[i];
    }
    if(predTemp<pred){
      pred = predTemp;
      hydRatOut = hydRat;
    }
  }
  std::vector<double> scatVals;
  for(int i=0;i<scatPhases[0].size();i++){
    double scatk=0.0;
    // here we get the scatteirng formuale (decay with q)
    double aminoScat = gaussianFixed(qvals[i]);
    double solScat = gaussianHydFixed(qvals[i]);
    for(int j=0;j<scatPhases.size();j++){
      scatk =  scatk + scatPhases[j][i]*(aminoScat*aminoScat*molDistNo[j] +hydRatOut*hydRatOut*solScat*solScat*solDistNo[j] + hydRatOut*solScat*aminoScat*solMolDistNo[j]);
    }
    scatVals.push_back(scatk);
  }
    // now take the log
  std::vector<std::pair<double,double> > scatterData;
  std::ofstream myfile;
  myfile.open(filename);
  std::vector<double> logdifs;
  double logDifMean=0.0;
  int noMean=0;
  for(int i=0;i<scatVals.size();i++){
    //std::cout<<std::log(scatVals[i])+(std::log(experimentalIntensity[0])-std::log(scatVals[0]))<<" "<<std::log(experimentalIntensity[i])<<"\n";
    double logScatDif= std::log(scatVals[i]) - std::log(experimentalIntensity[i]);
    logdifs.push_back(logScatDif);
    if(qvals[i]<0.1){
      logDifMean =  logDifMean + logScatDif;
      noMean++;
    }
  }
  logDifMean = logDifMean/double(noMean);
  //logDifMean =  std::log(scatVals[0]) - std::log(experimentalIntensity[0]);
  pred = 0.0;
  for(int i=0;i<scatVals.size();i++){
    myfile<<qvals[i]<<" "<<scatVals[i]<<" "<<std::log(scatVals[i])- logDifMean<<" "<<std::log(experimentalIntensity[i])<<"\n";
  }
  // finally add the percentage combination values
  for(int i=0;i<percentageCombinations.size();i++){
    if(i==(percentageCombinations.size()-1)){
      myfile<<percentageCombinations[i]<<"\n";
    }else{
      myfile<<percentageCombinations[i]<<" ";
    }
  }
  myfile.close();
}




void experimentalData::writeRawMolScatteringToFileMultiple(std::vector<std::vector<double> > &molDists,std::vector<std::vector<double> > &solDists,std::vector<std::vector<double> > &solMolDists,std::vector<int> &molSz,std::vector<int> &solSz,std::vector<double> &percentageCombinations,const char* filename){
  /*****************************************

	   molecule distance binning
 
  ******************************************/
  std::vector<double> molDistNo(distBins.size(),0);
  for(int j=0;j<molDists.size();j++){
    int k=0;
    for(int i=0;i<distBins.size();i++){
      std::pair<double,double> binPr = distBins[i];
      int startK=k;
      //std::cout<<binPr.first<<" "<<binPr.second<<" "<<molDists[k]<<"\n";
      while(molDists[j][k]>binPr.first &&molDists[j][k]<=binPr.second){
	//std::cout<<molDists[k]<<"\n";
	k++;
      }
      int endK =k;
      if(i!=0){
	molDistNo[i]=molDistNo[i] + percentageCombinations[j]*(2.0*double(endK-startK));
      }else{
	//if we are considering the first bin it must include the zero distances e,eg, the square of hte number of molecules (why suqe 
	molDistNo[i]=molDistNo[i]+ percentageCombinations[j]*(2.0*double(endK-startK)+double(molSz[j]));
      }
    }
    //std::cout<<"mol numbers "<<0.5*(binPr.first+binPr.second)<<" "<<2.0*(endK-startK)<<"\n";
  }
  std::vector<double> scatVals;
  for(int i=0;i<scatPhases[0].size();i++){
      double scatk=0.0;
      // here we get the scatteirng formuale (decay with q)
      double aminoScat = gaussianFixed(qvals[i]);
      double solScat = gaussianHydFixed(qvals[i]);
      for(int j=0;j<scatPhases.size();j++){
	//  std::cout<<"phases "<<scatPhases[j][i]<<" "<<molDistNo[j]<<" "<<solDistNo[j]<<" "<<solMolDistNo[j]<<"\n";
	scatk =  scatk + scatPhases[j][i]*(aminoScat*aminoScat*molDistNo[j]);
      }
      //std::cout<<i<<" "<<scatk<<"\n";
      scatVals.push_back(scatk);
    }
    // now take the log
    std::vector<double> logdifs;
    double logDifMean=0.0;
    int noMean=0;
    for(int i=0;i<scatVals.size();i++){
      //std::cout<<std::log(scatVals[i])<<" "<<std::log(scatVals[i])+(std::log(experimentalIntensity[0])-std::log(scatVals[0]))<<" "<<std::log(experimentalIntensity[i])<<"\n";
      double logScatDif= std::log(scatVals[i]) - std::log(experimentalIntensity[i]);
      logdifs.push_back(logScatDif);
      if(qvals[i]<0.1){
	logDifMean =  logDifMean + logScatDif;
	noMean++;
      }
    }
    // now take the log
  std::vector<std::pair<double,double> > scatterData;
  std::ofstream myfile;
  myfile.open(filename);
  for(int i=0;i<scatVals.size();i++){
    myfile<<qvals[i]<<" "<<std::log(scatVals[i])<<"\n";
  }
  // finally add the percentage combination values
  for(int i=0;i<percentageCombinations.size();i++){
    if(i==(percentageCombinations.size()-1)){
      myfile<<percentageCombinations[i]<<"\n";
    }else{
      myfile<<percentageCombinations[i]<<" ";
    }
  }
  myfile.close();
}
