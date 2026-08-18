#include "randomMolGen.h"

randomMol::randomMol(double &rminIn,double &rmaxIn,double &lminIn){
  rmin =rminIn;rmax = rmaxIn;
  std::ifstream myfile;
  myfile.open("probabilityInterpolation/mixedLinkerYfromXCumDist.dat");
  double val,x,y,den;
  std::string output;
  int xindex=0;
  int n=0;
  nGridSize=100;
  /*point tan1(10000.0,10000.0,10000.0);point norm1(10000.0,10000.0,10000.0);
  point norm2(10000.0,10000.0,10000.0);
  dummyFrame.push_back(tan1);dummyFrame.push_back(norm1);dummyFrame.push_back(norm2);*/
  std::vector<double> subVecMixedYFX(nGridSize,0.0);
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
      ss>>val;  
      x=val;
      ss>>val;
      y=val;
      ss>>val;
      den=val;
      subVecMixedYFX[n]=val;
      if(n==nGridSize-1){
        n=-1;
        xindex++;
        mixedYfromXInterpolants.push_back(subVecMixedYFX);
      };
      n++;
    }
     mixedYfromXInterpolants.push_back(subVecMixedYFX);
  }
  myfile.close();
  myfile.open("probabilityInterpolation/negBetaLinkerYfromXCumDist.dat");
  std::vector<double> subVecNegBetaYFX(nGridSize,0.0);
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      y=val;
      ss>>val;
      den=val;
      subVecNegBetaYFX[n]=val;
       if(n==nGridSize-1){
        n=-1;
        xindex++;
        negBetaYfromXInterpolants.push_back(subVecNegBetaYFX);
      };
      n++;
    }
     negBetaYfromXInterpolants.push_back(subVecNegBetaYFX);
  }
  myfile.close();
  myfile.open("probabilityInterpolation/posBetaLinkerYfromXCumDist.dat");
  std::vector<double> subVecPosBetaYFX(nGridSize,0.0);
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      y=val;
      ss>>val;
      den=val;
      subVecPosBetaYFX[n]=val;
       if(n==nGridSize-1){
        n=-1;
        xindex++;
        posBetaYfromXInterpolants.push_back(subVecPosBetaYFX);
      };
      n++;
    }
     posBetaYfromXInterpolants.push_back(subVecPosBetaYFX);
  }
  myfile.close();
  myfile.open("probabilityInterpolation/alphaLinkerYfromXCumDist.dat");
  std::vector<double> subVecAlphaYFX(nGridSize,0.0);
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      y=val;
      ss>>val;
      den=val;
      subVecAlphaYFX[n]=val;
      if(n==nGridSize-1){
        n=-1;
        xindex++;
        alphaYfromXInterpolants.push_back(subVecAlphaYFX);
      };
      n++;
    }
     alphaYfromXInterpolants.push_back(subVecAlphaYFX);
  }
  myfile.close();
  myfile.open("probabilityInterpolation/alphaStrandYfromXCumDist.dat");
  std::vector<double> subVecAlphaStrandYFX(nGridSize,0.0);
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      y=val;
      ss>>val;
      den=val;
      subVecAlphaStrandYFX[n]=val;
      if(n==nGridSize-1){
        n=-1;
        xindex++;
        alphaStrandYfromXInterpolants.push_back(subVecAlphaStrandYFX);
      };
      n++;
    }
     alphaStrandYfromXInterpolants.push_back(subVecAlphaStrandYFX);
  }
  myfile.close();
  myfile.open("probabilityInterpolation/posBetaStrandYfromXCumDist.dat");
  std::vector<double> subVecPosBetaStrandYFX(nGridSize,0.0);
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      y=val;
      ss>>val;
      den=val;
      subVecPosBetaStrandYFX[n]=val;
      if(n==nGridSize-1){
	n=-1;
        xindex++;
        posBetaStrandYfromXInterpolants.push_back(subVecPosBetaStrandYFX);
      };
      n++;
    }
     posBetaStrandYfromXInterpolants.push_back(subVecPosBetaStrandYFX);
  }
  myfile.close();
  myfile.open("probabilityInterpolation/negBetaStrandYfromXCumDist.dat");
  std::vector<double> subVecNegBetaStrandYFX(nGridSize,0.0);
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      y=val;
      ss>>val;
      den=val;
      subVecNegBetaStrandYFX[n]=val;
      if(n==nGridSize-1){
        n=-1;
        xindex++;
        negBetaStrandYfromXInterpolants.push_back(subVecNegBetaStrandYFX);
      };
      n++;
    }
     negBetaStrandYfromXInterpolants.push_back(subVecNegBetaStrandYFX);
  }
  myfile.close();

  myfile.open("probabilityInterpolation/mixedStrandYfromXCumDist.dat");
  std::vector<double> subVecMixedStrandYFX(nGridSize,0.0);
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      y=val;
      ss>>val;
      den=val;
      subVecMixedStrandYFX[n]=val;
      if(n==nGridSize-1){
        n=-1;
        xindex++;
        mixedStrandYfromXInterpolants.push_back(subVecMixedStrandYFX);
      };
      n++;
    }
     mixedStrandYfromXInterpolants.push_back(subVecMixedStrandYFX);
  }
  myfile.close();
  myfile.open("probabilityInterpolation/joinHelixToLinkerYfromXCumDist.dat");
  std::vector<double> subVecJHTOLYFX(nGridSize,0.0);
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      y=val;
      ss>>val;
      den=val;
      subVecJHTOLYFX[n]=val;
      if(n==nGridSize-1){
        n=-1;
        xindex++;
        helixToLinkerYfromXInterpolants.push_back(subVecJHTOLYFX);
      };
      n++;
    }
     helixToLinkerYfromXInterpolants.push_back(subVecJHTOLYFX);
  }
  myfile.close();
  myfile.open("probabilityInterpolation/joinHelixToStrandYfromXCumDist.dat");
  std::vector<double> subVecJHTOLStrandYFX(nGridSize,0.0);
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      y=val;
      ss>>val;
      den=val;
      subVecJHTOLStrandYFX[n]=val;
      if(n==nGridSize-1){
        n=-1;
        xindex++;
        helixToStrandYfromXInterpolants.push_back(subVecJHTOLStrandYFX);
      };
      n++;
    }
     helixToStrandYfromXInterpolants.push_back(subVecJHTOLStrandYFX);
  }
  myfile.close();
  myfile.open("probabilityInterpolation/joinLinkerToHelixYfromXCumDist.dat");
  std::vector<double> subVecJLTOHYFX(nGridSize,0.0);
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      y=val;
      ss>>val;
      den=val;
      subVecJLTOHYFX[n]=val;
      if(n==nGridSize-1){
        n=-1;
        xindex++;
        linkerToHelixYfromXInterpolants.push_back(subVecJLTOHYFX);
      };
      n++;
    }
     linkerToHelixYfromXInterpolants.push_back(subVecJLTOHYFX);
  }
  myfile.close();
  myfile.open("probabilityInterpolation/joinStrandToHelixYfromXCumDist.dat");
  std::vector<double> subVecJSTOHYFX(nGridSize,0.0);
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      y=val;
      ss>>val;
      den=val;
      subVecJSTOHYFX[n]=val;
      if(n==nGridSize-1){
        n=-1;
        xindex++;
        strandToHelixYfromXInterpolants.push_back(subVecJSTOHYFX);
      };
      n++;
    }
     strandToHelixYfromXInterpolants.push_back(subVecJSTOHYFX);
  }
  myfile.close();
   myfile.open("probabilityInterpolation/joinStrandToLinkerYfromXCumDist.dat");
  std::vector<double> subVecJSTOLYFX(nGridSize,0.0);
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      y=val;
      ss>>val;
      den=val;
      subVecJSTOLYFX[n]=val;
      if(n==nGridSize-1){
        n=-1;
        xindex++;
        strandToLinkerYfromXInterpolants.push_back(subVecJSTOLYFX);
      };
      n++;
    }
     strandToLinkerYfromXInterpolants.push_back(subVecJSTOLYFX);
  }
  myfile.close();
  myfile.open("probabilityInterpolation/joinLinkerToStrandYfromXCumDist.dat");
  std::vector<double> subVecJLTOSYFX(nGridSize,0.0);
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      y=val;
      ss>>val;
      den=val;
      subVecJLTOSYFX[n]=val;
      if(n==nGridSize-1){
        n=-1;
        xindex++;
        linkerToStrandYfromXInterpolants.push_back(subVecJLTOSYFX);
      };
      n++;
    }
     linkerToStrandYfromXInterpolants.push_back(subVecJLTOSYFX);
  }
  myfile.close();
  myfile.open("probabilityInterpolation/mixedLinkerxCumDist.dat");
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      den=val;
      mixedXCumInterpolants.push_back(val);
    }
  }else{
    std::cout<<"failed to open mixedLinkerxCumDist.dat\n";
  }
  myfile.close();
  myfile.open("probabilityInterpolation/mixedStrandxCumDist.dat");
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      den=val;
      mixedStrandXCumInterpolants.push_back(val);
    }
  }else{
    std::cout<<"failed to open mixedStrandxCumDist.dat\n";
  }
  myfile.close();
  myfile.open("probabilityInterpolation/negBetaLinkerxCumDist.dat");
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      den=val;
      negBetaXCumInterpolants.push_back(val);
    }
  }else{
    std::cout<<"failed to open negBetaLinkerxCumDist.dat\n";
  }
  myfile.close();
   myfile.open("probabilityInterpolation/negBetaStrandxCumDist.dat");
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      den=val;
      negBetaStrandXCumInterpolants.push_back(val);
    }
  }else{
    std::cout<<"failed to open negBetaLinkerxCumDist.dat\n";
  }
  myfile.close();
  myfile.open("probabilityInterpolation/posBetaLinkerxCumDist.dat");
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      den=val;
      posBetaXCumInterpolants.push_back(val);
    }
  }else{
    std::cout<<"failed to open posBetaLinkerxCumDist.dat\n";
  }
  myfile.close();
  myfile.open("probabilityInterpolation/posBetaStrandxCumDist.dat");
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      den=val;
      posBetaStrandXCumInterpolants.push_back(val);
    }
  }else{
    std::cout<<"failed to open posBetaStrandxCumDist.dat\n";
  }
  myfile.close();
  myfile.open("probabilityInterpolation/alphaLinkerxCumDist.dat");
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      den=val;
      alphaXCumInterpolants.push_back(val);
    }
  }else{
    std::cout<<"failed to open alphaLinkerxCumDist.dat\n";
  }
  myfile.close();
  myfile.open("probabilityInterpolation/alphaStrandxCumDist.dat");
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      den=val;
      alphaStrandXCumInterpolants.push_back(val);
    }
  }else{
    std::cout<<"failed to open alphaStrandxCumDist.dat\n";
  }
  myfile.close();
  myfile.open("probabilityInterpolation/joinLinkerToHelixxCumDist.dat");
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      den=val;
      linkerToHelixXCumInterpolants.push_back(val);
    }
  }else{
    std::cout<<"failed to open joinLinkerToHelixxCumDist.dat\n";
  }
  myfile.close();
  myfile.open("probabilityInterpolation/joinStrandToHelixxCumDist.dat");
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      den=val;
      strandToHelixXCumInterpolants.push_back(val);
    }
  }else{
    std::cout<<"failed to open joinStrandToHelixxCumDist.dat\n";
  }
  myfile.close();
  myfile.open("probabilityInterpolation/joinHelixToLinkerxCumDist.dat");
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      den=val;
      helixToLinkerXCumInterpolants.push_back(val);
    }
  }else{
    std::cout<<"failed to open joinHelixToLinkerxCumDist.dat\n";
  }
  myfile.close();
  myfile.open("probabilityInterpolation/joinHelixToStrandxCumDist.dat");
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      den=val;
      helixToStrandXCumInterpolants.push_back(val);
    }
  }else{
    std::cout<<"failed to open joinHelixToStrandxCumDist.dat\n";
  }
  myfile.close();
  myfile.open("probabilityInterpolation/joinStrandToLinkerxCumDist.dat");
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      den=val;
      strandToLinkerXCumInterpolants.push_back(val);
    }
  }else{
    std::cout<<"failed to open joinHelixToLinkerxCumDist.dat\n";
  }
  myfile.close();
  myfile.open("probabilityInterpolation/joinLinkerToStrandxCumDist.dat");
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      den=val;
      linkerToStrandXCumInterpolants.push_back(val);
    }
  }else{
    std::cout<<"failed to open joinHelixToStrandxCumDist.dat\n";
  }
  myfile.close();
  // set up random number generator
  std::uniform_real_distribution<double> distribution(0.0,1.0);
  // set the domain bounds
  kapL=0.0;kapU=0.7;tauL=-0.7;tauU=0.7;
  nkap=nGridSize;ntau=nGridSize;dkap = 0.7;dtau = 1.4;
  pAlpha=0.8794108622276772;pposBeta=0.002247621969929464;pnegBeta=0.0599877262940994;pMixed=1.0-(pAlpha+pposBeta+pnegBeta);
  pAlphaStrand =0.270498;pposBetaStrand=0.0442899;pnegBetaStrand=0.449854;pMixedStrand=0.235359;
  lmin = lminIn;
  alphaKap = 0.3812;alphaTau=0.1465;alphaAl=4.24153;alphaCl=3.8003802317894965;
  betaKap = 0.523;betaTau=0.5458;betaAl=4.58415;betaCl=3.770814833345433;
}

void randomMol::setParams(double &rminIn,double &rmaxIn,double &lminIn){
  rmin =rminIn;rmax = rmaxIn;
   std::ifstream myfile;
  myfile.open("probabilityInterpolation/mixedLinkerYfromXCumDist.dat");
  double val,x,y,den;
  std::string output;
  int xindex=0;
  int n=0;
  nGridSize = 100;
  std::vector<double> subVecMixedYFX(nGridSize,0.0);
  /*point tan1(10000.0,10000.0,10000.0);point norm1(10000.0,10000.0,10000.0);
  point norm2(10000.0,10000.0,10000.0);
  dummyFrame.push_back(tan1);dummyFrame.push_back(norm1);dummyFrame.push_back(norm2);*/
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      y=val;
      ss>>val;
      den=val;
      subVecMixedYFX[n]=val;
      if(n==nGridSize-1){
        n=-1;
        xindex++;
        mixedYfromXInterpolants.push_back(subVecMixedYFX);
      };
      n++;
    }
     mixedYfromXInterpolants.push_back(subVecMixedYFX);
  }
  myfile.close();
  myfile.open("probabilityInterpolation/negBetaLinkerYfromXCumDist.dat");
  std::vector<double> subVecNegBetaYFX(nGridSize,0.0);
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      y=val;
      ss>>val;
      den=val;
      subVecNegBetaYFX[n]=val;
       if(n==nGridSize-1){
        n=-1;
        xindex++;
        negBetaYfromXInterpolants.push_back(subVecNegBetaYFX);
      };
      n++;
    }
     negBetaYfromXInterpolants.push_back(subVecNegBetaYFX);
  }
  myfile.close();
  myfile.open("probabilityInterpolation/posBetaLinkerYfromXCumDist.dat");
  std::vector<double> subVecPosBetaYFX(nGridSize,0.0);
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      y=val;
      ss>>val;
      den=val;
      subVecPosBetaYFX[n]=val;
       if(n==nGridSize-1){
        n=-1;
        xindex++;
        posBetaYfromXInterpolants.push_back(subVecPosBetaYFX);
      };
      n++;
    }
     posBetaYfromXInterpolants.push_back(subVecPosBetaYFX);
  }
  myfile.close();
  myfile.open("probabilityInterpolation/alphaLinkerYfromXCumDist.dat");
  std::vector<double> subVecAlphaYFX(nGridSize,0.0);
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      y=val;
      ss>>val;
      den=val;
      subVecAlphaYFX[n]=val;
      if(n==nGridSize-1){
        n=-1;
        xindex++;
        alphaYfromXInterpolants.push_back(subVecAlphaYFX);
      };
      n++;
    }
     alphaYfromXInterpolants.push_back(subVecAlphaYFX);
  }
  myfile.close();
  myfile.open("probabilityInterpolation/alphaStrandYfromXCumDist.dat");
  std::vector<double> subVecAlphaStrandYFX(nGridSize,0.0);
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      y=val;
      ss>>val;
      den=val;
      subVecAlphaStrandYFX[n]=val;
      if(n==nGridSize-1){
        n=-1;
        xindex++;
        alphaStrandYfromXInterpolants.push_back(subVecAlphaStrandYFX);
      };
      n++;
    }
     alphaStrandYfromXInterpolants.push_back(subVecAlphaStrandYFX);
  }
  myfile.close();
  myfile.open("probabilityInterpolation/posBetaStrandYfromXCumDist.dat");
  std::vector<double> subVecPosBetaStrandYFX(nGridSize,0.0);
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      y=val;
      ss>>val;
      den=val;
      subVecPosBetaStrandYFX[n]=val;
      if(n==nGridSize-1){
	n=-1;
        xindex++;
        posBetaStrandYfromXInterpolants.push_back(subVecPosBetaStrandYFX);
      };
      n++;
    }
     posBetaStrandYfromXInterpolants.push_back(subVecPosBetaStrandYFX);
  }
  myfile.close();
  myfile.open("probabilityInterpolation/negBetaStrandYfromXCumDist.dat");
  std::vector<double> subVecNegBetaStrandYFX(nGridSize,0.0);
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      y=val;
      ss>>val;
      den=val;
      subVecNegBetaStrandYFX[n]=val;
      if(n==nGridSize-1){
        n=-1;
        xindex++;
        negBetaStrandYfromXInterpolants.push_back(subVecNegBetaStrandYFX);
      };
      n++;
    }
     negBetaStrandYfromXInterpolants.push_back(subVecNegBetaStrandYFX);
  }
  myfile.close();

  myfile.open("probabilityInterpolation/mixedStrandYfromXCumDist.dat");
  std::vector<double> subVecMixedStrandYFX(nGridSize,0.0);
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      y=val;
      ss>>val;
      den=val;
      subVecMixedStrandYFX[n]=val;
      if(n==nGridSize-1){
        n=-1;
        xindex++;
        mixedStrandYfromXInterpolants.push_back(subVecMixedStrandYFX);
      };
      n++;
    }
     mixedStrandYfromXInterpolants.push_back(subVecMixedStrandYFX);
  }
  myfile.close();
  myfile.open("probabilityInterpolation/joinHelixToLinkerYfromXCumDist.dat");
  std::vector<double> subVecJHTOLYFX(nGridSize,0.0);
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      y=val;
      ss>>val;
      den=val;
      subVecJHTOLYFX[n]=val;
      if(n==nGridSize-1){
        n=-1;
        xindex++;
        helixToLinkerYfromXInterpolants.push_back(subVecJHTOLYFX);
      };
      n++;
    }
     helixToLinkerYfromXInterpolants.push_back(subVecJHTOLYFX);
  }
  myfile.close();
  myfile.open("probabilityInterpolation/joinHelixToStrandYfromXCumDist.dat");
  std::vector<double> subVecJHTOLStrandYFX(nGridSize,0.0);
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      y=val;
      ss>>val;
      den=val;
      subVecJHTOLStrandYFX[n]=val;
      if(n==nGridSize-1){
        n=-1;
        xindex++;
        helixToStrandYfromXInterpolants.push_back(subVecJHTOLStrandYFX);
      };
      n++;
    }
     helixToStrandYfromXInterpolants.push_back(subVecJHTOLStrandYFX);
  }
  myfile.close();
  myfile.open("probabilityInterpolation/joinLinkerToHelixYfromXCumDist.dat");
  std::vector<double> subVecJLTOHYFX(nGridSize,0.0);
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      y=val;
      ss>>val;
      den=val;
      subVecJLTOHYFX[n]=val;
      if(n==nGridSize-1){
        n=-1;
        xindex++;
        linkerToHelixYfromXInterpolants.push_back(subVecJLTOHYFX);
      };
      n++;
    }
     linkerToHelixYfromXInterpolants.push_back(subVecJLTOHYFX);
  }
  myfile.close();
  myfile.open("probabilityInterpolation/joinStrandToHelixYfromXCumDist.dat");
  std::vector<double> subVecJSTOHYFX(nGridSize,0.0);
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      y=val;
      ss>>val;
      den=val;
      subVecJSTOHYFX[n]=val;
      if(n==nGridSize-1){
        n=-1;
        xindex++;
        strandToHelixYfromXInterpolants.push_back(subVecJSTOHYFX);
      };
      n++;
    }
     strandToHelixYfromXInterpolants.push_back(subVecJSTOHYFX);
  }
  myfile.close();
   myfile.open("probabilityInterpolation/joinStrandToLinkerYfromXCumDist.dat");
  std::vector<double> subVecJSTOLYFX(nGridSize,0.0);
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      y=val;
      ss>>val;
      den=val;
      subVecJSTOLYFX[n]=val;
      if(n==nGridSize-1){
        n=-1;
        xindex++;
        strandToLinkerYfromXInterpolants.push_back(subVecJSTOLYFX);
      };
      n++;
    }
     strandToLinkerYfromXInterpolants.push_back(subVecJSTOLYFX);
  }
  myfile.close();
  myfile.open("probabilityInterpolation/joinLinkerToStrandYfromXCumDist.dat");
  std::vector<double> subVecJLTOSYFX(nGridSize,0.0);
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      y=val;
      ss>>val;
      den=val;
      subVecJLTOSYFX[n]=val;
      if(n==nGridSize-1){
        n=-1;
        xindex++;
        linkerToStrandYfromXInterpolants.push_back(subVecJLTOSYFX);
      };
      n++;
    }
     linkerToStrandYfromXInterpolants.push_back(subVecJLTOSYFX);
  }
  myfile.close();
  myfile.open("probabilityInterpolation/mixedLinkerxCumDist.dat");
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
      ss>>val;  
      x=val;
      ss>>val;
      den=val;
      mixedXCumInterpolants.push_back(val);
    }
  }else{
    std::cout<<"failed to open mixedLinkerxCumDist.dat\n";
  }
  myfile.close();
  myfile.open("probabilityInterpolation/mixedStrandxCumDist.dat");
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      den=val;
      mixedStrandXCumInterpolants.push_back(val);
    }
  }else{
    std::cout<<"failed to open mixedStrandxCumDist.dat\n";
  }
  myfile.close();
  myfile.open("probabilityInterpolation/negBetaLinkerxCumDist.dat");
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      den=val;
      negBetaXCumInterpolants.push_back(val);
    }
  }else{
    std::cout<<"failed to open negBetaLinkerxCumDist.dat\n";
  }
  myfile.close();
   myfile.open("probabilityInterpolation/negBetaStrandxCumDist.dat");
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      den=val;
      negBetaStrandXCumInterpolants.push_back(val);
    }
  }else{
    std::cout<<"failed to open negBetaLinkerxCumDist.dat\n";
  }
  myfile.close();
  myfile.open("probabilityInterpolation/posBetaLinkerxCumDist.dat");
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      den=val;
      posBetaXCumInterpolants.push_back(val);
    }
  }else{
    std::cout<<"failed to open posBetaLinkerxCumDist.dat\n";
  }
  myfile.close();
  myfile.open("probabilityInterpolation/posBetaStrandxCumDist.dat");
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      den=val;
      posBetaStrandXCumInterpolants.push_back(val);
    }
  }else{
    std::cout<<"failed to open posBetaStrandxCumDist.dat\n";
  }
  myfile.close();
  myfile.open("probabilityInterpolation/alphaLinkerxCumDist.dat");
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      den=val;
      alphaXCumInterpolants.push_back(val);
    }
  }else{
    std::cout<<"failed to open alphaLinkerxCumDist.dat\n";
  }
  myfile.close();
  myfile.open("probabilityInterpolation/alphaStrandxCumDist.dat");
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      den=val;
      alphaStrandXCumInterpolants.push_back(val);
    }
  }else{
    std::cout<<"failed to open alphaStrandxCumDist.dat\n";
  }
  myfile.close();
  myfile.open("probabilityInterpolation/joinLinkerToHelixxCumDist.dat");
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      den=val;
      linkerToHelixXCumInterpolants.push_back(val);
    }
  }else{
    std::cout<<"failed to open joinLinkerToHelixxCumDist.dat\n";
  }
  myfile.close();
  myfile.open("probabilityInterpolation/joinStrandToHelixxCumDist.dat");
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      den=val;
      strandToHelixXCumInterpolants.push_back(val);
    }
  }else{
    std::cout<<"failed to open joinStrandToHelixxCumDist.dat\n";
  }
  myfile.close();
  myfile.open("probabilityInterpolation/joinHelixToLinkerxCumDist.dat");
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      den=val;
      helixToLinkerXCumInterpolants.push_back(val);
    }
  }else{
    std::cout<<"failed to open joinHelixToLinkerxCumDist.dat\n";
  }
  myfile.close();
   myfile.open("probabilityInterpolation/joinHelixToStrandxCumDist.dat");
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      den=val;
      helixToStrandXCumInterpolants.push_back(val);
    }
  }else{
    std::cout<<"failed to open joinHelixToStrandxCumDist.dat\n";
  }
  myfile.close();
  myfile.open("probabilityInterpolation/joinStrandToLinkerxCumDist.dat");
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      den=val;
      strandToLinkerXCumInterpolants.push_back(val);
    }
  }else{
    std::cout<<"failed to open joinStrandToLinkerxCumDist.dat\n";
  }
  myfile.close();
  myfile.open("probabilityInterpolation/joinLinkerToStrandxCumDist.dat");
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      std::getline(myfile,output);
      std::stringstream ss(output);
	    ss>>val;  
      x=val;
      ss>>val;
      den=val;
      linkerToStrandXCumInterpolants.push_back(val);
    }
  }else{
    std::cout<<"failed to open joinlinkerToStrandxCumDist.dat\n";
  }
  myfile.close();
  // set up random number generator
  std::uniform_real_distribution<double> distribution(0.0,1.0);
  // set the domain bounds
  kapL=0.0;kapU=0.7;tauL=-0.7;tauU=0.7;
  nkap=nGridSize;ntau=nGridSize;dkap = 0.7;dtau = 1.4;
  pAlpha=0.8794108622276772;pposBeta=0.002247621969929464;pnegBeta=0.0599877262940994;pMixed=1.0-(pAlpha+pposBeta+pnegBeta);
  pAlphaStrand =0.270498;pposBetaStrand=0.0442899;pnegBetaStrand=0.449854;pMixedStrand=0.235359;
  lmin = lminIn;
  alphaKap = 0.3812;alphaTau=0.1465;alphaAl=4.24153;alphaCl=3.8003802317894965;
  betaKap = 0.523;betaTau=0.5458;betaAl=4.58415;betaCl=3.770814833345433;
}



std::pair<double,double> randomMol::kapTau(point &pt1,point &pt2,point &pt3,point &pt4){
  std::pair<double,double> kapTau;
  point pt1Mean = (pt1+pt2)*0.5;point pt2Mean = (pt2+pt3)*0.5;point pt3Mean = (pt3+pt4)*0.5;
  point pdif1 = pt1Mean-pt3Mean;point pdif2 = pt2Mean-pt3Mean;
  double dot =pdif1.dotprod(pdif2);
  double len1= pdif1.length();double len2=pdif2.length();
  double va = std::acos(dot/len1/len2);
  double area = 0.5*std::sin(va);
  point pdif3= pt2Mean-pt1Mean;
  double plen = pdif3.length();
  kapTau.first = 4.0*area/plen;
  point e1 = pt2-pt1;point e2= pt3-pt2;point e3 = pt4-pt3;
  double l1 = 0.5*(e2.length()+e3.length());
  e1.normalise();e2.normalise();e3.normalise();
  point b1 = binorm(e1,e2);point b2 = binorm(e2,e3);
  dot = b1.dotprod(b2);
  va = std::acos(dot);
  double direc = (b1.cross(b2)).dotprod(e2);
  if(direc<0){
    kapTau.second = -1.0*2.0*std::sin(va/2.0)/l1;
  }else{
    kapTau.second = 2.0*std::sin(va/2.0)/l1;
  }
  return kapTau;
}

point randomMol::binorm(point &T1,point &T2){
  point cp = T1.cross(T2);
  double cpsize = cp.length();
  if(cpsize>0.0000000001){
    cp.normalise();
  }else{
    point vert(0.0,0.0,1.0);
    if(std::abs(T1.dotprod(vert))<0.99999999){
      point vert(0.0,0.0,1.0);
       point cp = vert.cross(T1);
       cp.normalise();
    }else{
      point vert(0.0,1.0,0.0);
      cp = vert.cross(T1);
      cp.normalise();
    }
  }
  return cp;
}




point randomMol::getNextPoint(point &pt1,point &pt2,point &pt3,double &r,double &kap,double &tau){
  std::vector<double> theSols;
  std::vector<double> theSols2;
  point e1 = pt2-pt1;
  double le1 = e1.length();
  point e2 = pt3-pt2;double le2 = e2.length();
  e1.normalise();e2.normalise();
  point vert(0.0,0.0,1.0);
  point N1,N2;
  //e1.printPoint();
  //e2.printPoint();
  if(std::abs(e1.dotprod(vert))<0.99999999){
    N1 = vert.cross(e2);
    N1.normalise();
    N2 = e2.cross(N1);
    N2.normalise();
  }else{
    point vert2(0.0,1.0,0.0);
    N1 = vert2.cross(e2);
    N1.normalise();
    N2 = e2.cross(N1);
    N2.normalise();
  }
  double va = 0.25*(le2+r)*tau;
  double beta = 2.0*va*std::sqrt(1-va*va);
  point b1 = binorm(e1,e2);
  point cp1 = b1.cross(N2);point cp2 =b1.cross(N1);
  //cp1.printPoint();cp2.printPoint();
  double a = cp1.dotprod(cp1);double b = cp2.dotprod(cp2);double c = (-2.0)*cp1.dotprod(cp2);
  double d = 0.5*(a-b);double e = 0.5*c;double f = 0.5*(a+b);double R = std::sqrt(d*d + e*e);
  double RHS = (beta*beta-f)/R;
  double phi1,phi2;
  //std::cout<<"why fail "<<std::abs(RHS)<<" "<<R<<"\n";
  if(std::abs(RHS)<1.0){
    if(std::abs(e)>0.000001){ 
      double psi = std::atan2(d,e);
      double as = asin(RHS);
      phi1 = 0.5*(as-psi);phi2 = 0.5*(3.14159265359-as-psi);
    }else{
      phi1 = 0.5*std::acos(RHS);phi2 = 2.0*3.14159265359-phi1;
    }
  }else{
    //cannot solve
    phi1=10000.0;phi2=10000.0;
  }
  if(phi1== 10000.0 && phi2 ==10000.0){
    //cannot solve for these values
    point pout(10000.0,10000.0,10000.0);
    return pout;
  }else{
    // now solve for theta
    double r2 = pt2.eDist(pt3);
    point pc = pt3 -(pt1+pt2)*0.5;
    point cpc1 = N1.cross(pc);point cpc2=N2.cross(pc);point cpcTan = e2.cross(pc);
    double TCPsq = cpcTan.dotprod(cpcTan);double N1CPsq = cpc1.dotprod(cpc1);
    double N2CPsq = cpc2.dotprod(cpc2);double N1N2CP = N1.dotprod(cpc2);
    double N1CPN2CP = cpc1.dotprod(cpc2);double N2N1CP=N2.dotprod(cpc1);double N2TCP = N2.dotprod(cpcTan);
    double N2CPTCP = cpc2.dotprod(cpcTan);double N1CPTCP =cpc1.dotprod(cpcTan);double N1TCP = N1.dotprod(cpcTan);
    double rsq = r*r;
    double cosphi = std::cos(phi1);double sinphi = std::sin(phi1);
    double em2sqno = (rsq + r2*r2)*0.25;double em2sqcos = r2*r*0.5;
    double ejsqno = 0.25*rsq + pc.dotprod(pc);double ejsqcos = e2.dotprod(pc)*r;double ejsqsin= r*(N1.dotprod(pc)*cosphi+ N2.dotprod(pc)*sinphi);
    double areaTermno = r2*r2*TCPsq*0.25;
    double areaTermsinsq = 0.25*rsq*(N1CPsq + N2N1CP*r2 + 0.25*r2*r2)*cosphi*cosphi + 0.5*rsq*N1CPN2CP*cosphi*sinphi+0.25*rsq*(N2CPsq - N1N2CP*r2 + 0.25*r2*r2)*sinphi*sinphi;
    double areaTermcossq = rsq*0.25*TCPsq;double areaTermcos = 0.5*r*r2*TCPsq;
    double areaTermsin = 0.25*r*r2*(2.0*N1CPTCP + N2TCP*r2)*cosphi + 0.25*r*r2*(2.0*N2CPTCP-r2*N1TCP)*sinphi;
    double areaTermsincos = 0.25*rsq*(2.0*N1CPTCP + r2*N2TCP)*cosphi + 0.25*rsq*(2.0*N2CPTCP-r2*N1TCP)*sinphi;
    point em1sqP = (pt3+pt2)*0.5 -(pt1+pt2)*0.5;double em1Sq = em1sqP.dotprod(em1sqP);
    double kapTermnc = kap*kap*em1Sq*em2sqno*ejsqno;double kapTermcos = kap*kap*em1Sq*(em2sqcos*ejsqno+em2sqno*ejsqcos);double kapTermsin = kap*kap*em1Sq*em2sqno*ejsqsin;double kapTermsincos = kap*kap*em1Sq*em2sqcos*ejsqsin;
    double kapTermcossq = kap*kap*em1Sq*em2sqcos*ejsqcos;
    double noTerm = 4.0*areaTermno-kapTermnc;double sinTerm = 4.0*areaTermsin-kapTermsin;double cosTerm = 4.0*areaTermcos-kapTermcos;double sincosTerm = 4.0*areaTermsincos-kapTermsincos;double sinsqTerm =4.0*areaTermsinsq;double cossqTerm = 4.0*areaTermcossq-kapTermcossq;
    // now solve for phi first check if we can
    double val = noTerm +  cosTerm + cossqTerm;
    double prevval=val;
    std::vector<int> bi;
    int n=59;
    for(int i=1;i<=n;i++){
      double the =2.0*3.14159265359*i/(n);
      double st = std::sin(the);double ct =std::cos(the);
      val = noTerm + sinTerm*st + cosTerm*ct + sincosTerm*st*ct+sinsqTerm*st*st+cossqTerm*ct*ct;
      if(val*prevval<=0.0){
	bi.push_back(i);
      }
      prevval=val;
    }
    if(bi.size()>0){
      for(int i=0;i<bi.size();i++){
	double sol = 3.14159265359*(bi[i]-1)/(n)+3.14159265359*(bi[i])/(n);
	for(int j=0;j<=4;j++){
	  double st = std::sin(sol);double ct =std::cos(sol);
	  double f = noTerm + sinTerm*st + cosTerm*ct + sincosTerm*st*ct+sinsqTerm*st*st+cossqTerm*ct*ct;
	  double fdash = sinTerm*ct - cosTerm*st + sincosTerm*(ct*ct-st*st)+ sinsqTerm*2.0*st*ct - cossqTerm*2.0*st*ct; 
	  sol =  sol-f/fdash;   
	}
	  theSols.push_back(sol);
      }
    }
    //now for second phi
    cosphi = std::cos(phi2);sinphi = std::sin(phi2);
    ejsqsin= r*(N1.dotprod(pc)*cosphi+ N2.dotprod(pc)*sinphi);
    areaTermsinsq = 0.25*rsq*(N1CPsq + N2N1CP*r2 + 0.25*r2*r2)*cosphi*cosphi + 0.5*rsq*N1CPN2CP*cosphi*sinphi+0.25*rsq*(N2CPsq - N1N2CP*r2 + 0.25*r2*r2)*sinphi*sinphi;
    areaTermsin = 0.25*r*r2*(2.0*N1CPTCP + N2TCP*r2)*cosphi + 0.25*r*r2*(2.0*N2CPTCP-r2*N1TCP)*sinphi;
    areaTermsincos = 0.25*rsq*(2.0*N1CPTCP + r2*N2TCP)*cosphi + 0.25*rsq*(2.0*N2CPTCP-r2*N1TCP)*sinphi;
    kapTermnc = kap*kap*em1Sq*em2sqno*ejsqno;kapTermcos = kap*kap*em1Sq*(em2sqcos*ejsqno+em2sqno*ejsqcos);kapTermsin = kap*kap*em1Sq*em2sqno*ejsqsin;kapTermsincos = kap*kap*em1Sq*em2sqcos*ejsqsin;
    kapTermcossq = kap*kap*em1Sq*em2sqcos*ejsqcos;
    noTerm = 4.0*areaTermno-kapTermnc;sinTerm = 4.0*areaTermsin-kapTermsin;cosTerm = 4.0*areaTermcos-kapTermcos;sincosTerm = 4.0*areaTermsincos-kapTermsincos;sinsqTerm =4.0*areaTermsinsq;cossqTerm = 4.0*areaTermcossq-kapTermcossq;
    // now solve for phi first check if we can
    val = noTerm +  cosTerm + cossqTerm;
    prevval=val;
    std::vector<int> bi2;
    for(int i=1;i<=n;i++){
      double the =2.0*3.14159265359*i/(n);
      double st = std::sin(the);double ct =std::cos(the);
      val = noTerm + sinTerm*st + cosTerm*ct + sincosTerm*st*ct+sinsqTerm*st*st+cossqTerm*ct*ct;
      if(val*prevval<=0.0){
	bi2.push_back(i);
      }
      prevval=val;
    }
    if(bi2.size()>0){
      //newtons method to find solutions
      for(int i=0;i<bi.size();i++){
	double sol = 3.14159265359*(bi2[i]-1)/(n)+3.14159265359*(bi2[i])/(n);
	for(int j=0;j<=4;j++){
	  double st = std::sin(sol);double ct =std::cos(sol);
	  double f = noTerm + sinTerm*st + cosTerm*ct + sincosTerm*st*ct+sinsqTerm*st*st+cossqTerm*ct*ct;
	  double fdash = sinTerm*ct - cosTerm*st + sincosTerm*(ct*ct-st*st)+ sinsqTerm*2.0*st*ct - cossqTerm*2.0*st*ct; 
	  sol =  sol-f/fdash;   
	}
	theSols2.push_back(sol);
      }
    }    
    std::vector<point> pt4List;
    for(int i=0;i<theSols.size();i++){
      double ct = std::cos(theSols[i]);double st= std::sin(theSols[i]);
      double cp = std::cos(phi1);double sp= std::sin(phi1);
      point tan = e2*ct + N1*st*cp + N2*st*sp;
      point pt4 = pt3 + tan*r;
      std::pair<double,double> ktpair = kapTau(pt1,pt2,pt3,pt4);
      if(std::abs(ktpair.first-kap)<0.01&&std::abs(ktpair.second-tau)<0.01){
        pt4List.push_back(pt4);
      }
    }
    for(int i=0;i<theSols2.size();i++){
      double ct = std::cos(theSols2[i]);double st= std::sin(theSols2[i]);
      double cp = std::cos(phi2);double sp= std::sin(phi2);
      point tan = e2*ct + N1*st*cp + N2*st*sp;
      point pt4 = pt3 + tan*r;
      std::pair<double,double> ktpair = kapTau(pt1,pt2,pt3,pt4);
      if(std::abs(ktpair.first-kap)<0.01&&std::abs(ktpair.second-tau)<0.01){
        pt4List.push_back(pt4);
      }
    }
    if(pt4List.size()>0){
      int outIndex;
      double maxDist=0.0;
      for(int i=0;i<pt4List.size();i++){
	      point pd = pt4List[i]-pt2;
          if(pd.length()>=maxDist){
	      maxDist=pd.length();
	      outIndex=i;
	    }	
      }
      return pt4List[outIndex];
    }else{
      //no solution
      point pout(10000.0,10000.0,10000.0);
      return pout;
    }
  }
}



point randomMol::getNextPointTest(point &pt1,point &pt2,point &pt3,double &r,double &kap,double &tau){
  std::vector<double> theSols;
  std::vector<double> theSols2;
  pt1.printPoint();
  pt2.printPoint();
  pt3.printPoint();
  point e1 = pt2-pt1;
  double le1 = e1.length();
  point e2 = pt3-pt2;double le2 = e2.length();
  e1.normalise();e2.normalise();
  e1.printPoint();
  e2.printPoint();
  point vert(0.0,0.0,1.0);
  point N1,N2;
  if(std::abs(e1.dotprod(vert))<0.99999999){
    N1 = vert.cross(e2);
    N1.normalise();
    N2 = e2.cross(N1);
    N2.normalise();
  }else{
    point vert2(0.0,1.0,0.0);
    N1 = vert2.cross(e2);
    N1.normalise();
    N2 = e2.cross(N1);
    N2.normalise();
  }
  double va = 0.25*(le2+r)*tau;
  std::cout<<"va is "<<va<<"\n";
  double beta = 2.0*va*std::sqrt(1-va*va);
  point b1 = binorm(e1,e2);
  point cp1 = b1.cross(N2);point cp2 =b1.cross(N1);
  double a = cp1.dotprod(cp1);double b = cp2.dotprod(cp2);double c = (-2.0)*cp1.dotprod(cp2);
  double d = 0.5*(a-b);double e = 0.5*c;double f = 0.5*(a+b);double R = std::sqrt(d*d + e*e);
  double RHS = (beta*beta-f)/R;
  std::cout<<RHS<<" "<<R<<" rhs ? \n";
  double phi1,phi2;
  if(std::abs(RHS)<1.0){
    if(std::abs(e)>0.000001){ 
      double psi = std::atan2(d,e);
      double as = asin(RHS);
      phi1 = 0.5*(as-psi);phi2 = 0.5*(3.14159265359-as-psi);
    }else{
      phi1 = 0.5*std::acos(RHS);phi2 = 2.0*3.14159265359-phi1;
    }
  }else{
    //cannot solve
    phi1=10000.0;phi2=10000.0;
  }
  if(phi1== 10000.0 && phi2 ==10000.0){
    //cannot solve for these values
    point pout(10000.0,10000.0,10000.0);
    return pout;
  }else{
    // now solve for theta
    double r2 = pt2.eDist(pt3);
    point pc = pt3 -(pt1+pt2)*0.5;
    point cpc1 = N1.cross(pc);point cpc2=N2.cross(pc);point cpcTan = e2.cross(pc);
    double TCPsq = cpcTan.dotprod(cpcTan);double N1CPsq = cpc1.dotprod(cpc1);
    double N2CPsq = cpc2.dotprod(cpc2);double N1N2CP = N1.dotprod(cpc2);
    double N1CPN2CP = cpc1.dotprod(cpc2);double N2N1CP=N2.dotprod(cpc1);double N2TCP = N2.dotprod(cpcTan);
    double N2CPTCP = cpc2.dotprod(cpcTan);double N1CPTCP =cpc1.dotprod(cpcTan);double N1TCP = N1.dotprod(cpcTan);
    double rsq = r*r;
    double cosphi = std::cos(phi1);double sinphi = std::sin(phi1);
    double em2sqno = (rsq + r2*r2)*0.25;double em2sqcos = r2*r*0.5;
    double ejsqno = 0.25*rsq + pc.dotprod(pc);double ejsqcos = e2.dotprod(pc)*r;double ejsqsin= r*(N1.dotprod(pc)*cosphi+ N2.dotprod(pc)*sinphi);
    double areaTermno = r2*r2*TCPsq*0.25;
    double areaTermsinsq = 0.25*rsq*(N1CPsq + N2N1CP*r2 + 0.25*r2*r2)*cosphi*cosphi + 0.5*rsq*N1CPN2CP*cosphi*sinphi+0.25*rsq*(N2CPsq - N1N2CP*r2 + 0.25*r2*r2)*sinphi*sinphi;
    double areaTermcossq = rsq*0.25*TCPsq;double areaTermcos = 0.5*r*r2*TCPsq;
    double areaTermsin = 0.25*r*r2*(2.0*N1CPTCP + N2TCP*r2)*cosphi + 0.25*r*r2*(2.0*N2CPTCP-r2*N1TCP)*sinphi;
    double areaTermsincos = 0.25*rsq*(2.0*N1CPTCP + r2*N2TCP)*cosphi + 0.25*rsq*(2.0*N2CPTCP-r2*N1TCP)*sinphi;
    point em1sqP = (pt3+pt2)*0.5 -(pt1+pt2)*0.5;double em1Sq = em1sqP.dotprod(em1sqP);
    double kapTermnc = kap*kap*em1Sq*em2sqno*ejsqno;double kapTermcos = kap*kap*em1Sq*(em2sqcos*ejsqno+em2sqno*ejsqcos);double kapTermsin = kap*kap*em1Sq*em2sqno*ejsqsin;double kapTermsincos = kap*kap*em1Sq*em2sqcos*ejsqsin;
    double kapTermcossq = kap*kap*em1Sq*em2sqcos*ejsqcos;
    double noTerm = 4.0*areaTermno-kapTermnc;double sinTerm = 4.0*areaTermsin-kapTermsin;double cosTerm = 4.0*areaTermcos-kapTermcos;double sincosTerm = 4.0*areaTermsincos-kapTermsincos;double sinsqTerm =4.0*areaTermsinsq;double cossqTerm = 4.0*areaTermcossq-kapTermcossq;
    // now solve for phi first check if we can
    double val = noTerm +  cosTerm + cossqTerm;
    double prevval=val;
    std::vector<int> bi;
    int n=59;
    for(int i=1;i<=n;i++){
      double the =2.0*3.14159265359*i/(n);
      double st = std::sin(the);double ct =std::cos(the);
      val = noTerm + sinTerm*st + cosTerm*ct + sincosTerm*st*ct+sinsqTerm*st*st+cossqTerm*ct*ct;
      if(val*prevval<=0.0){
	bi.push_back(i);
      }
      prevval=val;
    }
    if(bi.size()>0){
      for(int i=0;i<bi.size();i++){
	double sol = 3.14159265359*(bi[i]-1)/(n)+3.14159265359*(bi[i])/(n);
	for(int j=0;j<=4;j++){
	  double st = std::sin(sol);double ct =std::cos(sol);
	  double f = noTerm + sinTerm*st + cosTerm*ct + sincosTerm*st*ct+sinsqTerm*st*st+cossqTerm*ct*ct;
	  double fdash = sinTerm*ct - cosTerm*st + sincosTerm*(ct*ct-st*st)+ sinsqTerm*2.0*st*ct - cossqTerm*2.0*st*ct; 
	  sol =  sol-f/fdash;   
	}
	  theSols.push_back(sol);
      }
    }
    //now for second phi
    cosphi = std::cos(phi2);sinphi = std::sin(phi2);
    ejsqsin= r*(N1.dotprod(pc)*cosphi+ N2.dotprod(pc)*sinphi);
    areaTermsinsq = 0.25*rsq*(N1CPsq + N2N1CP*r2 + 0.25*r2*r2)*cosphi*cosphi + 0.5*rsq*N1CPN2CP*cosphi*sinphi+0.25*rsq*(N2CPsq - N1N2CP*r2 + 0.25*r2*r2)*sinphi*sinphi;
    areaTermsin = 0.25*r*r2*(2.0*N1CPTCP + N2TCP*r2)*cosphi + 0.25*r*r2*(2.0*N2CPTCP-r2*N1TCP)*sinphi;
    areaTermsincos = 0.25*rsq*(2.0*N1CPTCP + r2*N2TCP)*cosphi + 0.25*rsq*(2.0*N2CPTCP-r2*N1TCP)*sinphi;
    kapTermnc = kap*kap*em1Sq*em2sqno*ejsqno;kapTermcos = kap*kap*em1Sq*(em2sqcos*ejsqno+em2sqno*ejsqcos);kapTermsin = kap*kap*em1Sq*em2sqno*ejsqsin;kapTermsincos = kap*kap*em1Sq*em2sqcos*ejsqsin;
    kapTermcossq = kap*kap*em1Sq*em2sqcos*ejsqcos;
    noTerm = 4.0*areaTermno-kapTermnc;sinTerm = 4.0*areaTermsin-kapTermsin;cosTerm = 4.0*areaTermcos-kapTermcos;sincosTerm = 4.0*areaTermsincos-kapTermsincos;sinsqTerm =4.0*areaTermsinsq;cossqTerm = 4.0*areaTermcossq-kapTermcossq;
    // now solve for phi first check if we can
    val = noTerm +  cosTerm + cossqTerm;
    prevval=val;
    std::vector<int> bi2;
    for(int i=1;i<=n;i++){
      double the =2.0*3.14159265359*i/(n);
      double st = std::sin(the);double ct =std::cos(the);
      val = noTerm + sinTerm*st + cosTerm*ct + sincosTerm*st*ct+sinsqTerm*st*st+cossqTerm*ct*ct;
      if(val*prevval<=0.0){
	bi2.push_back(i);
      }
      prevval=val;
    }
    if(bi2.size()>0){
      //newtons method to find solutions
      for(int i=0;i<bi.size();i++){
	double sol = 3.14159265359*(bi2[i]-1)/(n)+3.14159265359*(bi2[i])/(n);
	for(int j=0;j<=4;j++){
	  double st = std::sin(sol);double ct =std::cos(sol);
	  double f = noTerm + sinTerm*st + cosTerm*ct + sincosTerm*st*ct+sinsqTerm*st*st+cossqTerm*ct*ct;
	  double fdash = sinTerm*ct - cosTerm*st + sincosTerm*(ct*ct-st*st)+ sinsqTerm*2.0*st*ct - cossqTerm*2.0*st*ct; 
	  sol =  sol-f/fdash;   
	}
	theSols2.push_back(sol);
      }
    }    
    std::vector<point> pt4List;
    for(int i=0;i<theSols.size();i++){
      double ct = std::cos(theSols[i]);double st= std::sin(theSols[i]);
      double cp = std::cos(phi1);double sp= std::sin(phi1);
      point tan = e2*ct + N1*st*cp + N2*st*sp;
      point pt4 = pt3 + tan*r;
      std::pair<double,double> ktpair = kapTau(pt1,pt2,pt3,pt4);
      if(std::abs(ktpair.first-kap)<0.01&&std::abs(ktpair.second-tau)<0.01){
        pt4List.push_back(pt4);
      }
    }
    for(int i=0;i<theSols2.size();i++){
      double ct = std::cos(theSols2[i]);double st= std::sin(theSols2[i]);
      double cp = std::cos(phi2);double sp= std::sin(phi2);
      point tan = e2*ct + N1*st*cp + N2*st*sp;
      point pt4 = pt3 + tan*r;
      std::pair<double,double> ktpair = kapTau(pt1,pt2,pt3,pt4);
      if(std::abs(ktpair.first-kap)<0.01&&std::abs(ktpair.second-tau)<0.01){
        pt4List.push_back(pt4);
      }
    }
    if(pt4List.size()>0){
      int outIndex;
      double maxDist=0.0;
      for(int i=0;i<pt4List.size();i++){
	      point pd = pt4List[i]-pt2;
          if(pd.length()>=maxDist){
	      maxDist=pd.length();
	      outIndex=i;
	    }	
      }
      return pt4List[outIndex];
    }else{
      //no solution
      point pout(10000.0,10000.0,10000.0);
      return pout;
    }
  }
}



double randomMol::interpVel2D(double &x0,double &y0,double &x,double &y,std::vector<std::vector<double> > &vx,double &d1,double &d2,int &nx,int &ny){
  //find the x and y indicies to the left of the actual value
  double iv = ((x -x0)/d1)*double(nx-1);
  double jv = ((y -y0)/d2)*double(ny-1);
  int i = int(floor(iv));
  int j = int(floor(jv));
  if(std::abs(i+1-iv)<0.00001){
    i=i+1;
  }
  if(std::abs(j+1-jv)<0.00001){
    j=j+1;
  }  
  if(i>= nx-1){
    i = nx-2;
  }
  if(j>= ny-1){
    j = ny-2;
  }
  if(i< 0){
    i = 0;
  }

  if(j<0){
    j = 0;
  }
  int ip = i+1;
  int jp = j+1;
  // get the distance from the left x/y grid coordinate
  double xd = (x-((d1*i)/(nx-1) + x0))*(nx-1)/d1;
  double yd = (y-((d2*j)/(ny-1) + y0))*(ny-1)/d2;;
  //grab the required velocity grid terms
  double outVal = vx[i][j] + (vx[ip][j]-vx[i][j])*xd + (vx[i][jp]-vx[i][j])*yd + (vx[ip][jp]+vx[i][j]-vx[ip][j]-vx[i][jp])*xd*yd;
  return outVal;
}

double randomMol::interpVel1D(double &x0,double &x,std::vector<double> &vx,double &d1,int &nx){
  //find the x and y indicies to the left of the actual value
  double iv = ((x -x0)/d1)*double(nx-1);
  int i = int(floor(iv));
  if(std::abs(i+1-iv)<0.00001){
    i=i+1;
  }
  if(i>= nx-1){
    i = nx-2;
  }
  if(i< 0){
    i = 0;
  }
  int ip = i+1;
  // get the distance from the left x/y grid coordinate
  double xd = (x-((d1*i)/(nx-1) + x0))*(nx-1)/d1;
  //grab the required velocity grid terms
  double outVal = vx[i] + (vx[ip]-vx[i])*xd;
  return outVal;
}

double randomMol::cumDistMixed(double x){
  return interpVel1D(kapL,x,mixedXCumInterpolants,dkap,nkap);
}

double randomMol::cumDistNegBeta(double x){
  return interpVel1D(kapL,x,negBetaXCumInterpolants,dkap,nkap);
}

double randomMol::cumDistPosBeta(double x){
  return interpVel1D(kapL,x,posBetaXCumInterpolants,dkap,nkap);
}

double randomMol::cumDistAlpha(double x){
  return interpVel1D(kapL,x,alphaXCumInterpolants,dkap,nkap);
}

double randomMol::cumDistMixedStrand(double x){
  return interpVel1D(kapL,x,mixedStrandXCumInterpolants,dkap,nkap);
}

double randomMol::cumDistNegBetaStrand(double x){
  return interpVel1D(kapL,x,negBetaStrandXCumInterpolants,dkap,nkap);
}

double randomMol::cumDistPosBetaStrand(double x){
  return interpVel1D(kapL,x,posBetaStrandXCumInterpolants,dkap,nkap);
}

double randomMol::cumDistAlphaStrand(double x){
  return interpVel1D(kapL,x,alphaStrandXCumInterpolants,dkap,nkap);
}

double randomMol::cumDistLinkerToHelix(double x){
  return interpVel1D(kapL,x,linkerToHelixXCumInterpolants,dkap,nkap);
}

double randomMol::cumDistHelixToLinker(double x){
  return interpVel1D(kapL,x,helixToLinkerXCumInterpolants,dkap,nkap);
}

double randomMol::cumDistStrandToHelix(double x){
  return interpVel1D(kapL,x,strandToHelixXCumInterpolants,dkap,nkap);
}

double randomMol::cumDistHelixToStrand(double x){
  return interpVel1D(kapL,x,helixToStrandXCumInterpolants,dkap,nkap);
}

double randomMol::cumDistStrandToLinker(double x){
  return interpVel1D(kapL,x,strandToLinkerXCumInterpolants,dkap,nkap);
}

double randomMol::cumDistLinkerToStrand(double x){
  return interpVel1D(kapL,x,linkerToStrandXCumInterpolants,dkap,nkap);
}



/*Note need to make sure kapVal is set */

double randomMol::cumXGYMixed(double Y,double kapValIn){
  return interpVel2D(kapL,tauL,kapValIn,Y,mixedYfromXInterpolants,dkap,dtau,nkap,ntau);
}

double randomMol::cumXGYNegBeta(double Y,double kapValIn){
  return interpVel2D(kapL,tauL,kapValIn,Y,negBetaYfromXInterpolants,dkap,dtau,nkap,ntau);
}

double randomMol::cumXGYPosBeta(double Y,double kapValIn){
  return interpVel2D(kapL,tauL,kapValIn,Y,posBetaYfromXInterpolants,dkap,dtau,nkap,ntau);
}

double randomMol::cumXGYAlpha(double Y,double kapValIn){
  return interpVel2D(kapL,tauL,kapValIn,Y,alphaYfromXInterpolants,dkap,dtau,nkap,ntau);
}

double randomMol::cumXGYMixedStrand(double Y,double kapValIn){
  return interpVel2D(kapL,tauL,kapValIn,Y,mixedStrandYfromXInterpolants,dkap,dtau,nkap,ntau);
}

double randomMol::cumXGYNegBetaStrand(double Y,double kapValIn){
  return interpVel2D(kapL,tauL,kapValIn,Y,negBetaStrandYfromXInterpolants,dkap,dtau,nkap,ntau);
}

double randomMol::cumXGYPosBetaStrand(double Y,double kapValIn){
  return interpVel2D(kapL,tauL,kapValIn,Y,posBetaStrandYfromXInterpolants,dkap,dtau,nkap,ntau);
}

double randomMol::cumXGYAlphaStrand(double Y,double kapValIn){
  return interpVel2D(kapL,tauL,kapVal,Y,alphaStrandYfromXInterpolants,dkap,dtau,nkap,ntau);
}

double randomMol::cumXGYLinkerToHelix(double Y,double kapValIn){
  return interpVel2D(kapL,tauL,kapValIn,Y,linkerToHelixYfromXInterpolants,dkap,dtau,nkap,ntau);
}

double randomMol::cumXGYHelixToLinker(double Y,double kapValIn){
  return interpVel2D(kapL,tauL,kapValIn,Y,helixToLinkerYfromXInterpolants,dkap,dtau,nkap,ntau);
}

double randomMol::cumXGYStrandToHelix(double Y,double kapValIn){
  return interpVel2D(kapL,tauL,kapValIn,Y,strandToHelixYfromXInterpolants,dkap,dtau,nkap,ntau);
}

double randomMol::cumXGYHelixToStrand(double Y,double kapValIn){
  return interpVel2D(kapL,tauL,kapVal,Y,helixToStrandYfromXInterpolants,dkap,dtau,nkap,ntau);
}

double randomMol::cumXGYStrandToLinker(double Y,double kapValIn){
  return interpVel2D(kapL,tauL,kapValIn,Y,strandToLinkerYfromXInterpolants,dkap,dtau,nkap,ntau);
}

double randomMol::cumXGYLinkerToStrand(double Y,double kapValIn){
  return interpVel2D(kapL,tauL,kapValIn,Y,linkerToStrandYfromXInterpolants,dkap,dtau,nkap,ntau);
}


double randomMol::bisection(double (randomMol::*func)(double),double &val,int n,double rangeLower,double rangeUpper,double &prevVal){
  double lowerVal = rangeLower;double upperVal =rangeUpper;
  double midVal = 0.5*(lowerVal +upperVal);
  int j=1;
  while(j<=n){
    if((*this.*func)(midVal)>val){
      upperVal=midVal;
    }else{
      lowerVal=midVal;
    }
    midVal = 0.5*(upperVal +lowerVal);
    j++;
  }
  
  double a = (*this.*func)(lowerVal);double b = (*this.*func)(upperVal);
  double t = (val-a)/(b-a);
  if(std::isnan(t)!=true && std::isinf(t)!=true){
    double t = (val-a)/(b-a);
    return lowerVal + t*(upperVal-lowerVal);
  }else{
    return prevVal;
  }
}

double randomMol::bisectionDoubleArg(double (randomMol::*func)(double,double),double &val,int n,double rangeLower,double rangeUpper,double secondArg,double &prevVal){
  double lowerVal = rangeLower;double upperVal =rangeUpper;
  double midVal = 0.5*(lowerVal +upperVal);
  int j=1;
  while(j<=n){
    if((*this.*func)(midVal,secondArg)>val){
      upperVal=midVal;
    }else{
      lowerVal=midVal;
    }
    midVal = 0.5*(upperVal +lowerVal);
    j++;
  }
  double a = (*this.*func)(lowerVal,secondArg);double b = (*this.*func)(upperVal,secondArg);
   double t = (val-a)/(b-a);
  if(std::isnan(t)!=true&& std::isinf(t)!=true){
    double t = (val-a)/(b-a);
    return lowerVal + t*(upperVal-lowerVal);
  }else{
      return prevVal;
  }
}



std::pair<double,double> randomMol::getKapTauProb(int &type,std::default_random_engine &generator){
  std::pair<double,double> ktpair;//need to define
  double ranUnitary =distribution(generator);
  double (randomMol::*func)(double);
  double (randomMol::*func2)(double,double);
  // in case of the bisection failing, set fixed values
  double tauVal;
  if(type<=6){
    // mixed set to alpha value
    kapVal = alphaKap;
    tauVal = alphaTau;
  }else{
    kapVal = betaKap;
    tauVal = betaTau;
  }
  int nGridSize=100;
  double xMin = 0.0;double xMax = 0.7;
  double yMin = -0.7;double yMax=0.7;
  if(type==1){
    func = &randomMol::cumDistMixed;
    kapVal = bisection(func,ranUnitary,nGridSize,xMin,xMax,kapVal);
    func2 = &randomMol::cumXGYMixed;
    ranUnitary = distribution(generator);
    tauVal = bisectionDoubleArg(func2,ranUnitary,nGridSize,yMin,yMax,kapVal,tauVal);
  }else if(type==2){
    func = &randomMol::cumDistNegBeta;
    kapVal = bisection(func,ranUnitary,nGridSize,xMin,xMax,kapVal);
    func2 = &randomMol::cumXGYNegBeta;
    ranUnitary = distribution(generator);
    tauVal = bisectionDoubleArg(func2,ranUnitary,nGridSize,yMin,yMax,kapVal,tauVal);
  }else if(type==3){
    func = &randomMol::cumDistPosBeta;
    kapVal = bisection(func,ranUnitary,nGridSize,xMin,xMax,kapVal);
    func2 = &randomMol::cumXGYPosBeta;
    ranUnitary = distribution(generator);
    tauVal = bisectionDoubleArg(func2,ranUnitary,nGridSize,yMin,yMax,kapVal,tauVal);
  }else if(type==4){
    func = &randomMol::cumDistAlpha;
    kapVal = bisection(func,ranUnitary,nGridSize,xMin,xMax,kapVal);
    func2 = &randomMol::cumXGYAlpha;
    ranUnitary = distribution(generator);
    tauVal = bisectionDoubleArg(func2,ranUnitary,nGridSize,yMin,yMax,kapVal,tauVal);
  }else if(type==5){
    func = &randomMol::cumDistLinkerToHelix;
    kapVal = bisection(func,ranUnitary,nGridSize,xMin,xMax,kapVal);
    func2 = &randomMol::cumXGYLinkerToHelix;
    ranUnitary = distribution(generator);
    tauVal = bisectionDoubleArg(func2,ranUnitary,nGridSize,yMin,yMax,kapVal,tauVal);
  }else if(type==6){
    func = &randomMol::cumDistHelixToLinker;
    kapVal = bisection(func,ranUnitary,nGridSize,xMin,xMax,kapVal);
    func2 = &randomMol::cumXGYHelixToLinker;
    ranUnitary = distribution(generator);
    tauVal = bisectionDoubleArg(func2,ranUnitary,nGridSize,yMin,yMax,kapVal,tauVal);
  }else if(type==7){
    func = &randomMol::cumDistMixedStrand;
    kapVal = bisection(func,ranUnitary,nGridSize,xMin,xMax,kapVal);
    func2 = &randomMol::cumXGYMixedStrand;
    ranUnitary = distribution(generator);
    tauVal = bisectionDoubleArg(func2,ranUnitary,nGridSize,yMin,yMax,kapVal,tauVal);
  }else if(type==8){
    func = &randomMol::cumDistNegBetaStrand;
    kapVal = bisection(func,ranUnitary,nGridSize,xMin,xMax,kapVal);
    func2 = &randomMol::cumXGYNegBetaStrand;
    ranUnitary = distribution(generator);
    tauVal = bisectionDoubleArg(func2,ranUnitary,nGridSize,yMin,yMax,kapVal,tauVal);
  }else if(type==9){
    func = &randomMol::cumDistPosBetaStrand;
    kapVal = bisection(func,ranUnitary,nGridSize,xMin,xMax,kapVal);
    func2 = &randomMol::cumXGYPosBetaStrand;
    ranUnitary = distribution(generator);
    tauVal = bisectionDoubleArg(func2,ranUnitary,nGridSize,yMin,yMax,kapVal,tauVal);
  }else if(type==10){
    func = &randomMol::cumDistAlphaStrand;
    kapVal = bisection(func,ranUnitary,nGridSize,xMin,xMax,kapVal);
    func2 = &randomMol::cumXGYAlphaStrand;
    ranUnitary = distribution(generator);
    tauVal = bisectionDoubleArg(func2,ranUnitary,nGridSize,yMin,yMax,kapVal,tauVal);
  }else if(type==11){
    func = &randomMol::cumDistStrandToHelix;
    kapVal = bisection(func,ranUnitary,nGridSize,xMin,xMax,kapVal);
    func2 = &randomMol::cumXGYStrandToHelix;
    ranUnitary = distribution(generator);
    tauVal = bisectionDoubleArg(func2,ranUnitary,nGridSize,yMin,yMax,kapVal,tauVal);
  }else if(type==12){
    func = &randomMol::cumDistHelixToStrand;
    kapVal = bisection(func,ranUnitary,nGridSize,xMin,xMax,kapVal);
    func2 = &randomMol::cumXGYHelixToStrand;
    ranUnitary = distribution(generator);
    tauVal = bisectionDoubleArg(func2,ranUnitary,nGridSize,yMin,yMax,kapVal,tauVal);
  }else if(type==13){
    func = &randomMol::cumDistLinkerToStrand;
    kapVal = bisection(func,ranUnitary,nGridSize,xMin,xMax,kapVal);
    func2 = &randomMol::cumXGYLinkerToStrand;
    ranUnitary = distribution(generator);
    tauVal = bisectionDoubleArg(func2,ranUnitary,nGridSize,yMin,yMax,kapVal,tauVal);
  }else{
    func = &randomMol::cumDistStrandToLinker;
    kapVal = bisection(func,ranUnitary,nGridSize,xMin,xMax,kapVal);
    func2 = &randomMol::cumXGYStrandToLinker;
    ranUnitary = distribution(generator);
    tauVal = bisectionDoubleArg(func2,ranUnitary,nGridSize,yMin,yMax,kapVal,tauVal);
  }
  
  ktpair.first = kapVal;ktpair.second=tauVal;
  return ktpair;
}

std::pair<double,double> randomMol::getProbFromKapTau(int &type,std::default_random_engine &generator,std::pair<double,double> &kapTauPair,double probChange){
  std::pair<double,double> ktpair;//need to define
  std::uniform_real_distribution<double> pdist(-probChange,probChange);
  double delProb = pdist(generator);
  double kapValPrev = kapTauPair.first;
  double tauValPrev = kapTauPair.second;
  double  tauVal;
  int nGridSize=100;
  double prob;
  double (randomMol::*func)(double);
  double (randomMol::*func2)(double,double);
  double xMin = 0.0;double xMax = 0.7;
  double yMin = -0.7;double yMax=0.7;
  if(type==1){
    prob = cumDistMixed(kapValPrev);
    prob = prob + delProb*prob;
    func = &randomMol::cumDistMixed;
    kapVal = bisection(func,prob,nGridSize,xMin,xMax,kapValPrev);
    prob = cumXGYMixed(tauValPrev,kapValPrev);
    prob = prob + delProb*prob;
    func2 = &randomMol::cumXGYMixed;
    tauVal = bisectionDoubleArg(func2,prob,nGridSize,yMin,yMax,kapVal,tauValPrev);
  }else if(type==2){
    prob = cumDistNegBeta(kapValPrev);
    prob = prob + delProb*prob;
    func = &randomMol::cumDistNegBeta;
    kapVal = bisection(func,prob,nGridSize,xMin,xMax,kapValPrev);
    prob = cumXGYNegBeta(tauValPrev,kapValPrev);
    prob = prob + delProb*prob;
    func2 = &randomMol::cumXGYNegBeta;
    tauVal = bisectionDoubleArg(func2,prob,nGridSize,yMin,yMax,kapVal,tauValPrev);
  }else if(type==3){
    prob = cumDistPosBeta(kapValPrev);
    prob = prob + delProb*prob;
    func = &randomMol::cumDistPosBeta;
    kapVal = bisection(func,prob,nGridSize,xMin,xMax,kapValPrev);
    prob = cumXGYPosBeta(tauValPrev,kapValPrev);
    prob = prob + delProb*prob;
    func2 = &randomMol::cumXGYPosBeta;
    tauVal = bisectionDoubleArg(func2,prob,nGridSize,yMin,yMax,kapVal,tauValPrev);
  }else if(type==4){
    prob = cumDistAlpha(kapValPrev);
    prob = prob + delProb*prob;
    func = &randomMol::cumDistAlpha;
    kapVal = bisection(func,prob,nGridSize,xMin,xMax,kapValPrev);
    prob = cumXGYAlpha(tauValPrev,kapValPrev);
    prob = prob + delProb*prob;
    func2 = &randomMol::cumXGYAlpha;
    tauVal = bisectionDoubleArg(func2,prob,nGridSize,yMin,yMax,kapVal,tauValPrev);
  }else if(type==5){
    prob = cumDistLinkerToHelix(kapValPrev);
    prob = prob + delProb*prob;
    func = &randomMol::cumDistLinkerToHelix;
    kapVal = bisection(func,prob,nGridSize,xMin,xMax,kapValPrev);
    prob = cumXGYLinkerToHelix(tauValPrev,kapValPrev);
    prob = prob + delProb*prob;
    func2 = &randomMol::cumXGYLinkerToHelix;
    tauVal = bisectionDoubleArg(func2,prob,nGridSize,yMin,yMax,kapVal,tauValPrev);
  }else if(type==6){
    prob = cumDistHelixToLinker(kapValPrev);
    prob = prob + delProb*prob;
    func = &randomMol::cumDistHelixToLinker;
    kapVal = bisection(func,prob,nGridSize,xMin,xMax,kapValPrev);
    prob = cumXGYHelixToLinker(tauValPrev,kapValPrev);
    prob = prob + delProb*prob;
    func2 = &randomMol::cumXGYHelixToLinker;
    tauVal = bisectionDoubleArg(func2,prob,nGridSize,yMin,yMax,kapVal,tauValPrev);
  }else if(type==7){
    prob = cumDistMixedStrand(kapValPrev);
    prob = prob + delProb*prob;
    func = &randomMol::cumDistMixedStrand;
    kapVal = bisection(func,prob,nGridSize,xMin,xMax,kapValPrev);
    prob = cumXGYMixedStrand(tauValPrev,kapValPrev);
    prob = prob + delProb*prob;
    func2 = &randomMol::cumXGYMixedStrand;
    tauVal = bisectionDoubleArg(func2,prob,nGridSize,yMin,yMax,kapVal,tauValPrev);
  }else if(type==8){
    prob = cumDistNegBetaStrand(kapValPrev);
    prob = prob + delProb*prob;
    func = &randomMol::cumDistNegBetaStrand;
    kapVal = bisection(func,prob,nGridSize,xMin,xMax,kapValPrev);
    prob = cumXGYNegBetaStrand(tauValPrev,kapValPrev);
    prob = prob + delProb*prob;
    func2 = &randomMol::cumXGYNegBetaStrand;
    tauVal = bisectionDoubleArg(func2,prob,nGridSize,yMin,yMax,kapVal,tauValPrev);
  }else if(type==9){
    prob = cumDistPosBetaStrand(kapValPrev);
    prob = prob + delProb*prob;
    func = &randomMol::cumDistPosBetaStrand;
    kapVal = bisection(func,prob,nGridSize,xMin,xMax,kapValPrev);
    prob = cumXGYPosBetaStrand(tauValPrev,kapValPrev);
    prob = prob + delProb*prob;
    func2 = &randomMol::cumXGYPosBetaStrand;
    tauVal = bisectionDoubleArg(func2,prob,nGridSize,yMin,yMax,kapVal,tauValPrev);
  }else if(type==10){
    prob = cumDistAlphaStrand(kapValPrev);
    //std::cout<<"in 10 "<<prob<<"\n";
    prob = prob + delProb*prob;
    // std::cout<<"in 10 "<<prob<<"\n";
    func = &randomMol::cumDistAlphaStrand;
    //std::cout<<"in 10 "<<kapValPrev<<"\n";
    kapVal = bisection(func,prob,nGridSize,xMin,xMax,kapValPrev);
    //std::cout<<"in 10 "<<kapVal<<"\n";
    prob = cumXGYAlphaStrand(tauValPrev,kapValPrev);
    prob = prob + delProb*prob;
    func2 = &randomMol::cumXGYAlphaStrand;
    tauVal = bisectionDoubleArg(func2,prob,nGridSize,yMin,yMax,kapVal,tauValPrev);
  }else if(type==11){
    prob = cumDistStrandToHelix(kapValPrev);
    prob = prob + delProb*prob;
    func = &randomMol::cumDistStrandToHelix;
    kapVal = bisection(func,prob,nGridSize,xMin,xMax,kapValPrev);
    prob = cumXGYStrandToHelix(tauValPrev,kapValPrev);
    prob = prob + delProb*prob;
    func2 = &randomMol::cumXGYStrandToHelix;
    tauVal = bisectionDoubleArg(func2,prob,nGridSize,yMin,yMax,kapVal,tauValPrev);
  }else if(type==12){
    prob = cumDistHelixToStrand(kapValPrev);
    prob = prob + delProb*prob;
    func = &randomMol::cumDistHelixToStrand;
    kapVal = bisection(func,prob,nGridSize,xMin,xMax,kapValPrev);
    prob = cumXGYHelixToStrand(tauValPrev,kapValPrev);
    prob = prob + delProb*prob;
    func2 = &randomMol::cumXGYHelixToStrand;
    tauVal = bisectionDoubleArg(func2,prob,nGridSize,yMin,yMax,kapVal,tauValPrev);
  }else if(type==13){
    prob = cumDistLinkerToStrand(kapValPrev);
    prob = prob + delProb*prob;
    func = &randomMol::cumDistLinkerToStrand;
    kapVal = bisection(func,prob,nGridSize,xMin,xMax,kapValPrev);
    prob = cumXGYLinkerToStrand(tauValPrev,kapValPrev);
    prob = prob + delProb*prob;
    func2 = &randomMol::cumXGYLinkerToStrand;
    tauVal = bisectionDoubleArg(func2,prob,nGridSize,yMin,yMax,kapVal,tauValPrev);
  }else{
    prob = cumDistStrandToLinker(kapValPrev);
    prob = prob + delProb*prob;
    func = &randomMol::cumDistStrandToLinker;
    kapVal = bisection(func,prob,nGridSize,xMin,xMax,kapValPrev);
    prob = cumXGYStrandToLinker(tauValPrev,kapValPrev);
    prob =  prob + delProb*prob;
    func2 = &randomMol::cumXGYStrandToLinker;
    tauVal = bisectionDoubleArg(func2,prob,nGridSize,yMin,yMax,kapVal,tauValPrev);
  }
  ktpair.first = kapVal;ktpair.second=tauVal;
  return ktpair;
}



point randomMol::parallelTransport(point &tan1,point &tan2,point &norm1){
  point bvec = binorm(tan1,tan2);
  bvec.normalise();
  if(std::abs(tan1.dotprod(tan2))<0.9999999){
    double costhe = tan1.dotprod(tan2);
    double sinthe = std::sqrt(1-costhe*costhe);
    return norm1*costhe + (bvec.cross(norm1))*sinthe + bvec*(bvec.dotprod(norm1))*(1.0-costhe);
  }else{
    return norm1;
  }
}

std::vector<point> randomMol::getTanForHelix(double &k,double &t,double &l,double &xv,double &yv,double &zv,point &p){
  std::vector<point> frame;
  double ktsqrd = k*k+t*t;double ktsqrdRt = std::sqrt(ktsqrd);
  double a = p.getZ()-zv;
  double b=(l*t*t/ktsqrd +k*k*std::sin(l*ktsqrdRt)/std::pow(ktsqrd,1.5));
  double c = k*(-1.0+ std::cos(l*ktsqrdRt))/ktsqrd;
  double the1;
  if(b*b-a*a+c*c>0.0){
  if(std::abs(c)>0.0000000001){
     double rtTerm = std::sqrt((b*b-a*a+c*c)*c*c);
     double yterm = (b*b*a/(b*b+c*c)-a + (b*rtTerm)/(b*b+c*c))/c;
     the1 = std::atan2(yterm,(-b*a-rtTerm)/(b*b+c*c));
  }else{
    the1 = -std::atan(-a/b);
  }
  double a2 = p.getY()-yv;double a1=p.getX()-xv;
  double b2 = k*t*(l/ktsqrd - std::sin(l*ktsqrdRt)/std::pow(ktsqrd,1.5));
  double c1 = -c*std::cos(the1) + b*std::sin(the1);
  double phi1,phi2;
  double arg = (a1*b2 -a2*c1)/(b2*b2+c1*c1);
  double sarg = std::asin(arg);
  phi1 = 3.14159265359-sarg;phi2 =sarg;
  point tan1(std::sin(the1)*std::cos(phi1),std::sin(the1)*std::sin(phi1),std::cos(the1));
  point norm1(std::cos(the1)*std::cos(phi1),std::cos(the1)*std::sin(phi1),-std::sin(the1));
  point tan2(std::sin(the1)*std::cos(phi2),std::sin(the1)*std::sin(phi2),std::cos(the1));
  point norm2(std::cos(the1)*std::cos(phi2),std::cos(the1)*std::sin(phi2),-std::sin(the1));
  point binorm1 = tan1.cross(norm1);point binorm2 = tan2.cross(norm2);
  polyHelix ph;
  ph.updateYvecGeneral(k,t,tan1,norm1,binorm1,p,l);
  point p1 = ph.getCoord();
  ph.updateYvecGeneral(k,t,tan2,norm2,binorm2,p,l);
  point p2 = ph.getCoord();
  point pdesired(xv,yv,zv);
  double ed1 = p1.eDist(pdesired);double ed2 = p2.eDist(pdesired);
  if(ed1<ed2){
    frame.push_back(tan1);frame.push_back(norm1);frame.push_back(binorm1);
  }else{
    frame.push_back(tan2);frame.push_back(norm2);frame.push_back(binorm2);
  }
  return frame;
  }else{
    point tan1(10000.0,10000.0,10000.0);point norm1(10000.0,10000.0,10000.0);
    point norm2(10000.0,10000.0,10000.0);
    frame.push_back(tan1);frame.push_back(norm1);frame.push_back(norm2);
    return frame;
  }
}



std::vector<point> randomMol::makeRandomSectionWithDist(point &stPt,point &stTan,point &stNorm,int &npts,std::default_random_engine &generator,int &strandOrLinker,bool &suceeded){
  // decide the type of section
  std::vector<point> ptsOut;
  double typeProb = distribution(generator);int type;
  if(strandOrLinker==0){
    if(typeProb<=pAlpha){
      type =1;
    }else if(pAlpha<typeProb<=pAlpha+pposBeta){
      type =2;
    }else if(pAlpha+pposBeta<typeProb<=pAlpha+pposBeta+pnegBeta){
      type =3;
    }else{
      type=4;
    }
  }else{
    if(typeProb<=pAlphaStrand){
      type =8;
    }else if(pAlphaStrand<typeProb<=pAlphaStrand+pposBetaStrand){
      type =8;
    }else if(pAlphaStrand+pposBetaStrand<typeProb<=pAlphaStrand+pposBetaStrand+pnegBetaStrand){
      type =9;
    }else{
      type=9;
    }
  }
  //get the first three lengths (this allows us the set their theta values
  std::uniform_real_distribution<double> rDistribution(rmin,rmax);
  double l1=rDistribution(generator);double l2 =rDistribution(generator);double l3 =rDistribution(generator);
  point stBinorm = stTan.cross(stNorm);
  point pt1,pt2,pt3,pt4;
  pt1=stPt;
  pt2 = stPt + stTan*l1;
  double theMax1 =3.14159265359-std::acos((l1*l1+l2*l2-lmin*lmin)/(2.0*l1*l2));
  double theMax2 =3.14159265359-std::acos((l2*l2+l3*l3-lmin*lmin)/(2.0*l2*l3));
  std::uniform_real_distribution<double> theDist1(0,theMax1);
  std::uniform_real_distribution<double> theDist2(0,theMax2);
  std::uniform_real_distribution<double> phiDist(0.0,6.28318530718);
  point tanCurr;
  bool valid=false;
  if(npts>3){
    // lmin is the minimum allowed distance of non neigbouring c alpahs
    point tan1,tan2,tan3,tan4,norm2,norm3,norm4,binrom2,binorm3,binorm4,currnorm;
    int kv=0;
    while(valid==false &&kv<100){
      double the1 = theDist1(generator);double the2=theDist2(generator);
      double phi1 = phiDist(generator);double phi2=phiDist(generator);
      tan2 = stTan*std::sin(the1)+stNorm*std::cos(the1)*std::cos(phi1)+stBinorm*std::cos(the1)*std::sin(phi1);
      point norm2 = parallelTransport(stTan,tan2,stNorm);
      point binorm2 = tan2.cross(norm2);
      pt3 = pt2 +tan2*l2;
      std::pair<double,double> ktpair =getKapTauProb(type,generator);
      pt4 =getNextPoint(stPt,pt2,pt3,l3,ktpair.first,ktpair.second);
      if(pt4.getX()!=10000.0){
	// apply the non local check
	double t4 =pt1.eDist(pt3);double t5 =pt1.eDist(pt4);double t6 =pt2.eDist(pt4);
	if(t4>lmin && t5>lmin && t6>lmin){
	  valid = true;
	  ptsOut.push_back(stPt);ptsOut.push_back(pt2);ptsOut.push_back(pt3);ptsOut.push_back(pt4);
	}
	suceeded=true;
      }
      kv++;
    }
    if(suceeded==true){
      for(int i=2;i<=npts-3;i++){
	valid=false;
	pt1=ptsOut[ptsOut.size()-3];pt2=ptsOut[ptsOut.size()-2];pt3=ptsOut[ptsOut.size()-1];
	kv=0;
	if(suceeded==true){
	  while(valid==false&&kv<100){
	    kv++;
	    std::pair<double,double> ktpair =getKapTauProb(type,generator);
	    l3 =rDistribution(generator);
	    pt4 =getNextPoint(pt1,pt2,pt3,l3,ktpair.first,ktpair.second);
	    if(pt4.getX()!=10000.0){
	      std::vector<double> dsts; 
	      for(int k=0;k<ptsOut.size()-1;k++){
		dsts.push_back(pt4.eDist(ptsOut[k]));
	      }
	      if(*std::min_element(std::begin(dsts),std::end(dsts))>lmin){
		valid = true;
		ptsOut.push_back(pt4);
		suceeded=true;
	      }else{
		suceeded=false;
	      }
	    }
	  }
	}else{
	  pt4.setX(10000.0);pt4.setY(10000.0);pt4.setZ(10000.0);
	  ptsOut.push_back(pt4);
	}
      }
    }else{
      // here we have not been able to find a new next point
      suceeded=false;
      for(int i=2;i<=npts-3;i++){
	pt4.setX(10000.0);pt4.setY(10000.0);pt4.setZ(10000);
	ptsOut.push_back(pt4);
      }
    }
  }else{
    double the1 = theDist1(generator);double the2=theDist2(generator);
    double phi1 = phiDist(generator);double phi2=phiDist(generator);
    point tan2 = stTan*std::sin(the1)+stNorm*std::cos(the1)*std::cos(phi1)+stBinorm*std::cos(the1)*std::sin(phi1);
    point norm2 = parallelTransport(stTan,tan2,stNorm);
    point binorm2 = tan2.cross(norm2);
    pt3 = pt2 +tan2*l2;
    if(npts==2){
      ptsOut.push_back(pt1);
      ptsOut.push_back(pt2);
    }else if(npts==3){
      ptsOut.push_back(pt1);
      ptsOut.push_back(pt2);
      ptsOut.push_back(pt3);
    }else{
      ptsOut.push_back(pt1);
    }
  }
  return ptsOut;
}
  
std::vector<point> randomMol::alterSectionWithDist(std::vector<point> currPts,int &strandOrLinker,double radOfSearch,std::default_random_engine &generator,bool &suceeded){
  // decide the type of section
  double typeProb = distribution(generator);int type;
  if(strandOrLinker==0){
    if(typeProb<=pAlpha){
      type =1;
    }else if(pAlpha<typeProb<=pAlpha+pposBeta){
      type =2;
    }else if(pAlpha+pposBeta<typeProb<=pAlpha+pposBeta+pnegBeta){
      type =3;
    }else{
      type=4;
    }
  }else{
    if(typeProb<=pAlphaStrand){
      type =8;
    }else if(pAlphaStrand<typeProb<=pAlphaStrand+pposBetaStrand){
      type =8;
    }else if(pAlphaStrand+pposBetaStrand<typeProb<=pAlphaStrand+pposBetaStrand+pnegBetaStrand){
      type =9;
    }else{
      type=9;
    }
  }
  //std::cout<<"in random change alg "<<type<<"\n";
  int npts = currPts.size();double l3;
  /*for(int i=0;i<npts;i++){
    currPts[i].printPoint();
    }*/
  std::vector<point> ptsOut;
  //use the original first three points
  std::uniform_real_distribution<double> rDistribution(rmin,rmax);
  double newK,newT;
  //std::cout<<"number of points "<<npts<<" "<<currPts.size()<<"\n";
  if(npts>3){
    //use the original first three points
     point pt1,pt2,pt3,pt4;
     pt1 =currPts[0];pt2 = currPts[1];pt3 = currPts[2];
     double l3;
     point tanCurr;
     bool valid=false;
    // lmin is the minimum allowed distance of non neigbouring c alpahs
    point tan1,tan2,tan3,tan4,norm2,norm3,norm4,binrom2,binorm3,binorm4,currnorm;
    int kv=0;
    std::pair<double,double> ktpair = kapTau(pt1,pt2,pt3,currPts[3]);
    suceeded=false;
    while(valid==false &&kv<100){
      // get the old kt pair
      //std::cout<<kv<<" "<<4000<<" in loop "<<type<<" "<<ktpair.first<<" "<<ktpair.second<<"\n";
      std::pair<double,double> ktpairNew = getProbFromKapTau(type,generator,ktpair,radOfSearch);
      //std::cout<<"new pair "<<ktpairNew.first<<" "<<ktpairNew.second<<"\n";
      l3 =rDistribution(generator);
      newK = ktpairNew.first;newT=ktpairNew.second;
      pt4 =getNextPoint(pt1,pt2,pt3,l3,newK,newT);
      if(pt4.getX()!=10000.0){
	// apply the non local check
	double t4 =pt1.eDist(pt3);double t5 =pt1.eDist(pt4);double t6 =pt2.eDist(pt4);
	if(t4>lmin && t5>lmin && t6>lmin){
	  valid = true;
	  ptsOut.push_back(pt1);ptsOut.push_back(pt2);ptsOut.push_back(pt3);ptsOut.push_back(pt4);
	  suceeded=true;
	}
      }
      kv++;
    }
    if(suceeded==true){
      for(int i=2;i<=currPts.size()-3;i++){
	if(suceeded==true){
	  valid=false;
	  pt1=ptsOut[ptsOut.size()-3];pt2=ptsOut[ptsOut.size()-2];pt3=ptsOut[ptsOut.size()-1];pt4=currPts[ptsOut.size()];
	  kv=0;
	  while(valid==false&&kv<100){
	    kv++;
	    std::pair<double,double> ktpair =kapTau(pt1,pt2,pt3,pt4);
	    ktpair = getProbFromKapTau(type,generator,ktpair,radOfSearch);
	    l3 =rDistribution(generator);
	    newK = ktpair.first;newT=ktpair.second;
	    pt4 =getNextPoint(pt1,pt2,pt3,l3,newK,newT);
	    if(pt4.getX()!=10000.0){
	      std::vector<double> dsts; 
	      for(int k=0;k<ptsOut.size()-1;k++){
		dsts.push_back(pt4.eDist(ptsOut[k]));
	      }
	      if(*std::min_element(std::begin(dsts),std::end(dsts))>lmin){
		valid = true;
		ptsOut.push_back(pt4);
		suceeded=true;
	      }else{
		suceeded=false;
	      }
	    }else{
	      suceeded=false;
	    }
	  }
	}
      }
    }else{
      // here we have not been able to find a new net point retreat back to the old points
      //suceeded=false;
      ptsOut=currPts;
    }
  }else{
    suceeded=true;
    if(npts==2){
      ptsOut.push_back(currPts[0]);
      ptsOut.push_back(currPts[1]);
    }else if(npts==3){
      ptsOut.push_back(currPts[0]);
      ptsOut.push_back(currPts[1]);
      ptsOut.push_back(currPts[2]);
    }else{
      ptsOut.push_back(currPts[0]);
    }
  }
  if(suceeded==false){
    ptsOut=currPts;
  }
  /*for(int i=0;i<npts;i++){
    ptsOut[i].printPoint();
    }*/
  suceeded=true;
  return ptsOut;
}
   

  

std::vector<point> randomMol::blendLoopToHelix(point &ptloop1,point &ptloop2,point &ptloop3,double k,double t,double al,double cl,int &nHel,std::default_random_engine &generator,bool &suceeded,int &index){
  std::vector<point> ptsOut;std::vector<point> tanOut;
  std::vector<point> normOut;std::vector<point> binormOut;
  std::uniform_real_distribution<double> rDistribution(rmin,rmax);
  bool valid=false;
  int type = 5;
  point pt4;
  int kv=0;
  std::vector<point> frame;
  std::vector<point> dummyframe;
  point tan1(10000.0,10000.0,10000.0);point norm1(10000.0,10000.0,10000.0);
  point norm2(10000.0,10000.0,10000.0);
  dummyframe.push_back(tan1);dummyframe.push_back(norm1);dummyframe.push_back(norm2);
  while(valid ==false&&kv<1000){
    kv++;
    std::pair<double,double> ktpair =getKapTauProb(type,generator);
    //std::cout<<"loop to helix "<<kv<<" "<<ktpair.first<<" "<<ktpair.second<<"\n";
    //pt4 =getNextPointTest(ptloop1,ptloop2,ptloop3,cl,ktpair.first,ktpair.second);
    pt4 =getNextPoint(ptloop1,ptloop2,ptloop3,cl,ktpair.first,ktpair.second);
    if(pt4.getX()!=10000.0){
      double t1 = pt4.eDist(ptloop1);double t2 = pt4.eDist(ptloop2);
      if(t1>lmin && t2>lmin){
	double xv =pt4.getX();double yv =pt4.getY();double zv =pt4.getZ();
	frame = getTanForHelix(k,t,al,xv,yv,zv,ptloop3);
	if(frame[0].getX()!=10000.0){
	  valid = true;
	  suceeded=true;
	}
      }
    }else{
      suceeded=false;
      frame = dummyframe;
    }
  }
  polyHelix ph;
  for(int i=1;i<=nHel;i++){
    double len = al*i;
    ph.updateYvecGeneral(k,t,frame[0],frame[1],frame[2],ptloop3,len);
    ptsOut.push_back(ph.getCoord());
  }
  //frameSet[index] = ph.getFrame();
  return ptsOut;
}

std::vector<point> randomMol::blendLoopToHelixAlter(point &ptloop1,point &ptloop2,point &ptloop3,point &ptloop4,double k,double t,double al,double cl,int &nHel,std::default_random_engine &generator,double radOfSearch,bool &suceeded,int &index){
  std::vector<point> ptsOut;std::vector<point> tanOut;
  std::vector<point> normOut;std::vector<point> binormOut;
  std::uniform_real_distribution<double> rDistribution(rmin,rmax);
  bool valid=false;
  int type = 5;
  point pt4;
  int kv=0;
  std::vector<point> frame;
  std::vector<point> dummyframe;
  point tan1(10000.0,10000.0,10000.0);point norm1(10000.0,10000.0,10000.0);
  point norm2(10000.0,10000.0,10000.0);
  dummyframe.push_back(tan1);dummyframe.push_back(norm1);dummyframe.push_back(norm2);
  bool suc=false;
  while(valid ==false&&kv<100){
    kv++;
    std::pair<double,double> ktpair = kapTau(ptloop1,ptloop2,ptloop3,ptloop4);
    //std::cout<<ktpair.first<<" "<<ktpair.second<<"\n";
    ktpair = getProbFromKapTau(type,generator,ktpair,radOfSearch);
    //std::cout<<ktpair.first<<" "<<ktpair.second<<"\n";
     double newK = ktpair.first;double newT=ktpair.second; 
    pt4 =getNextPoint(ptloop1,ptloop2,ptloop3,cl,newK,newT);
    //pt4.printPoint();
    if(pt4.getX()!=10000.0){
      double t1 = pt4.eDist(ptloop1);double t2 = pt4.eDist(ptloop2);
      //std::cout<<t1<<" "<<t2<<" "<<lmin<<"\n";
      if(t1>lmin && t2>lmin){
	double xv =pt4.getX();double yv =pt4.getY();double zv =pt4.getZ();
	frame = getTanForHelix(k,t,al,xv,yv,zv,ptloop3);
	if(frame[0].getX()!=10000.0){
	  valid = true;
	  suc=true;
	}
      }
    }else{
      suc=false;
      frame = dummyframe;
    }
  }
  //std::cout<<"step1 done "<<frame.size()<<"\n";
  if(suc==false){
    //std::cout<<"here in faliure ?\n";
    double xv =ptloop4.getX();double yv =ptloop4.getY();double zv =ptloop4.getZ();
    frame = getTanForHelix(k,t,al,xv,yv,zv,ptloop3);
    pt4=ptloop4;
  }
  polyHelix ph;
  for(int i=1;i<=nHel;i++){
    double len = al*i;
    ph.updateYvecGeneral(k,t,frame[0],frame[1],frame[2],ptloop3,len);
    ptsOut.push_back(ph.getCoord());
  }
  // frameSet[index] = ph.getFrame();
  suceeded=true;
  return ptsOut;
}


std::vector<point> randomMol::blendLoopWithStrand(point &ptloop1,point &ptloop2,point &ptloop3,int &nHel,std::default_random_engine &generator,int strandOrLoop,bool &suceeded){
  std::vector<point> ptsOut;
  std::uniform_real_distribution<double> rDistribution(rmin,rmax);
  bool valid=false;
  int type;
  if(strandOrLoop ==0){
    // strand to loop
    type = 14;
  }else{
    type = 13; 
  }
  point pt4;
  int kv=0;
  std::vector<point> frame;
  while(valid ==false &&kv<1000){
    kv++;
    std::pair<double,double> ktpair =getKapTauProb(type,generator);
    //std::cout<<ktpair.first<<" "<<ktpair.second<<"\n";
    double l3 =rDistribution(generator);
    pt4 =getNextPoint(ptloop1,ptloop2,ptloop3,l3,ktpair.first,ktpair.second);
    if(pt4.getX()!=10000.0){
      double t1 = pt4.eDist(ptloop1);double t2 = pt4.eDist(ptloop2);
      if(t1>lmin && t2>lmin){
	  valid = true;
	  suceeded=true;
      }else{
	suceeded=false;
      }
    }
  }
  point tan1 = pt4-ptloop3;
  double l1=tan1.length();double l2 =rDistribution(generator);
  double theMax1 =3.14159265359-std::acos((l1*l1+l2*l2-lmin*lmin)/(2.0*l1*l2));;
  std::uniform_real_distribution<double> theDist1(0,theMax1);
  std::uniform_real_distribution<double> phiDist(0.0,6.28318530718);
  double the1 = theDist1(generator);double phi1 = phiDist(generator);
  tan1.normalise();
  point vert(0.0,0.0,1.0);
  point norm1;
  if(std::abs(tan1.dotprod(vert))>0.0001){
    point vert(0.0,0.0,1.0);
    norm1 = vert.cross(tan1);
    norm1.normalise();
  }else{
    point vert(0.0,1.0,0.0);
    norm1 = vert.cross(tan1);
    norm1.normalise();
  }
  point binorm1 = tan1.cross(norm1);
  binorm1.normalise();
  point tan2 = tan1*std::sin(the1)+norm1*std::cos(the1)*std::cos(phi1)+binorm1*std::cos(the1)*std::sin(phi1);
  point norm2 = parallelTransport(tan1,tan2,norm1);
  if(suceeded==true){
    ptsOut= makeRandomSectionWithDist(pt4,tan2,norm2,nHel,generator,strandOrLoop,suceeded);
  }else{
    for(int i=0;i<nHel;i++){
      point p(10000.0,10000.0,10000.0);
      ptsOut.push_back(p);
    }
  }
  return ptsOut;
}

std::vector<point>  randomMol::getFrameHelix(std::vector<point> &coords){
  std::vector<point> outputFrame;
  int sz =coords.size();
  std::vector<point> lines;point mean(0.0,0.0,0.0);
  for(int j=0;j<sz;j++){
    mean = mean + coords[j];
  }
  mean = mean*(1.0/double(sz));
  point endDif = coords[sz-1]-mean;
  for(int j=0;j<sz;j++){
    lines.push_back(mean - coords[j]);
    if(lines[j].dotprod(endDif)<0.0){
      lines[j] = lines[j]*(-1.0);
    }
  }
  point direc(0.0,0.0,0.0);
  for(int j=0;j<sz;j++){
    direc = direc + lines[j];
  }
  double hlfDist= 0.5*(mean.eDist(coords[0])+mean.eDist(coords[sz-1]));
  direc.normalise();
  outputFrame.push_back(direc);
  point endPt = mean+direc*hlfDist;
  point approxNorm = coords[coords.size()-1]-endPt;
  approxNorm.normalise();
  //std::cout<<"norm comp \n";
  //approxNorm.printPoint();
  point norm = approxNorm - direc*(approxNorm.dotprod(direc));
  norm.normalise();
  //norm.printPoint();
  outputFrame.push_back(norm);
  point binorm = direc.cross(norm);
  outputFrame.push_back(binorm);
  return outputFrame;
}

std::vector<point> randomMol::blendHelixToHelix(std::vector<point> &prevSec,double k,double t,double al,double cl,int &nHel,std::default_random_engine &generator,bool &suceeded,int &index){
  std::vector<point> ptsOut;
  /*std::uniform_real_distribution<double> rDistribution(rmin,rmax);
  bool valid=false;
  int type = 10;
  point pt4;
  std::vector<point> frame;
  std::vector<point> dummyframe;
  point tan1(10000.0,10000.0,10000.0);point norm1(10000.0,10000.0,10000.0);
  point norm2(10000.0,10000.0,10000.0);
  dummyframe.push_back(tan1);dummyframe.push_back(norm1);dummyframe.push_back(norm2);
  std::uniform_real_distribution<double> kapdel(-0.005,0.005);
  std::uniform_real_distribution<double> taudel(-0.005,0.005);
  int kv=0;
  while(valid ==false&&kv<1000){
    kv++;
    //std::pair<double,double> ktpair =getKapTauProb(type,generator);
    std::pair<double,double> ktpair;
    ktpair.first = 0.47 + kapdel(generator);
    ktpair.second = 0.23 + taudel(generator);
    //pt4 =getNextPointTest(ptloop1,ptloop2,ptloop3,cl,ktpair.first,ktpair.second);
    pt4 =getNextPoint(ptloop1,ptloop2,ptloop3,cl,ktpair.first,ktpair.second);
    if(pt4.getX()!=10000.0){
      double t1 = pt4.eDist(ptloop1);double t2 = pt4.eDist(ptloop2);
      if(t1>lmin && t2>lmin){
      	double xv =pt4.getX();double yv =pt4.getY();double zv =pt4.getZ();
      	frame = getTanForHelix(k,t,al,xv,yv,zv,ptloop3);
    	if(frame[0].getX()!=10000.0){
	  valid = true;
	  suceeded=true;
	}else{
	  suceeded=false;
	}
      }
    }else{
      frame = dummyframe;
    }
    kv++;
    }
  polyHelix ph;
  for(int i=1;i<=nHel;i++){
    double len = al*i;
    ph.updateYvecGeneral(k,t,frame[0],frame[1],frame[2],ptloop3,len);
    ptsOut.push_back(ph.getCoord());
  }
  */
  std::uniform_real_distribution<double> angVar(-0.2,0.2);
  std::uniform_real_distribution<double> tanVar(-0.2,0.2);
  std::vector<point> newFrame = getFrameHelix(prevSec);
  // first rotate notmal and binormal
  double ang = angVar(generator);
  point newNorm = newFrame[1]*std::cos(ang) + newFrame[2]*std::sin(ang);
  point newBinorm =newFrame[1]*(-1.0)*std::sin(ang) + newFrame[2]*std::cos(ang);
  newFrame[1] = newNorm;
  newFrame[2] = newBinorm;
  double xtanCh = newFrame[0].getX()+tanVar(generator);
  double ytanCh = newFrame[0].getY()+tanVar(generator);
  double ztanCh = newFrame[0].getZ()+tanVar(generator);
  point newTan(xtanCh,ytanCh,ztanCh);
  newTan.normalise();
  // have to make the norm and binorm vectors normal to this, firts prject the normal
  point normProj = newNorm - newTan*newNorm.dotprod(newTan);
  normProj.normalise();
  newBinorm = newTan.cross(normProj);
  newBinorm.normalise();
  newFrame[0]=newTan;
  newFrame[1]=normProj;
  newFrame[2]=newBinorm;
  double sqprod = k*k+t*t;
  double helR = k/sqprod;
  double helc = t/sqprod;
  double delt = al/std::sqrt(helR*helR + helc*helc);
  point stpt = prevSec[prevSec.size()-1]-newFrame[1]*helR;
  for(int i=1;i<=nHel;i++){
   double len = delt*i;
   if(t>0.0){
     // right handed
     point newnorm = (newFrame[1]*std::cos(len)+newFrame[2]*std::sin(len))*helR;
     ptsOut.push_back(stpt  + newFrame[0]*helc*len + newnorm);
   }else{
     // left handed
     point newnorm = (newFrame[1]*std::cos(-len)+newFrame[2]*std::sin(-len))*helR;
     ptsOut.push_back(stpt  + newFrame[0]*helc*len + newnorm);
   }
  }
  //frameSet[index] = ph.getFrame();
  return ptsOut;
}
 
std::vector<point> randomMol::blendHelixToHelixAlter(point &ptloop1,point &ptloop2,point &ptloop3,point &ptloop4,double k,double t,double al,double cl,int &nHel,std::default_random_engine &generator,double radOfSearch,bool &suceeded,int &index){
  std::vector<point> ptsOut;
  std::uniform_real_distribution<double> rDistribution(rmin,rmax);
  bool valid=false;
  int type = 6;
  point pt4;
  std::vector<point> frame;
  std::vector<point> dummyframe;
  point tan1(10000.0,10000.0,10000.0);point norm1(10000.0,10000.0,10000.0);
  point norm2(10000.0,10000.0,10000.0);
  dummyframe.push_back(tan1);dummyframe.push_back(norm1);dummyframe.push_back(norm2);
  int kv=0;
  while(valid ==false&&kv<100){
    kv++;
    std::pair<double,double> ktpair = kapTau(ptloop1,ptloop2,ptloop3,ptloop4);
    ktpair = getProbFromKapTau(type,generator,ktpair,radOfSearch);
    double l3 =rDistribution(generator);
    double newK = ktpair.first;double newT=ktpair.second;
    pt4 =getNextPoint(ptloop1,ptloop2,ptloop3,cl,newK,newT);
    if(pt4.getX()!=10000.0){
      double t1 = pt4.eDist(ptloop1);double t2 = pt4.eDist(ptloop2);
      if(t1>lmin && t2>lmin){
      	double xv =pt4.getX();double yv =pt4.getY();double zv =pt4.getZ();
      	frame = getTanForHelix(k,t,al,xv,yv,zv,ptloop3);
    	if(frame[0].getX()!=10000.0){
	  valid = true;
	  suceeded=true;
	}else{
	  suceeded=false;
	}
      }
    }else{
      frame = dummyframe;
    }
    kv++;
  }
  if(suceeded==false){
    double xv =ptloop4.getX();double yv =ptloop4.getY();double zv =ptloop4.getZ();
    frame = getTanForHelix(k,t,al,xv,yv,zv,ptloop3);
    pt4=ptloop4;
  }
  polyHelix ph;
  for(int i=1;i<=nHel;i++){
    double len = al*i;
    ph.updateYvecGeneral(k,t,frame[0],frame[1],frame[2],ptloop3,len);
    ptsOut.push_back(ph.getCoord());
  }
  //frameSet[index] = ph.getFrame();
  suceeded=true;
  return ptsOut;
}


std::vector<point> randomMol::blendHelixToLoop(point &pthel1,point &pthel2,point &pthel3,int &nl,std::default_random_engine &generator,int strandOrLoop,bool &suceeded){
  std::uniform_real_distribution<double> rDistribution(rmin,rmax);
  bool valid=false;
  int type;
  std::vector<point> ptsOut;
  if(strandOrLoop==0){
    type = 6;
  }else{
    type = 12;
  }
  point pt4;
  int kv=0;
  while(valid ==false &&kv<1000){
    kv++;
    std::pair<double,double> ktpair =getKapTauProb(type,generator);
    double l3 =rDistribution(generator);
    // std::cout<<l3<<" "<<ktpair.first<<" "<<ktpair.second<<" "<<type<<"\n";
    pt4 =getNextPoint(pthel1,pthel2,pthel3,l3,ktpair.first,ktpair.second);
    //pt4.printPoint();
    if(pt4.getX()!=10000.0){
      double t1 = pt4.eDist(pthel1);double t2 = pt4.eDist(pthel2);
      if(t1>lmin && t2>lmin){
	valid = true;
	suceeded=true;
      }else{
	suceeded=false;
      }
    }
    kv++;
  }
  point tan1 = pt4-pthel3;
  double l1=tan1.length();double l2 =rDistribution(generator);
  double theMax1 =3.14159265359-std::acos((l1*l1+l2*l2-lmin*lmin)/(2.0*l1*l2));;
  std::uniform_real_distribution<double> theDist1(0,theMax1);
  std::uniform_real_distribution<double> phiDist(0.0,6.28318530718);
  double the1 = theDist1(generator);double phi1 = phiDist(generator);
  tan1.normalise();
  point vert(0.0,0.0,1.0);
  point norm1;
  if(std::abs(tan1.dotprod(vert))>0.0001){
    point vert(0.0,0.0,1.0);
    norm1 = vert.cross(tan1);
    norm1.normalise();
  }else{
    point vert(0.0,1.0,0.0);
    norm1 = vert.cross(tan1);
    norm1.normalise();
  }
  point binorm1 = tan1.cross(norm1);
  binorm1.normalise();
  point tan2 = tan1*std::sin(the1)+norm1*std::cos(the1)*std::cos(phi1)+binorm1*std::cos(the1)*std::sin(phi1);
  point norm2 = parallelTransport(tan1,tan2,norm1);
  if(suceeded==true){
    ptsOut= makeRandomSectionWithDist(pt4,tan2,norm2,nl,generator,strandOrLoop,suceeded);
    }else{
    for(int i=0;i<nl;i++){
      point p(10000.0,10000.0,10000.0);
      ptsOut.push_back(p);
    }
  }
  return ptsOut;
}


 std::vector<std::vector<point> >  randomMol::makeRandomMolecule(std::vector<std::pair<std::string,int> > &molDat,point &sp,bool &suceeded){
  std::random_device rdev{};
  unsigned int seed =std::chrono::steady_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  point origTan(0.0,0.0,1.0);point origNorm(0.0,1.0,0.0);point origBinorm = origTan.cross(origNorm);
  std::vector<std::vector<point> > allPts;
  // make the first section
  polyHelix ph;
  std::vector<point> pts;
  // fill an empty frame set
  /*for(int i=0;i<molDat.size();i++){
     frameSet.push_back(dummyFrame);
     };*/ 
  if(molDat[0].first=="Helix"){
    for(int i=0;i<molDat[0].second;i++){
      double len = alphaAl*i;
      ph.updateYvecGeneral(alphaKap,alphaTau,origTan,origNorm,origBinorm,sp,len);
      pts.push_back(ph.getCoord());
    }
    //frameSet[0]=ph.getFrame();
  }else if(molDat[0].first=="strand"){
    int strandOrLoop =1;
    pts = makeRandomSectionWithDist(sp,origTan,origNorm,molDat[0].second,generator,strandOrLoop,suceeded);
  }else{
    int strandOrLoop =0;
    pts = makeRandomSectionWithDist(sp,origTan,origNorm,molDat[0].second,generator,strandOrLoop,suceeded);
  }
  allPts.push_back(pts);
  // now loop over the rest
  std::string prevType = molDat[0].first;
  std::string type;
  for(int i=1;i<molDat.size();i++){
    type = molDat[i].first;
    //std::cout<<prevType<<" "<<type<<" "<<molDat[i-1].second<<" "<<molDat[i].second<<"\n";
    if(type == prevType){
      if(type=="Helix"){
        std::vector<point> prevSec = allPts[i-1];
        allPts.push_back(blendHelixToHelix(prevSec,alphaKap,alphaTau,alphaAl,alphaCl,molDat[i].second,generator,suceeded,i));
      }else{
        std::vector<point> prevSec = allPts[i-1];
        allPts.push_back(blendHelixToHelix(prevSec,betaKap,betaTau,betaAl,betaCl,molDat[i].second,generator,suceeded,i));
      }
    }else{
      if(prevType=="loop"&&type=="Helix"){
        std::vector<point> prevSec = allPts[i-1];
        //std::cout<<"in here ?\n";
        if(prevSec.size()==2){
          std::vector<point> prevPrevSec = allPts[i-2];
          allPts.push_back(blendLoopToHelix(prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],alphaKap,alphaTau,alphaAl,alphaCl,molDat[i].second,generator,suceeded,i));
        }else if(prevSec.size()==1){
          std::vector<point> prevPrevSec = allPts[i-2];
          if(prevPrevSec.size()==1){
          std::cout<<"problem two length 1 sections in a row\n";
          }
          //std::cout<<"here\n";
          allPts.push_back(blendLoopToHelix(prevPrevSec[prevPrevSec.size()-2],prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-1],alphaKap,alphaTau,alphaAl,alphaCl,molDat[i].second,generator,suceeded,i));
        }
        else{
          allPts.push_back(blendLoopToHelix(prevSec[prevSec.size()-3],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],alphaKap,alphaTau,alphaAl,alphaCl,molDat[i].second,generator,suceeded,i));
        }
      }else if(prevType=="loop"&&type=="Strand"){
        std::vector<point> prevSec = allPts[i-1];
        if(prevSec.size()==2){
          std::vector<point> prevPrevSec = allPts[i-2];
	  int loopOrStrand=1;
          allPts.push_back(blendLoopWithStrand(prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],molDat[i].second,generator,loopOrStrand,suceeded));
        }else if(prevSec.size()==1){
          std::vector<point> prevPrevSec = allPts[i-2];
	  int loopOrStrand=1;
           if(prevPrevSec.size()==1){
          std::cout<<"problem two length 1 sections in a row\n";
          }
	   allPts.push_back(blendLoopWithStrand(prevPrevSec[prevPrevSec.size()-2],prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-1],molDat[i].second,generator,loopOrStrand,suceeded));
        }
        else{
	  int loopOrStrand=1;
          allPts.push_back(blendLoopWithStrand(prevSec[prevSec.size()-3],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],molDat[i].second,generator,loopOrStrand,suceeded));
        }
      }
      else if(prevType=="Strand"&&type=="loop"){
        std::vector<point> prevSec = allPts[i-1];
        if(prevSec.size()==2){
          std::vector<point> prevPrevSec = allPts[i-2];
	  int loopOrStrand=0;
          allPts.push_back(blendLoopWithStrand(prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],molDat[i].second,generator,loopOrStrand,suceeded));
        }else if(prevSec.size()==1){
          std::vector<point> prevPrevSec = allPts[i-2];
          if(prevPrevSec.size()==1){
          std::cout<<"problem two length 1 sections in a row\n";
          }
	  int loopOrStrand=0;
          allPts.push_back(blendLoopWithStrand(prevPrevSec[prevPrevSec.size()-2],prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-1],molDat[i].second,generator,loopOrStrand,suceeded));
        }
        else{
	  int loopOrStrand=0;
          allPts.push_back(blendLoopWithStrand(prevSec[prevSec.size()-3],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],molDat[i].second,generator,loopOrStrand,suceeded));
        }
      }
      else{
	// Here we have a helix to loop transition
        std::vector<point> prevSec = allPts[i-1];
	if(type=="Strand"){
	  int loopOrStrand=1;
	  if(prevSec.size()>=3){
	    allPts.push_back(blendHelixToLoop(prevSec[prevSec.size()-3],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],molDat[i].second,generator,loopOrStrand,suceeded));
	  }else{
	    std::vector<point> prevPrevSec = allPts[i-2];
	    allPts.push_back(blendHelixToLoop(prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],molDat[i].second,generator,loopOrStrand,suceeded));
	  }
	}else{
	  int loopOrStrand=0;
	  if(prevSec.size()>=3){
	    allPts.push_back(blendHelixToLoop(prevSec[prevSec.size()-3],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],molDat[i].second,generator,loopOrStrand,suceeded));
	  }else{
	    std::vector<point> prevPrevSec = allPts[i-2];
	    allPts.push_back(blendHelixToLoop(prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],molDat[i].second,generator,loopOrStrand,suceeded));
	  }
	}
      }
    }
    prevType=type;
  }
  return allPts;  
}
  
  void randomMol::tanToSec(point &sp,std::vector<point> &prevSec){
  std::vector<point> tempSec;
  tempSec.push_back(sp);
  point newpt;
  for(int i=1;i<prevSec.size();i++){
    if(i==1){
      newpt = sp + (prevSec[i]-prevSec[i-1]);
      tempSec.push_back(newpt);
    }else{
      newpt = newpt + (prevSec[i]-prevSec[i-1]);
      tempSec.push_back(newpt);
    }
  }
  prevSec = tempSec;
}

double randomMol::getMaxDistChange(std::vector<point> &oldPts,std::vector<point> &newPoints){
  double maxDist = 0.0;
  for(int i=0;i<oldPts.size();i++){
    double dist = oldPts[i].eDist(newPoints[i]);
    if(dist>maxDist){
      maxDist = dist;
    }
  }
  return maxDist;
}


 double randomMol::reshapeMol(std::vector<std::pair<std::string,int> > &molDat,std::vector<std::vector<point> > &molPos,int &index,point &sp,bool &suceeded){
  std::random_device rdev{};
  std::default_random_engine generator{rdev()};
  std::vector<point> newSection;
  double maxDist;
  std::string type = molDat[index].first;
  //std::cout<<"section type "<<type<<" "<<molPos[index].size()<<"\n";
  if(index>0){
    std::string prevType = molDat[index-1].first;
    if(type == prevType){
      if(type=="Helix"){
        std::vector<point> prevSec = molPos[index-1];
        newSection = blendHelixToHelix(prevSec,alphaKap,alphaTau,alphaAl,alphaCl,molDat[index].second,generator,suceeded,index);
      }else if(type == "strand"){
        std::vector<point> prevSec = molPos[index-1];
        newSection = blendHelixToHelix(prevSec,betaKap,betaTau,betaAl,betaCl,molDat[index].second,generator,suceeded,index);
      }
    }else{
      if(prevType=="loop"&&type=="Helix"){
        std::vector<point> prevSec = molPos[index-1];
        if(prevSec.size()==2){
          std::vector<point> prevPrevSec = molPos[index-2];
          newSection = blendLoopToHelix(prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],alphaKap,alphaTau,alphaAl,alphaCl,molDat[index].second,generator,suceeded,index);
        }else if(prevSec.size()==1){
          std::vector<point> prevPrevSec = molPos[index-2];
          if(prevPrevSec.size()==1){
          std::cout<<"problem two length 1 sections in a row\n";
          }
          newSection = blendLoopToHelix(prevPrevSec[prevPrevSec.size()-2],prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-1],alphaKap,alphaTau,alphaAl,alphaCl,molDat[index].second,generator,suceeded,index);
        }else{
          newSection = blendLoopToHelix(prevSec[prevSec.size()-3],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],alphaKap,alphaTau,alphaAl,alphaCl,molDat[index].second,generator,suceeded,index);
        }
      }else if(prevType=="loop"&&type=="Strand"){
        std::vector<point> prevSec = molPos[index-1];
        if(prevSec.size()==2){
          std::vector<point> prevPrevSec = molPos[index-2];
	  int loopOrStrand=1;
          newSection=blendLoopWithStrand(prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],molDat[index].second,generator,loopOrStrand,suceeded);
        }else if(prevSec.size()==1){
        std::vector<point> prevPrevSec = molPos[index-2];
	      int loopOrStrand=1;
        if(prevPrevSec.size()==1){
          std::cout<<"problem two length 1 sections in a row\n";
        }
        newSection=blendLoopWithStrand(prevPrevSec[prevPrevSec.size()-2],prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-1],molDat[index].second,generator,loopOrStrand,suceeded);
        }
        else{
	  int loopOrStrand=1;
          newSection=blendLoopWithStrand(prevSec[prevSec.size()-3],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],molDat[index].second,generator,loopOrStrand,suceeded);
        }
      }
      else if(prevType=="Strand"&&type=="loop"){
        std::vector<point> prevSec = molPos[index-1];
        if(prevSec.size()==2){
          std::vector<point> prevPrevSec = molPos[index-2];
	  int loopOrStrand=0;
          newSection=blendLoopWithStrand(prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],molDat[index].second,generator,loopOrStrand,suceeded);
        }else if(prevSec.size()==1){
            std::vector<point> prevPrevSec = molPos[index-2];
            if(prevPrevSec.size()==1){
          std::cout<<"problem two length 1 sections in a row\n";
        }
           int loopOrStrand=0;
	   newSection=blendLoopWithStrand(prevPrevSec[prevPrevSec.size()-2],prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-1],molDat[index].second,generator,loopOrStrand,suceeded);
        }
        else{
	  int loopOrStrand=0;
          newSection=blendLoopWithStrand(prevSec[prevSec.size()-3],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],molDat[index].second,generator,loopOrStrand,suceeded);
        }
      }
      else{
	// Here we have a helix to loop transition
	//	std::cout<<" hel to loop ? \n";
        std::vector<point> prevSec = molPos[index-1];
	if(type=="Strand"){
	  int loopOrStrand=1;
	  if(prevSec.size()>=3){
	    newSection=blendHelixToLoop(prevSec[prevSec.size()-3],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],molDat[index].second,generator,loopOrStrand,suceeded);
	  }else{
	    std::vector<point> prevPrevSec = molPos[index-2];
	    newSection=blendHelixToLoop(prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],molDat[index].second,generator,loopOrStrand,suceeded);
	  }
	}else{
	  int loopOrStrand=0;
	  //std::cout<<"here helix to loop \n";
	  if(prevSec.size()>=3){
	    newSection=blendHelixToLoop(prevSec[prevSec.size()-3],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],molDat[index].second,generator,loopOrStrand,suceeded);
	  }else{
	    std::vector<point> prevPrevSec = molPos[index-2];
	   newSection =blendHelixToLoop(prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],molDat[index].second,generator,loopOrStrand,suceeded);
	  }
	}
      }
    }
  }else{
    // here we are at the first index
      point origTan(0.0,0.0,1.0);point origNorm(0.0,1.0,0.0);point origBinorm = origTan.cross(origNorm);
      polyHelix ph;
      std::vector<point> pts;
      if(molDat[0].first=="Helix"){
	for(int i=0;i<molDat[0].second;i++){
	   double len = alphaAl*i;
	   ph.updateYvecGeneral(alphaKap,alphaTau,origTan,origNorm,origBinorm,sp,len);
	   newSection.push_back(ph.getCoord());
	}
      }else if(molDat[0].first=="Strand"){
	int loopOrStrand = 1; 
	newSection = makeRandomSectionWithDist(sp,origTan,origNorm,molDat[0].second,generator,loopOrStrand,suceeded);
      }else{
	int loopOrStrand = 0;
	newSection = makeRandomSectionWithDist(sp,origTan,origNorm,molDat[0].second,generator,loopOrStrand,suceeded);
      }
  }
  maxDist = getMaxDistChange(molPos[index],newSection);
  point lastmol = molPos[index][molPos[index].size()-1];
  if(suceeded==true){
    molPos[index] = newSection;
  }
  // now we need to find the new positions of the rest of the molecule after index
  //   if the current section is a loop we will need to blend into the next type
  if(suceeded==true&&index<(molDat.size()-1)){
    int newIndex = index+1;
    if(type=="loop" || type=="Strand"){
      std::string nextType = molDat[index+1].first;
      std::vector<point> newNewSection;
      if(nextType=="Helix"){
	std::vector<point> prevSec = newSection;
	if(prevSec.size()==2){
	  std::vector<point> prevPrevSec = molPos[index-1];
	  newNewSection = blendLoopToHelix(prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],alphaKap,alphaTau,alphaAl,alphaCl,molDat[index+1].second,generator,suceeded,newIndex);
	}else if(prevSec.size()==1){
	  std::vector<point> prevPrevSec = molPos[index-1];
	  if(prevPrevSec.size()==1){
	    std::cout<<"problem two length 1 sections in a row\n";
	  }
	  newNewSection = blendLoopToHelix(prevPrevSec[prevPrevSec.size()-2],prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-1],alphaKap,alphaTau,alphaAl,alphaCl,molDat[index+1].second,generator,suceeded,newIndex);
	}
	else{
	  newNewSection = blendLoopToHelix(prevSec[prevSec.size()-3],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],alphaKap,alphaTau,alphaAl,alphaCl,molDat[index+1].second,generator,suceeded,newIndex);
	}
      }else if(nextType=="Strand"){
	std::vector<point> prevSec = newSection;
	if(prevSec.size()==2){
	  std::vector<point> prevPrevSec = molPos[index-1];
	  int loopOrStrand =1;
	  newNewSection=blendLoopWithStrand(prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],molDat[index+1].second,generator,loopOrStrand,suceeded);
	}else if(prevSec.size()==1){
	  std::vector<point> prevPrevSec = molPos[index-1];
	  if(prevPrevSec.size()==1){
	    std::cout<<"problem two length 1 sections in a row\n";
	  }
	  int loopOrStrand =1;
	  newNewSection=blendLoopWithStrand(prevPrevSec[prevPrevSec.size()-2],prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-1],molDat[index+1].second,generator,loopOrStrand,suceeded);
	}
	else{
	  int loopOrStrand =1;
	  newNewSection =blendLoopWithStrand(prevSec[prevSec.size()-3],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],molDat[index+1].second,generator,loopOrStrand,suceeded);
	}
      }else{
	// here it is a strand to loop transistion
	std::vector<point> prevSec = newSection;
	if(prevSec.size()==2){
	  std::vector<point> prevPrevSec = molPos[index-1];
	  int loopOrStrand =0;
	  newNewSection=blendLoopWithStrand(prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],molDat[index+1].second,generator,loopOrStrand,suceeded);
	}else if(prevSec.size()==1){
	  std::vector<point> prevPrevSec = molPos[index-1];
	  if(prevPrevSec.size()==1){
	    std::cout<<"problem two length 1 sections in a row\n";
	  }
	  int loopOrStrand =0;
	  newNewSection=blendLoopWithStrand(prevPrevSec[prevPrevSec.size()-2],prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-1],molDat[index+1].second,generator,loopOrStrand,suceeded);
	}
	else{
	  int loopOrStrand =0;
	  newNewSection =blendLoopWithStrand(prevSec[prevSec.size()-3],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],molDat[index+1].second,generator,loopOrStrand,suceeded);
	}
      }
      double tempMaxDist = getMaxDistChange(molPos[index+1],newNewSection);
      if(tempMaxDist>maxDist){
	maxDist=tempMaxDist;
      }
      lastmol = molPos[index+1][molPos[index+1].size()-1];
      if(suceeded==true){
	molPos[index+1] = newNewSection;
	// now complete the molecule, if it is a helix update the frame state
	for(int j=index+2;j<molDat.size();j++){
	  if(j==index+2){
	    point lp = newNewSection[newNewSection.size()-1] + molPos[j][0]-lastmol;
	    lastmol = molPos[j][molPos[j].size()-1];
	    std::vector<point> prevSec = molPos[j];
	    tanToSec(lp,molPos[j]);
	    double tempMaxDist = getMaxDistChange(prevSec,molPos[j]);
	    if(tempMaxDist>maxDist){
	      maxDist=tempMaxDist;
	    }
	  }else{
	    point lp = molPos[j-1][molPos[j-1].size()-1]+molPos[j][0]-lastmol;
	    lastmol = molPos[j][molPos[j].size()-1];
	    std::vector<point> prevSec = molPos[j];
	    tanToSec(lp,molPos[j]);
	    double tempMaxDist = getMaxDistChange(prevSec,molPos[j]);
	    if(tempMaxDist>maxDist){
	      maxDist=tempMaxDist;
	    }
	  }
	}
      }
    }else{
      for(int j=index+1;j<molDat.size();j++){
	if(j==index+1){
	  point lp = newSection[newSection.size()-1] + molPos[j][0] -lastmol;
	  lastmol = molPos[j][molPos[j].size()-1];
	  std::vector<point> prevSec = molPos[j];
	  tanToSec(lp,molPos[j]);
	  double tempMaxDist = getMaxDistChange(prevSec,molPos[j]);
	  if(tempMaxDist>maxDist){
	    maxDist=tempMaxDist;
	  }
	}else{
	  point lp = molPos[j-1][molPos[j-1].size()-1]+ molPos[j][0] -lastmol;
	  lastmol = molPos[j][molPos[j].size()-1];
	  std::vector<point> prevSec = molPos[j];
	  tanToSec(lp,molPos[j]);
	  double tempMaxDist = getMaxDistChange(prevSec,molPos[j]);
	  if(tempMaxDist>maxDist){
	    maxDist=tempMaxDist;
	  }
	}
      }
    }
  }
  //now update all the relevant Tangents
  return maxDist;
 }


//if we are changing a helix then there is a run of helices


double randomMol::reshapeMolHelixSet(std::vector<std::pair<std::string,int> > &molDat,std::vector<std::vector<point> > &molPos,int &index,int noHelices,point &sp,bool &suceeded){
  std::random_device rdev{};
  std::default_random_engine generator{rdev()};
  std::vector<point> newSection;
  double maxDist;
  std::vector<point> prevSec = molPos[index-1];
  if(prevSec.size()==2){
    std::vector<point> prevPrevSec = molPos[index-2];
    newSection = blendLoopToHelix(prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],alphaKap,alphaTau,alphaAl,alphaCl,molDat[index].second,generator,suceeded,index);
  }else if(prevSec.size()==1){
    std::vector<point> prevPrevSec = molPos[index-2];
    if(prevPrevSec.size()==1){
      std::cout<<"problem two length 1 sections in a row\n";
    }
    newSection = blendLoopToHelix(prevPrevSec[prevPrevSec.size()-2],prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-1],alphaKap,alphaTau,alphaAl,alphaCl,molDat[index].second,generator,suceeded,index);
  }else{
    newSection = blendLoopToHelix(prevSec[prevSec.size()-3],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],alphaKap,alphaTau,alphaAl,alphaCl,molDat[index].second,generator,suceeded,index);
  }  
  maxDist = getMaxDistChange(molPos[index],newSection);
  point lastmol = molPos[index][molPos[index].size()-1];
  if(suceeded==true){
    molPos[index] = newSection;
  }
  // now we need to find the new positions of the rest of the molecule after index
  //   if the current section is a loop we will need to blend into the next type
  if(suceeded==true&&index<(molDat.size()-1)){
    int newIndex;
    std::vector<point> newNewSection;
    for(int l=1;l<=noHelices;l++){
      newIndex = index+l;
      std::vector<point> prevSec = molPos[newIndex-1];
      newNewSection = blendHelixToHelix(prevSec,alphaKap,alphaTau,alphaAl,alphaCl,molDat[newIndex].second,generator,suceeded,index);
      double tempMaxDist = getMaxDistChange(molPos[newIndex],newNewSection);
      if(tempMaxDist>maxDist){
	maxDist=tempMaxDist;
      }
      lastmol = molPos[newIndex][molPos[newIndex].size()-1];
      if(suceeded==true){
	molPos[newIndex] = newNewSection;
      }
    }
    if(suceeded==true){
      // now complete the molecule, if it is a helix update the frame state
      for(int j=newIndex+1;j<molDat.size();j++){
	if(j==newIndex+1){
	  point lp = newNewSection[newNewSection.size()-1] + molPos[j][0]-lastmol;
	  lastmol = molPos[j][molPos[j].size()-1];
	  std::vector<point> prevSec = molPos[j];
	  tanToSec(lp,molPos[j]);
	  double tempMaxDist = getMaxDistChange(prevSec,molPos[j]);
	  if(tempMaxDist>maxDist){
	    maxDist=tempMaxDist;
	  }
	}else{
	  point lp = molPos[j-1][molPos[j-1].size()-1]+molPos[j][0]-lastmol;
	  lastmol = molPos[j][molPos[j].size()-1];
	  std::vector<point> prevSec = molPos[j];
	  tanToSec(lp,molPos[j]);
	  double tempMaxDist = getMaxDistChange(prevSec,molPos[j]);
	  if(tempMaxDist>maxDist){
	    maxDist=tempMaxDist;
	  }
	}
      }
    }
  }
  //now update all the relevant Tangents
  return maxDist;
 }


 double randomMol::reshapeMolLoopThenHelixSet(std::vector<std::pair<std::string,int> > &molDat,std::vector<std::vector<point> > &molPos,int &index,int noHelices,point &sp,bool &suceeded){
  std::random_device rdev{};
  std::default_random_engine generator{rdev()};
  std::vector<point> newSection;
  double maxDist;
  std::string type = molDat[index].first;
  if(index>0){
    std::string prevType = molDat[index-1].first;
    //std::cout<<"index "<<index<<" "<<prevType<<" "<<type<<"\n";
    if(type == prevType){
      if(type=="Helix"){
        std::vector<point> prevSec = molPos[index-1];
        newSection = blendHelixToHelix(prevSec,alphaKap,alphaTau,alphaAl,alphaCl,molDat[index].second,generator,suceeded,index);
      }else if(type == "strand"){
        std::vector<point> prevSec = molPos[index-1];
        newSection = blendHelixToHelix(prevSec,betaKap,betaTau,betaAl,betaCl,molDat[index].second,generator,suceeded,index);
      }
    }else{
      if(prevType=="loop"&&type=="Helix"){
        std::vector<point> prevSec = molPos[index-1];
        if(prevSec.size()==2){
          std::vector<point> prevPrevSec = molPos[index-2];
          newSection = blendLoopToHelix(prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],alphaKap,alphaTau,alphaAl,alphaCl,molDat[index].second,generator,suceeded,index);
        }else if(prevSec.size()==1){
          std::vector<point> prevPrevSec = molPos[index-2];
          if(prevPrevSec.size()==1){
          std::cout<<"problem two length 1 sections in a row\n";
          }
          newSection = blendLoopToHelix(prevPrevSec[prevPrevSec.size()-2],prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-1],alphaKap,alphaTau,alphaAl,alphaCl,molDat[index].second,generator,suceeded,index);
        }else{
          newSection = blendLoopToHelix(prevSec[prevSec.size()-3],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],alphaKap,alphaTau,alphaAl,alphaCl,molDat[index].second,generator,suceeded,index);
        }
      }else if(prevType=="loop"&&type=="Strand"){
        std::vector<point> prevSec = molPos[index-1];
        if(prevSec.size()==2){
          std::vector<point> prevPrevSec = molPos[index-2];
	  int loopOrStrand=1;
          newSection=blendLoopWithStrand(prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],molDat[index].second,generator,loopOrStrand,suceeded);
        }else if(prevSec.size()==1){
        std::vector<point> prevPrevSec = molPos[index-2];
	      int loopOrStrand=1;
        if(prevPrevSec.size()==1){
          std::cout<<"problem two length 1 sections in a row\n";
        }
        newSection=blendLoopWithStrand(prevPrevSec[prevPrevSec.size()-2],prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-1],molDat[index].second,generator,loopOrStrand,suceeded);
        }
        else{
	  int loopOrStrand=1;
          newSection=blendLoopWithStrand(prevSec[prevSec.size()-3],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],molDat[index].second,generator,loopOrStrand,suceeded);
        }
      }
      else if(prevType=="Strand"&&type=="loop"){
        std::vector<point> prevSec = molPos[index-1];
        if(prevSec.size()==2){
          std::vector<point> prevPrevSec = molPos[index-2];
	  int loopOrStrand=0;
          newSection=blendLoopWithStrand(prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],molDat[index].second,generator,loopOrStrand,suceeded);
        }else if(prevSec.size()==1){
	  std::vector<point> prevPrevSec = molPos[index-2];
	  if(prevPrevSec.size()==1){
	    std::cout<<"problem two length 1 sections in a row\n";
	  }
	  int loopOrStrand=0;
	  newSection=blendLoopWithStrand(prevPrevSec[prevPrevSec.size()-2],prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-1],molDat[index].second,generator,loopOrStrand,suceeded);
        }
        else{
	  int loopOrStrand=0;
          newSection=blendLoopWithStrand(prevSec[prevSec.size()-3],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],molDat[index].second,generator,loopOrStrand,suceeded);
        }
      }
      else{
	// Here we have a helix to loop transition
        std::vector<point> prevSec = molPos[index-1];
	if(type=="Strand"){
	  int loopOrStrand=1;
	  newSection = blendHelixToLoop(prevSec[prevSec.size()-3],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],molDat[index].second,generator,loopOrStrand,suceeded);
	}else{
	  int loopOrStrand=0;
	  newSection = blendHelixToLoop(prevSec[prevSec.size()-3],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],molDat[index].second,generator,loopOrStrand,suceeded);
	}
      }
    }
  }else{
    // here we are at the first index
      point origTan(0.0,0.0,1.0);point origNorm(0.0,1.0,0.0);point origBinorm = origTan.cross(origNorm);
      polyHelix ph;
      std::vector<point> pts;
      if(molDat[0].first=="Helix"){
	for(int i=0;i<molDat[0].second;i++){
	   double len = alphaAl*i;
	   ph.updateYvecGeneral(alphaKap,alphaTau,origTan,origNorm,origBinorm,sp,len);
	   newSection.push_back(ph.getCoord());
	}
      }else if(molDat[0].first=="Strand"){
	int loopOrStrand = 1; 
	newSection = makeRandomSectionWithDist(sp,origTan,origNorm,molDat[0].second,generator,loopOrStrand,suceeded);
      }else{
	int loopOrStrand = 0;
	newSection = makeRandomSectionWithDist(sp,origTan,origNorm,molDat[0].second,generator,loopOrStrand,suceeded);
      }
  }
 // now we need to find the new positions of the rest of the molecule after index
  //   if the current section is a loop we will need to blend into the next type
  maxDist = getMaxDistChange(molPos[index],newSection);
  point lastmol = molPos[index][molPos[index].size()-1];
  if(suceeded==true){
    molPos[index] = newSection;
  }
  if(suceeded==true&&index<(molDat.size()-1)){
    int newIndex;
    std::vector<point> newNewSection;
    for(int l=1;l<=noHelices;l++){
      newIndex = index+l;
      if(l==1){
	// attached to a loop on th left
	std::vector<point> prevSec=molPos[newIndex-1];
	if(prevSec.size()==1){
          std::vector<point> prevPrevSec = molPos[newIndex-2];
          if(prevPrevSec.size()==1){
          std::cout<<"problem two length 1 sections in a row\n";
          }
          newNewSection = blendLoopToHelix(prevPrevSec[prevPrevSec.size()-2],prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-1],alphaKap,alphaTau,alphaAl,alphaCl,molDat[newIndex].second,generator,suceeded,index);
        }else{
          newNewSection = blendLoopToHelix(prevSec[prevSec.size()-3],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],alphaKap,alphaTau,alphaAl,alphaCl,molDat[newIndex].second,generator,suceeded,index);
        }
	double tempMaxDist = getMaxDistChange(molPos[newIndex],newNewSection);
	if(tempMaxDist>maxDist){
	  maxDist=tempMaxDist;
	}
	lastmol = molPos[newIndex][molPos[newIndex].size()-1];
	if(suceeded==true){
	  molPos[newIndex] = newNewSection;
	}
      }else{
	std::vector<point> prevSec = molPos[newIndex-1];
	newNewSection = blendHelixToHelix(prevSec,alphaKap,alphaTau,alphaAl,alphaCl,molDat[newIndex].second,generator,suceeded,index);
	double tempMaxDist = getMaxDistChange(molPos[newIndex],newNewSection);
	if(tempMaxDist>maxDist){
	  maxDist=tempMaxDist;
	}
	lastmol = molPos[newIndex][molPos[newIndex].size()-1];
	if(suceeded==true){
	  molPos[newIndex] = newNewSection;
	}
      }
    }
    if(suceeded==true){
      // now complete the molecule, if it is a helix update the frame state
      for(int j=newIndex+1;j<molDat.size();j++){
	if(j==newIndex+1){
	  point lp = newNewSection[newNewSection.size()-1] + molPos[j][0]-lastmol;
	  lastmol = molPos[j][molPos[j].size()-1];
	  std::vector<point> prevSec = molPos[j];
	  tanToSec(lp,molPos[j]);
	  double tempMaxDist = getMaxDistChange(prevSec,molPos[j]);
	  if(tempMaxDist>maxDist){
	    maxDist=tempMaxDist;
	  }
	}else{
	  point lp = molPos[j-1][molPos[j-1].size()-1]+molPos[j][0]-lastmol;
	  lastmol = molPos[j][molPos[j].size()-1];
	  std::vector<point> prevSec = molPos[j];
	  tanToSec(lp,molPos[j]);
	  double tempMaxDist = getMaxDistChange(prevSec,molPos[j]);
	  if(tempMaxDist>maxDist){
	    maxDist=tempMaxDist;
	  }
	}
      }
    }
  }
  return maxDist;
 }


 double randomMol::reshapeMolSmallVariation(std::vector<std::pair<std::string,int> > &molDat,std::vector<std::vector<point> > &molPos,int &index,point &sp,double variationSize,bool &suceeded){
  std::random_device rdev{};
  std::default_random_engine generator{rdev()};
  std::vector<point> newSection;
  double maxDist;
  std::string type = molDat[index].first;
  if(index>0){
    std::string prevType = molDat[index-1].first;
    //std::cout<<"INDEX NO "<<index<<" "<<prevType<<" "<<type<<" "<<index<<"\n";
    if(type == prevType){
      if(type=="Helix"){
        std::vector<point> prevSec = molPos[index-1];
         newSection = blendHelixToHelix(prevSec,alphaKap,alphaTau,alphaAl,alphaCl,molDat[index].second,generator,suceeded,index);
      }else if(type == "strand"){
        std::vector<point> prevSec = molPos[index-1];
        newSection = blendHelixToHelix(prevSec,betaKap,betaTau,betaAl,betaCl,molDat[index].second,generator,suceeded,index);
      }
    }else{
      if(prevType=="loop"&&type=="Helix"){
        std::vector<point> prevSec = molPos[index-1];
        if(prevSec.size()==2){
          std::vector<point> prevPrevSec = molPos[index-2];
	  newSection = blendLoopToHelixAlter(prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],molPos[index][0],alphaKap,alphaTau,alphaAl,alphaCl,molDat[index].second,generator,variationSize,suceeded,index);
        }else if(prevSec.size()==1){
          std::vector<point> prevPrevSec = molPos[index-2];
          if(prevPrevSec.size()==1){
          std::cout<<"problem two length 1 sections in a row\n";
          }
          newSection = blendLoopToHelixAlter(prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],molPos[index][0],alphaKap,alphaTau,alphaAl,alphaCl,molDat[index].second,generator,variationSize,suceeded,index);
        }else{
	  newSection = blendLoopToHelixAlter(prevSec[prevSec.size()-3],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],molPos[index][0],alphaKap,alphaTau,alphaAl,alphaCl,molDat[index].second,generator,variationSize,suceeded,index);
        }
      }else if(prevType=="loop"&&type=="Strand"){
        //std::vector<point> prevSec = molPos[index-1];
	      int loopOrStrand=1;
        newSection= alterSectionWithDist(molPos[index],loopOrStrand,variationSize,generator,suceeded);
      }
      else if(prevType=="Strand"&&type=="loop"){
	      int loopOrStrand=0;
       newSection= alterSectionWithDist(molPos[index],loopOrStrand,variationSize,generator,suceeded);
      }
      else{
	// Here we have a helix to loop transition
        std::vector<point> prevSec = molPos[index-1];
	if(type=="Strand"){
	  int loopOrStrand=1;
      newSection= alterSectionWithDist(molPos[index],loopOrStrand,variationSize,generator,suceeded);
	}else{
	  int loopOrStrand=0;
	  newSection= alterSectionWithDist(molPos[index],loopOrStrand,variationSize,generator,suceeded);
	}
      }
    }
  }else{
    // here we are at the first index
      point origTan(0.0,0.0,1.0);point origNorm(0.0,1.0,0.0);point origBinorm = origTan.cross(origNorm);
      polyHelix ph;
      std::vector<point> pts;
      if(molDat[0].first=="Helix"){
	for(int i=0;i<molDat[0].second;i++){
	   double len = alphaAl*i;
	   ph.updateYvecGeneral(alphaKap,alphaTau,origTan,origNorm,origBinorm,sp,len);
	   newSection.push_back(ph.getCoord());
	}
      }else if(molDat[0].first=="Strand"){
	int loopOrStrand = 1; 
	newSection= alterSectionWithDist(molPos[index],loopOrStrand,variationSize,generator,suceeded);
      }else{
	int loopOrStrand = 0;
	newSection= alterSectionWithDist(molPos[index],loopOrStrand,variationSize,generator,suceeded);
      }
  }
  maxDist = getMaxDistChange(molPos[index],newSection);
  point lastmol = molPos[index][molPos[index].size()-1];
  if(suceeded==true){
    molPos[index] = newSection;
  }
  // now we need to find the new positions of the rest of the molecule after index
  //   if the current section is a loop we will need to blend into the next type
  //std::cout<<"INDEX NO "<<index<<" "<<type<<" "<<molDat[index+1].first<<"\n";
   if(suceeded==true&&index<(molDat.size()-1)){
     int newIndex = index+1;
    if(type=="loop" || type=="Strand"){
      std::string nextType = molDat[index+1].first;
      std::vector<point> newNewSection;
      if(nextType=="Helix"){
	std::vector<point> prevSec = newSection;
	if(prevSec.size()==2){
	  std::vector<point> prevPrevSec = molPos[index-1];
	  newNewSection =blendLoopToHelixAlter(prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],molPos[index+1][0],alphaKap,alphaTau,alphaAl,alphaCl,molDat[index+1].second,generator,variationSize,suceeded,newIndex);
	}else if(prevSec.size()==1){
	  std::vector<point> prevPrevSec = molPos[index-1];
	  if(prevPrevSec.size()==1){
	    std::cout<<"problem two length 1 sections in a row\n";
	  }
	  newNewSection = blendLoopToHelixAlter(prevPrevSec[prevPrevSec.size()-1],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],molPos[index+1][0],alphaKap,alphaTau,alphaAl,alphaCl,molDat[index+1].second,generator,variationSize,suceeded,newIndex);
	}else{
	  newNewSection = blendLoopToHelixAlter(prevSec[prevSec.size()-3],prevSec[prevSec.size()-2],prevSec[prevSec.size()-1],molPos[index+1][0],alphaKap,alphaTau,alphaAl,alphaCl,molDat[index+1].second,generator,variationSize,suceeded,newIndex);
	}
      }else if(nextType=="Strand"){
	std::vector<point> prevSec = newSection;
	if(prevSec.size()==2){
	  std::vector<point> prevPrevSec = molPos[index-1];
	  int loopOrStrand =1;
	  newNewSection= alterSectionWithDist(molPos[index+1],loopOrStrand,variationSize,generator,suceeded);
	}else if(prevSec.size()==1){
	  std::vector<point> prevPrevSec = molPos[index-1];
	  if(prevPrevSec.size()==1){
	    std::cout<<"problem two length 1 sections in a row\n";
	  }
	  int loopOrStrand =1;
	   newNewSection= alterSectionWithDist(molPos[index+1],loopOrStrand,variationSize,generator,suceeded);
	}
	else{
	  int loopOrStrand =1;
	  newNewSection= alterSectionWithDist(molPos[index+1],loopOrStrand,variationSize,generator,suceeded);
	}
      }else{
	// here it is a strand to loop transistion
	std::vector<point> prevSec = newSection;
	if(prevSec.size()==2){
	  std::vector<point> prevPrevSec = molPos[index-1];
	  int loopOrStrand =0;
	  newNewSection= alterSectionWithDist(molPos[index+1],loopOrStrand,variationSize,generator,suceeded);
	}else if(prevSec.size()==1){
	  std::vector<point> prevPrevSec = molPos[index-1];
	  if(prevPrevSec.size()==1){
	    std::cout<<"problem two length 1 sections in a row\n";
	  }
	  int loopOrStrand =0;
	   newNewSection= alterSectionWithDist(molPos[index+1],loopOrStrand,variationSize,generator,suceeded);
	}
	else{
	  int loopOrStrand =0;
	  newNewSection= alterSectionWithDist(molPos[index+1],loopOrStrand,variationSize,generator,suceeded);
	}
      }
      double tempMaxDist = getMaxDistChange(molPos[index+1],newNewSection);
      if(tempMaxDist>maxDist){
	maxDist=tempMaxDist;
      }
      lastmol = molPos[index+1][molPos[index+1].size()-1];
      if(suceeded==true){
	molPos[index+1] = newNewSection;
	// now complete the molecule
	for(int j=index+2;j<molDat.size();j++){
	  if(j==index+2){
	    point lp = newNewSection[newNewSection.size()-1] + molPos[j][0]-lastmol;
	    lastmol = molPos[j][molPos[j].size()-1];
	    std::vector<point> prevSec = molPos[j];
	    tanToSec(lp,molPos[j]);
	    double tempMaxDist = getMaxDistChange(prevSec,molPos[j]);
	    if(tempMaxDist>maxDist){
	      maxDist=tempMaxDist;
	    }
	  }else{
	    point lp = molPos[j-1][molPos[j-1].size()-1]+molPos[j][0]-lastmol;
	    lastmol = molPos[j][molPos[j].size()-1];
	    std::vector<point> prevSec = molPos[j];
	    tanToSec(lp,molPos[j]);
	    double tempMaxDist = getMaxDistChange(prevSec,molPos[j]);
	    if(tempMaxDist>maxDist){
	      maxDist=tempMaxDist;
	    }
	  }
	}
      }
    }else{
      for(int j=index+1;j<molDat.size();j++){
	if(j==index+1){
	  point lp = newSection[newSection.size()-1] + molPos[j][0] -lastmol;
	  lastmol = molPos[j][molPos[j].size()-1];
	  std::vector<point> prevSec = molPos[j];
	  tanToSec(lp,molPos[j]);
	  double tempMaxDist = getMaxDistChange(prevSec,molPos[j]);
	  if(tempMaxDist>maxDist){
	    maxDist=tempMaxDist;
	  }
	}else{
	  point lp = molPos[j-1][molPos[j-1].size()-1]+ molPos[j][0] -lastmol;
	  lastmol = molPos[j][molPos[j].size()-1];
	  std::vector<point> prevSec = molPos[j];
	  tanToSec(lp,molPos[j]);
	  double tempMaxDist = getMaxDistChange(prevSec,molPos[j]);
	  if(tempMaxDist>maxDist){
	    maxDist=tempMaxDist;
	  }
	}
      }
    }
  }
  //now update all the relevant Tangents
  return maxDist;
 }

