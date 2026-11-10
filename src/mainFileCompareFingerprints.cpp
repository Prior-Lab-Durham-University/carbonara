#include "writheFP.h"
#include <string.h>
#include <string>
#include <iostream>



int main( int argc, const char* argv[] )
{
       if(argc >=2){
	bool check=false;
	std::ifstream myfile;
 	myfile.open(argv[1]);
	std::ifstream myfile2;
 	myfile2.open(argv[2]);
	std::string output;
	std::vector<point> points1;
	if (myfile.is_open()) { 
	  while(std::getline(myfile,output)){
	    points1.push_back(point(output));
	  }
	}else{
	std::cout<<"Curve data file 1 failed to open";
	}
	std::vector<point> points2;
	if (myfile2.is_open()) { 
	  while(std::getline(myfile2,output)){
	    points2.push_back(point(output));
	  }
	}else{
	std::cout<<"Curve data file 2 failed to open";
	}
	if((points1.size()>3) && (points2.size()>3)){
	    check =true;
	}
	myfile.close();myfile2.close();
	if(check){
	  writheFP w;
	  int len1 = points1.size();int len2 = points2.size();
	  double cutOff = std::atof(argv[3]);
	  std::pair<std::vector<std::vector<point> >,std::vector<std::vector<point> > >  fp1=w.getWritheFingerprints(points1);
	  std::pair<std::vector<std::vector<point> >,std::vector<std::vector<point> > >  fp2=w.getWritheFingerprints(points2);
	  std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > > fSetBest;
	  int len1Best=0;
	  for(int i=0;i<10;i++){
	    double cutofftest = cutOff - i*(cutOff-0.01)/9;
	    //std::cout<<cutofftest<<"\n";
	    std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > > fSet =w.compareFingerPrints(fp1.first,len1,fp2.first,len2,cutofftest);
	    int lenTot1=0;int lenTot2=0;
	    for(int i=0;i<fSet.size();i++){
	      if(lenTot1+fSet[i].first.second-fSet[i].first.first>0){
		lenTot1 = lenTot1+fSet[i].first.second-fSet[i].first.first+1;
	      }
	      if(lenTot1+fSet[i].second.second-fSet[i].second.first>0){
		lenTot2 = lenTot2+fSet[i].second.second-fSet[i].second.first+1;
	      }
	    }
	    //std::cout<<lenTot1<<"\n";
	    if(lenTot1>len1Best || i==0){
	      fSetBest = fSet;
	      len1Best = lenTot1;
	    }
	  }
	   int lenTot1=0;int lenTot2=0;
	  for(int i=0;i<fSetBest.size();i++){
	      std::cout<<fSetBest[i].first.first<<" "<<fSetBest[i].first.second<<" "<<fSetBest[i].second.first<<" "<<fSetBest[i].second.second<<" ";
	    if(lenTot1+fSetBest[i].first.second-fSetBest[i].first.first>0){
	      lenTot1 = lenTot1+fSetBest[i].first.second-fSetBest[i].first.first+1;
	    }
	    if(lenTot1+fSetBest[i].second.second-fSetBest[i].second.first>0){
	      lenTot2 = lenTot2+fSetBest[i].second.second-fSetBest[i].second.first+1;
	    }
	  }
	  std::cout<<double(lenTot1)/double(len1)<<" ";
	  std::cout<<double(lenTot2)/double(len2)<<" ";
	  std::cout<<"\n";
	  
	}else{
	  std::cout<<"need more than 3 points in the curve\n";
	}
       }else{
	  std::cout<<"must supply a file containing the curve data";
       }
       return 0;
}
