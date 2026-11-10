#include "getSections.h"

getSections::getSections(std::vector<point>& lst){
      list =lst;                                               
};

void getSections::checkTurn(point& n1,int n1pos,point& n2,int n2pos,point& n3,int n3pos){
      double prod = n1.getZ()*n2.getZ();
      if(prod<=0){
         listSections.push_back(n1pos);
         listSections.push_back(n2pos);
      }
}

std::vector<int> getSections::getListTurnPts(){
            listSections.push_back(0);
            int npts = list.size();
			for(int i=0;i<list.size()-2;i++){
				checkTurn(list[i],i,list[i+1],i+1,list[i+2],i+2);  
            }
            listSections.push_back(npts-1);
            return listSections;
}
