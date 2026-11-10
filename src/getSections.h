#ifndef GET_SECTIONS_H
#define GET_SECTIONS_H

#include "point.h"

class getSections{
      public: 
    	 getSections(std::vector<point>& lst);
    	 void checkTurn(point& n1,int n1pos,point& n2,int n2pos,point& n3,int n3pos);
         std::vector<int> getListTurnPts();
      private: 
         std::vector<point> list;
         double dif1,dif2,prod;
         std::vector<int> listSections /*think this needs to be intialised as empty*/;
      };
      
#endif
