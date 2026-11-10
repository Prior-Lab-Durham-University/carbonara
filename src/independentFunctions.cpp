#include <sstream>
#include <iostream>
#include <vector>

static const char* doubleToChar(double &din){
  std::stringstream ss;
  ss<<din;
  const char* str = ss.str().c_str();
  return str;      	       
}

static const char* intToChar(int &iin){
  std::stringstream ss;
  ss<<iin;
  const char* str = ss.str().c_str();
  return str;      	       
}


template< class T >
void reorder(std::vector<T> &v, std::vector<size_t> const &order ){   
  std::vector<T> vnew;
  for(int l=0;l<v.size();l++){
    vnew.push_back(v[order[l]]);
    }
  v = vnew;
}
