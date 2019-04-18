#ifndef DISCONTINUOUS_H
#define DISCONTINUOUS_H

#include <cassert>
#include <algorithm>

class DISCONTINUOUS{
 protected:
  vector<double> x, data;

 public:
  void resize(size_t nx){
    x.resize(nx);
    data.resize(nx);
  }

  size_t index(double xin) const{
    assert(x.size()>0);
    assert(xin >= x[0]);
    assert(xin <= x[x.size()-1]);
    size_t index = std::upper_bound(x.begin(), x.end(), xin) - x.begin()-1;
    assert(index<x.size()-1);
    return index;
  }
  
  double Derivative(double xin, int ind=-1) const{
    if(ind<0) ind = index(xin);
    double slope = (data[ind+1]-data[ind]) / (x[ind+1] - x[ind]);
    assert(slope==slope);
    return slope;
  }
  
  double operator()(double xin) const{
    int ind = index(xin);
    double slope = Derivative(xin, ind);
    double result = data[ind] + slope * (xin - x[ind]);
    assert(result==result);
    return result;
  }

  void SetData(const vector<double>& x, const vector<double>& data){
    this->x = x;
    this->data = data;
  }
  
  void Open(string filename){
    ifstream fin(filename.c_str());
    double xin, d; 
    while((fin>>xin) && (fin>>d)){
      x.push_back(xin);
      data.push_back(d);
    }
    fin.close();
    assert(x.size()>0);
    assert(data.size()==x.size());
  }

 void Open(string filename,char ignore){
   ifstream fin(filename.c_str());
   double xin, d; 
   string line;
   while(fin.peek()==ignore)
     getline(fin,line);
   while((fin>>xin) && (fin>>d)){
     x.push_back(xin);
     data.push_back(d);
   }
   fin.close();
   assert(x.size()>0);
   assert(data.size()==x.size());
 }
 
  double XMin() const{
    return x[0];
  }
  double XMax() const{
    return x[x.size()-1];
  }

  DISCONTINUOUS copy_logy() const{
    assert(this->x.size()>0);
    assert(this->x.size()==this->data.size());
    
    DISCONTINUOUS output;
    output.x = this->x;
    output.data.resize(this->x.size());
    for(size_t i=0; i<output.x.size(); i++)
      output.data[i] = log(this->data[i]);
    
    assert(output.x.size()>0);
    assert(output.x.size()==output.data.size());

    return output;
  }
};


#endif
