#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  //dati
  DATA_MATRIX(Y);  // Response matrix dim= n x d 
  
  //parameters
  PARAMETER_MATRIX(eta); //size n x (d + d x (d+1)/2) 
  
  Type out = 0.0;     // Negative log likelihood function
  Type aux = 0.0;
  Type eij = 0.0;

  int n = Y.rows();
  int d = Y.cols();
  int i;
  int j;
  int k;
  int count;
  vector<Type> r(d);
  for(i = 0; i < n; i++){
    for(j = 0; j < d; j++){
     r(j) = Y(i,j)-eta(i,j); 
    };
    eij = eta(i,d); 
    out += 0.5*eta(i,d) + 0.5*exp(-eta(i,d))*r(0)*r(0);
         
    count = 2*d;
    for(j = d+1; j < 2*d; j++){
     aux = 0.0; 
     eij = eta(i,j); 
     for(k = d; k < j; k++){
       aux +=  r(k-d)*eta(i,count);
       count += 1; 
     };
     out += 0.5*eij + 0.5*exp(-eij)*(aux+r(j-d))*(aux+r(j-d));
     };
  };
  
   return -out;
};
