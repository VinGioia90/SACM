#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  //dati
  DATA_MATRIX(Y);  // Response matrix dim= n x d

  //parameters
  PARAMETER_MATRIX(eta); //size n x (d + d x (d+1)/2)

  Type ll = 0.0;     // Negative log likelihood function

  int n = Y.rows();
  int d = Y.cols();
  int i;
  int j;
  int k;
  int count = 0;

  matrix<Type> nTheta(d,d); //negative log covariance
  matrix<Type> r(d,1);
  Type rexpThetar;

  for(i=0; i<n; i++){


   for(j = 0; j < d; j++){
    r(j,0) = Y(i,j)-eta(i,j);
    nTheta(j,j) = -eta(i,j+d);
   }

   count = 0;
   for(j = 1; j < d; j++){
      for(k = 0; k < j; k++){
        nTheta(j, k) = -eta(i, count + 2 * d);
        nTheta(k, j) = nTheta(j, k);
        count += 1;
      };
    };


   rexpThetar = (r.transpose()*atomic::expm(nTheta)*r)(0);
   ll += +0.5*nTheta.trace()  -0.5* rexpThetar ;
 };

   return ll;
};
