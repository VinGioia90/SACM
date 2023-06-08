// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export(name="d2_mcd_eta_row")]]
double d2_mcd_eta_row(const arma::rowvec& eta, const arma::rowvec& y, arma::mat& d2l,   arma::vec& z, arma::vec& w, arma::mat& G, arma::vec& t, Rcpp::List&  idx_jk){
  using namespace arma;
  //uint32_t n = y.n_rows;
  uint32_t d = y.n_elem;
  uint32_t q = eta.n_elem;

  double ee_s;
  double ee;
  double eee;
  double sj;
  uint32_t wj;
  uint32_t k1;
  double swj;
  uint32_t zj;
  double ee_r;
  Rcpp::IntegerVector idx_k;

  double aux1 = 0.0;

  mat out(q, q, fill::zeros);
  vec s(d, fill::zeros);
  rowvec r(d, fill::zeros);
  Rcpp::IntegerVector ik1;

  uint32_t j;
  uint32_t k;
  uint32_t l;

    r = y - eta(span(0, d - 1));

    // sum in brackets (computed one time and used below)
    s[0] = r[0];
    for(j = 1; j < d; j++){
      for(k = 0; k < j; k++){
        aux1 += r[k] * eta[G.at(j - 1, k)];
      }
      s[j] = aux1 + r[j];
      aux1 = 0.0;
    }

    for(j = 0; j < d; j++){
      sj = s[j];
      ee = -exp(-eta[j + d]);
      ee_s =  ee * sj;

      // Block (1,1) - mean vs mean
      out.at(j, j) = ee;

      for(k = j + 1; k < d; k++){
        eee = exp(-eta[k + d]) * eta[G.at(k - 1, j)];
        // Block (1,1) - mean vs mean
        out.at(j, k) = -eee;
        for(l = k + 1; l < d; l++){
          aux1 +=  -exp(-eta[l + d]) * eta[G.at(l - 1, j)] * eta[G.at(l - 1,k)];
        }

        // Block (1,1) - mean vs mean
        out.at(j, k) += aux1;
        aux1 = 0.0;

        // Block (1,2) - mean vs logD2
        out.at(j, k + d) = -eee * s[k] ;
        //Block (1,1) - mean vs mean
        out.at(j, j) -=  eee * eta[G.at(k - 1, j)];
      }


      // Block (1,2) - mean vs logD2
      out.at(j, j + d) = ee_s;

      // Block (2,2) - logD2 vs logD2
      out.at(j + d, j + d) = 0.5 * ee_s * sj;
    }

    for(k = 0; k < d * (d - 1)/2; k++){
      k1 = k + 2 * d;
      wj = w[k];
      swj = s[wj];
      zj = z[k];
      ee = exp(-eta[wj + d]);
      ee_r = ee * r[zj];
      //Block (1,3) - mean vs T
      for(j = 0; j < wj; j++){
        out.at(j, k1) = ee_r * eta[G.at(wj - 1, j)];
      }
      out.at(zj, k1) += ee * swj;
      out.at(wj, k1) = ee_r;

      //Block (2,3) - logD2 vs T
      out.at(wj + d, k1) = ee_r * swj;

      //Block (3,3) - T vs T
      for(j = 0; j < t[k]; j++){
        out.at(k1, k1 + j) = -ee_r * r(z[k+j]);
      }
    }

  return(1.0);
}
