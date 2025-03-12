// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export(name="d2_logm_eta_row")]]
double  d2_logm_eta_row(const arma::rowvec& eta, const arma::rowvec& y, arma::mat& res, arma::vec& z, arma::vec& w) {
  using namespace arma;

  uint32_t d = y.n_elem;
  uint32_t q = eta.n_elem;

  uint32_t j;
  uint32_t k;
  uint32_t kp;
  uint32_t l;
  uint32_t m;
  uint32_t c1;
  uint32_t zl;// auxiliary access to z elements
  uint32_t wl;// auxiliary access to w elements
  uint32_t zm;// auxiliary access to z elements
  uint32_t wm;// auxiliary access to w elements

  rowvec r(d, fill::zeros); //i-th residual vector
  mat Theta(d, d, fill::zeros); //i-th log Sigma matrix
  vec L(d, fill::zeros); //eigenvalues - Lambda in the paper
  mat U(d, d, fill::zeros); //eigenvectors

  vec s(d, fill::zeros);
  mat out(q, q, fill::zeros);
  mat loewner(d, d, fill::zeros);
  mat F(d, d, fill::zeros);
  mat Pi(d, d, fill::zeros);
  mat tildeU(d, d * (d - 1)/2, fill::zeros);
  mat tilde_loew(d, d, fill::zeros);
  mat A(d, d, fill::zeros);
  cube star_loew(d, d, d, fill::zeros);
  cube A1(d, d, d, fill::zeros);
  cube A2(d, d, d, fill::zeros);

    // i-th residuals vector
    r = y - eta(span(0, d - 1));

    // i-th logarithm of the covariance matrix
    c1 = 0;
    for(j = 1; j < d; j++){
      for(k = 0; k < j; k++){
        Theta.at(j, k) = eta[c1 + 2 * d];
        Theta.at(k, j) = Theta.at(j, k);
        c1 += 1;
      }
    }
    Theta.diag() = eta(span(d, 2 * d - 1));

    // eigen decomposition
    eig_sym(L, U, Theta);

    //loewner matrix
    for(j = 0; j < (d - 1); j++){
      for(k = j + 1; k < d; k++){
        loewner.at(j, k)= ((exp(-L[j])) - exp(-L[k]))/(L[k] - L[j]);
        loewner.at(k, j) = loewner.at(j, k);
      }
    }
    loewner.diag() =  exp(-L);


    for(j = 0; j < d; j++){
      for(k = 0; k < d; k++){
        if(k != j) {
          tilde_loew.at(j, k) = (loewner.at(j, k) - loewner.at(j, j))/(L[j] - L[k]);
          for(kp = 0; kp < d; kp++){
            if(kp != k and kp!= j) star_loew.at(j, k, kp) = (loewner.at(kp, k) - loewner.at(j, k))/(L[j] - L[kp]);
          }
        }
      }
    }

    s = U.t() * r.t();
    F = U.each_row() % s.t();
    Pi = F * loewner;
    A = F * tilde_loew.t();

    for(j = 0; j < d; j++){
      A1.slice(j) = F * (F.each_row() % tilde_loew.col(j).t()).t();
      A2.slice(j) = F  * (F * star_loew.slice(j)).t(); //star_loew.slice(j).t() is unuseful since star_loew is symmetric
    }


    for(l = 0; l < d * (d - 1)/2; l++){
      zl = z[l];
      wl = w[l];
      for(j = 0; j < d; j++){
        tildeU.at(j, l) = Pi.at(zl, j) * U.at(wl, j) + Pi.at(wl, j) * U.at(zl, j);
      }
    }

    //Block (1,1) + Block(1,2)
    for(l = 0; l < d; l++){
      for(m = l; m < d; m++){
        out(l, m) = -as_scalar(U.row(l) * ((U.row(m)).t() % loewner.diag()));
        for(j = 0; j < d; j++){
          out.at(l + d, m + d) -= U.at(l, j) * U.at(m, j) * (F.at(l, j) * F.at(m, j) * loewner.at(j, j)/2 + F.at(m, j) * A.at(l, j) + F.at(l, j) * A.at(m, j)+ A1.at(l, m, j) + A2.at(l, m, j));
        }
      }
      for(m = d; m < 2 * d; m++){
        out.at(l, m) = -as_scalar(U.row(l) * (U.row(m - d) % Pi.row(m - d)).t());
      }
      for(m = 2 * d; m < q; m++){
        zm = z[m - 2 * d];
        wm = w[m - 2 * d];
        out.at(l, m) = -as_scalar(U.row(l) * tildeU.col(m - 2 * d));
        for(j = 0; j < d; j++){
          out.at(l + d, m) -= U.at(l, j) * (U.at(wm, j) * (F.at(zm, j) * F.at(l, j) * loewner.at(j, j)/2 + F.at(zm, j) * A.at(l, j) + A.at(zm, j) * F.at(l, j) + A1.at(l, zm, j) + A2.at(l, zm, j)) +
                                            U.at(zm, j) * (F.at(wm, j) * F.at(l, j) * loewner.at(j, j)/2 + F.at(wm, j) * A.at(l, j) + A.at(wm, j) * F.at(l, j) + A1.at(l, wm, j) + A2.at(l, wm, j)));
        }
      }
    }

    //  Block(2,2)
    for(l = 2 * d; l < q; l++){
      zl = z[l - 2 * d];
      wl = w[l - 2 * d];
      for(m = l; m < q; m++){
        zm = z[m - 2 * d];
        wm = w[m - 2 * d];
        for(j = 0; j < d; j++){
          out.at(l, m) -= (U.at(wl, j) * (U.at(wm, j) * (F.at(zm, j) * F.at(zl, j) * loewner.at(j, j)/2 + F.at(zm, j) * A.at(zl, j) + A.at(zm, j) * F.at(zl, j) + A1.at(zl, zm, j) + A2.at(zm, zl, j)) +
                                          U.at(zm, j) * (F.at(wm, j) * F.at(zl, j) * loewner.at(j, j)/2 + F.at(wm, j) * A.at(zl, j) + A.at(wm, j) * F.at(zl, j) + A1.at(zl, wm, j) + A2.at(wm, zl, j)))+
                           U.at(zl, j) * (U.at(wm, j) * (F.at(zm, j) * F.at(wl, j) * loewner.at(j, j)/2 + F.at(zm, j) * A.at(wl, j) + A.at(zm, j) * F.at(wl, j) + A1.at(wl, zm, j) + A2.at(zm, wl, j))+
                                          U.at(zm, j) * (F.at(wm, j) * F.at(wl, j) * loewner.at(j, j)/2 + F.at(wm, j) * A.at(wl, j) + A.at(wm, j) * F.at(wl, j) + A1.at(wl, wm, j) + A2.at(wm, wl, j))));
        }
      }
    }

  return(1.0);
}
