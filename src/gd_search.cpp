#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double obj_fun(const arma::mat &A, const arma::vec &U2) {
  return arma::as_scalar(U2.t() * A * U2);
}

// [[Rcpp::export]]
double c_obj_fun(const arma::mat &M, const arma::vec &C) {
  arma::mat M_inv = arma::inv_sympd(M);
  return arma::as_scalar(C.t() * M_inv * C);
}

// [[Rcpp::export]]
arma::mat gen_m(const arma::mat &X, const arma::mat &A) {
  return X.t() * A * X;
}

// [[Rcpp::export]]
arma::vec gen_u(const arma::mat &M, const arma::mat &X, const arma::vec &C) {
  arma::mat M_inv = arma::inv_sympd(M);
  return C.t() * M_inv * X.t();
}

// [[Rcpp::export]]
double remove_one(const arma::mat &A, arma::uword i, const arma::vec &u) {
  arma::vec idx = arma::linspace(0, A.n_rows - 1, A.n_rows);
  arma::uvec uidx = arma::conv_to<arma::uvec>::from(idx);
  uidx.shed_row(i);

  double d = A(i, i);
  arma::vec b = A.submat(uidx, arma::uvec({i}));

  return obj_fun(A.submat(uidx, uidx) - b * b.t() / d, u(uidx));
}

// [[Rcpp::export]]
arma::mat remove_one_mat(const arma::mat &A, arma::uword i) {
  arma::vec idx = arma::linspace(0, A.n_rows - 1, A.n_rows);
  arma::uvec uidx = arma::conv_to<arma::uvec>::from(idx);
  uidx.shed_row(i);

  double d = A(i, i);
  arma::vec b = A.submat(uidx, arma::uvec({i}));

  return A.submat(uidx, uidx) - b * b.t() / d;
}

// [[Rcpp::export]]
double add_one(const arma::mat &A, double sigma_jj, const arma::vec &f,
               const arma::vec &u) {
  arma::mat A2(A.n_rows + 1, A.n_cols + 1, arma::fill::zeros);
  A2.submat(0, 0, A2.n_rows - 2, A2.n_cols - 2) = A;
  // for(arma::uword j =0; j < A2.n_rows-1; j++){
  //   A2(j,A2.n_cols-1) = 0;
  //   A2(A2.n_rows-1,j) = 0;
  // }
  A2(A2.n_rows - 1, A2.n_cols - 1) = 1 / sigma_jj;

  // step 3: compute K2_inv
  arma::vec u1 = arma::join_cols(f, arma::vec({0}));
  arma::vec v1(u1.n_elem, arma::fill::zeros);
  v1(v1.n_elem - 1) = 1;
  A2 = A2 -
       ((A2 * u1) * (v1.t() * A2)) / (1 + arma::as_scalar((v1.t() * A2) * u1));

  // step 4: compute K3_inv
  A2 = A2 -
       ((A2 * v1) * (u1.t() * A2)) / (1 + arma::as_scalar((u1.t() * A2) * v1));

  return obj_fun(A2, u);
}

// [[Rcpp::export]]
arma::mat add_one_mat(const arma::mat &A, double sigma_jj, const arma::vec &f) {
  arma::mat A2(A.n_rows + 1, A.n_cols + 1, arma::fill::zeros);
  // step 1: compute A*
  A2.submat(0, 0, A2.n_rows - 2, A2.n_cols - 2) = A;
  for (arma::uword j = 0; j < A2.n_rows - 1; j++) {
    A2(j, A2.n_cols - 1) = 0;
    A2(A2.n_rows - 1, j) = 0;
  }
  A2(A2.n_rows - 1, A2.n_cols - 1) = 1 / sigma_jj;

  // step 3: compute K2_inv
  arma::vec u1 = arma::join_cols(f, arma::vec({0}));
  arma::vec v1(u1.n_elem, arma::fill::zeros);
  v1(v1.n_elem - 1) = 1;
  A2 = A2 -
       ((A2 * u1) * (v1.t() * A2)) / (1 + arma::as_scalar((v1.t() * A2) * u1));

  // step 4: compute K3_inv
  A2 = A2 -
       ((A2 * v1) * (u1.t() * A2)) / (1 + arma::as_scalar((u1.t() * A2) * v1));

  return A2;
}

arma::uvec std_setdiff(arma::uvec &x, arma::uvec &y) {
  std::vector<int> a = arma::conv_to<std::vector<int> >::from(arma::sort(x));
  std::vector<int> b = arma::conv_to<std::vector<int> >::from(arma::sort(y));
  std::vector<int> out;

  std::set_difference(a.begin(), a.end(), b.begin(), b.end(),
                      std::inserter(out, out.end()));

  return arma::conv_to<arma::uvec>::from(out);
}

// [[Rcpp::export]]
double ChooseSwap(arma::uvec idx_in, arma::mat A, arma::mat sig, arma::vec u,
                  arma::uvec &out2, arma::mat &out3) {
  // generate the complete index
  arma::vec idx = arma::linspace(0, sig.n_rows - 1, sig.n_rows);
  arma::uvec uidx = arma::conv_to<arma::uvec>::from(idx);

  // get the index not included in complete index
  arma::uvec idx_out = std_setdiff(uidx, idx_in);

  // find one index from idx_in to remove
  // which results in largest val of remove_one()
  arma::vec u_idx_in = u.elem(idx_in);
  arma::vec val_out_vec(idx_in.n_elem, arma::fill::zeros);
  for (std::size_t i = 0; i < idx_in.n_elem; ++i) {
    val_out_vec(i) = remove_one(A, i, u_idx_in);
  }
  int idx_rm = val_out_vec.index_max();

  // compute the new A without the index of idx_rm
  arma::mat rm1A = remove_one_mat(A, idx_rm);

  // remove index idx_rm from idx_in
  idx_in.shed_row(idx_rm);

  // find one index from idx_out to add (swap)
  // which results in largest val of add_one()
  arma::vec val_in_vec(idx_out.n_elem, arma::fill::zeros);
  for (arma::uword i = 0; i < idx_out.n_elem; ++i) {
    arma::uword ii = idx_out(i);
    val_in_vec(i) =
          add_one(rm1A, sig(ii, ii), sig.submat(idx_in, arma::uvec({ii})),
                  u.elem(arma::join_cols(idx_in, arma::uvec({ii}))));
  }
  arma::uword idx_swap = idx_out(val_in_vec.index_max());

  // compute the new A with the index of idx_swap
  arma::mat newA = add_one_mat(rm1A, sig(idx_swap, idx_swap),
                               sig.submat(idx_in, arma::uvec({idx_swap})));
  idx_in = arma::join_cols(idx_in, arma::uvec({idx_swap}));

  out2 = idx_in;
  out3 = newA;
  return obj_fun(newA, u.elem(idx_in));
}

// [[Rcpp::export]]
arma::uvec Grad(arma::uvec idx_in, arma::mat A, arma::mat sig, arma::vec u,
                double tol = 1e-9, bool trace = true) {
  double new_val = obj_fun(A, u.elem(idx_in));
  double diff = 1.0;
  int i = 0;
  while (diff > tol) {
    double val = new_val;
    i = i + 1;

    arma::uvec out2;
    arma::mat out3;
    new_val = ChooseSwap(idx_in, A, sig, u, out2, out3);
    diff = new_val - val;
    if (diff > 0) {
      A = out3;
      idx_in = out2;
    }
    if (trace) Rcpp::Rcout << "\rIter: " << i << " " << diff << std::endl;
  }
  return (idx_in);
}

// [[Rcpp::export]]
double ChooseSwapRobust(arma::uword nlist, arma::uvec idx_in, const arma::mat &A_list,
                        const arma::mat &sig_list, const arma::vec &u_list, const arma::vec &weights,
                        arma::uvec &out2, arma::mat &out3) {
  const arma::uword u_nrows = u_list.n_rows / nlist;
  const arma::uword A_nrows = A_list.n_rows / nlist;

  // generate the complete index
  arma::vec idx = arma::linspace(0, u_nrows - 1, u_nrows);
  arma::uvec uidx = arma::conv_to<arma::uvec>::from(idx);

  // get the index not included in complete index
  arma::uvec idx_out = std_setdiff(uidx, idx_in);

  // find one index from idx_in to remove
  // which results in largest val of remove_one()
  arma::mat val_out_mat(idx_in.n_elem, nlist, arma::fill::zeros);
#pragma omp parallel for
  for (std::size_t j = 0; j < nlist; ++j) {
    arma::mat A = A_list.rows(j*A_nrows, (j+1)*A_nrows-1);
    arma::vec u = u_list.rows(j*u_nrows, (j+1)*u_nrows-1);
    arma::vec u_idx_in = u.elem(idx_in);
    for (std::size_t i = 0; i < idx_in.n_elem; ++i) {
      val_out_mat(i,j) = remove_one(A, i, u_idx_in);
    }
  }

  arma::vec val_out_vec = val_out_mat * weights;
  int idx_rm = val_out_vec.index_max();

  // compute the new A without the index of idx_rm
  arma::mat rm1A_list;
  for (std::size_t j = 0; j < nlist; ++j) {
    arma::mat A = A_list.rows(j*A_nrows, (j+1)*A_nrows-1);
    rm1A_list = arma::join_cols(rm1A_list, remove_one_mat(A, idx_rm));
  }

  // remove index idx_rm from idx_in
  idx_in.shed_row(idx_rm);

  // find one index from idx_out to add (swap)
  // which results in largest val of add_one()
  arma::mat val_in_mat(idx_out.n_elem,nlist,arma::fill::zeros);
#pragma omp parallel for
  for (arma::uword j = 0; j < nlist; ++j) {
    arma::mat sig = sig_list.rows(j*u_nrows, (j+1)*u_nrows-1);
    arma::vec u = u_list.subvec(j*u_nrows, (j+1)*u_nrows-1);
    arma::mat rm1A = rm1A_list.rows(j*(A_nrows-1), (j+1)*(A_nrows-1)-1);
    for (arma::uword i = 0; i < idx_out.n_elem; ++i) {
      arma::uword ii = idx_out(i);
      val_in_mat(i,j) =
        add_one(rm1A, sig(ii, ii), sig.submat(idx_in, arma::uvec({ii})),
                u.elem(arma::join_cols(idx_in, arma::uvec({ii}))));
    }
  }

  arma::vec val_in_vec = val_in_mat * weights;
  arma::uword indexmax = val_in_vec.index_max();
  arma::uword idx_swap = idx_out(indexmax);

  // compute the new A with the index of idx_swap
  arma::mat newA_list;
  for (std::size_t j = 0; j < nlist; ++j) {
    arma::mat sig = sig_list.rows(j*u_nrows, (j+1)*u_nrows-1);
    arma::mat rm1A = rm1A_list.rows(j*(A_nrows-1), (j+1)*(A_nrows-1)-1);
    newA_list = arma::join_cols(newA_list, add_one_mat(rm1A, sig(idx_swap, idx_swap),
                                                       sig.submat(idx_in, arma::uvec({idx_swap}))));
  }

  idx_in = arma::join_cols(idx_in, arma::uvec({idx_swap}));

  out2 = idx_in;
  out3 = newA_list;

  return val_in_vec(indexmax);
}

// [[Rcpp::export]]
arma::uvec GradRobust(arma::uword nlist, arma::uvec idx_in, arma::mat A_list,
                      arma::mat M_list, arma::vec C_list, arma::mat X_list,
                      arma::mat sig_list, arma::vec u_list, arma::vec weights,
                      double tol = 1e-9, bool trace = true) {
  const arma::uword u_nrows = u_list.n_rows / nlist;
  const arma::uword A_nrows = A_list.n_rows / nlist;
  const arma::uword M_nrows = M_list.n_rows / nlist;
  const arma::uword C_nrows = C_list.n_rows / nlist;
  const arma::uword X_nrows = X_list.n_rows / nlist;

  arma::vec new_val_vec(nlist, arma::fill::zeros);
  for (arma::uword j = 0; j < nlist; ++j) {
    arma::mat M = M_list.rows(j*M_nrows, (j+1)*M_nrows-1);
    arma::vec C = C_list.subvec(j*C_nrows, (j+1)*C_nrows-1);
    new_val_vec(j) = c_obj_fun(M, C);
  }
  double new_val = arma::dot(new_val_vec, weights);

  double diff = -1.0;
  int i = 0;
  // we now need diff to be negative
  while (diff < 0) {
    double val = new_val;
    i = i + 1;

    arma::uvec out2;
    arma::mat out3;

    new_val = ChooseSwapRobust(nlist, idx_in, A_list, sig_list, u_list, weights,
                               out2, out3);

    //we have to now recalculate all the lists of matrices for the new design proposed by the swap
    A_list = out3;
    for (arma::uword j = 0; j < nlist; ++j) {
      arma::mat X = X_list.rows(j*X_nrows, (j+1)*X_nrows-1);
      arma::mat A = A_list.rows(j*A_nrows, (j+1)*A_nrows-1);
      M_list.rows(j*M_nrows, (j+1)*M_nrows-1) = gen_m(X.rows(out2),A);

      //check positive semi-definite before proceeding
      // if(!check_psd(M_list[[j]]))stop(paste0("M not positive semi-definite. Column ",which(colSums(M_list[[j]])==0),
      //               " of design ",j," is not part of an optimal design."))

      arma::mat C = C_list.rows(j*C_nrows, (j+1)*C_nrows-1);
      arma::mat M = M_list.rows(j*M_nrows, (j+1)*M_nrows-1);
      arma::mat cM = C.t() * M.i();
      u_list.rows(j*u_nrows, (j+1)*u_nrows-1) = X * cM.t();
    }

    // calculate values for the new design - this is changed to new objective function
    for (arma::uword j = 0; j < nlist; ++j) {
      arma::mat M = M_list.rows(j*M_nrows, (j+1)*M_nrows-1);
      arma::vec C = C_list.subvec(j*C_nrows, (j+1)*C_nrows-1);
      new_val_vec(j) = c_obj_fun(M, C);
    }
    new_val = arma::dot(new_val_vec, weights);

    diff = new_val - val;
    // we are now looking for the smallest value rather than largest so diff<0
    if (diff < 0) {
      // A_list = out3;
      idx_in = out2;
    }
    if (trace) Rcpp::Rcout << "\rIter " << i << ": " << diff << std::endl;
  }
  return (idx_in);
}

// [[Rcpp::export]]
arma::uvec uvec_minus(const arma::uvec &v, arma::uword rm_idx) {
  arma::uword n = v.size();
  if (rm_idx == 0) return v.tail(n-1);
  if (rm_idx == n-1) return v.head(n-1);
  arma::uvec res(v.size()-1);
  res.head(rm_idx) = v.head(rm_idx);
  res.tail(n-1-rm_idx) = v.tail(n-1-rm_idx);
  return res;
}

class HillClimbing {
private:
  const arma::vec C_list_;
  const arma::mat X_list_;
  const arma::mat sig_list_;
  const arma::vec weights_;
  const arma::uword nlist_;

  const arma::uword C_nrows_;
  const arma::uword X_nrows_;
  const arma::uword sig_nrows_;

  const arma::uword n_;
  const arma::uword N_;

public:
  arma::uvec idx_in_;
  arma::uvec idx_out_;
  arma::vec new_val_vec_;
  arma::vec best_val_vec_;
  double val_;
  double new_val_;

  arma::mat A_list_;
  arma::mat A_list_sub_;
  arma::mat rm1A_list_;

  arma::mat M_list_;
  arma::vec u_list_;

  const arma::uword A_nrows_;
  const arma::uword rm1A_nrows_;

  const arma::uword M_nrows_;
  const arma::uword u_nrows_;

  const arma::uword nfix_;
  const arma::uword rd_mode_; // robust designs mode: 1 == weighted, 2 == minimax.

  bool trace_;

public:
  HillClimbing(arma::uvec idx_in, arma::vec C_list, arma::mat X_list, arma::mat sig_list, arma::vec weights,
               arma::uword nfix = 0, arma::uword rd_mode = 0, bool trace=false) :
  C_list_(C_list), X_list_(X_list), sig_list_(sig_list), weights_(weights), nlist_(weights_.n_elem),
  C_nrows_(C_list.n_rows / nlist_), X_nrows_(X_list.n_rows / nlist_), sig_nrows_(sig_list.n_rows / nlist_),
  n_(idx_in.n_elem), N_(sig_nrows_), idx_in_(idx_in), idx_out_(N_ - n_, arma::fill::zeros),
  new_val_vec_(nlist_, arma::fill::zeros), best_val_vec_(nlist_, arma::fill::zeros), val_(0.0), new_val_(0.0),
  A_list_(n_ * nlist_, n_, arma::fill::zeros), A_list_sub_(n_ * nlist_, n_, arma::fill::zeros),
  rm1A_list_((n_-1) * nlist_, n_-1, arma::fill::zeros),
  M_list_(C_nrows_ * nlist_, C_nrows_, arma::fill::zeros),
  u_list_(X_nrows_ * nlist_, arma::fill::zeros),
  A_nrows_(n_), rm1A_nrows_(n_-1), M_nrows_(C_nrows_), u_nrows_(X_nrows_),
  nfix_(nfix), rd_mode_(rd_mode), trace_(trace) {
    Update_A_list();
    Update_M_list();
    Update_u_list();
    cout << "rd_mode = " << rd_mode_ << endl;
  }

  arma::vec get_C(arma::uword i) const {
    return C_list_.subvec(i * C_nrows_, (i + 1) * C_nrows_ - 1);
  }
  arma::mat get_X(arma::uword i) const {
    return X_list_.rows(i * X_nrows_, (i + 1) * X_nrows_ - 1);
  }
  arma::mat get_X(arma::uword i, const arma::uvec &idx_vec) const {
    return X_list_.rows(i * X_nrows_ + idx_vec);
  }
  arma::mat get_sig(arma::uword i) const {
    return sig_list_.rows(i * sig_nrows_, (i + 1) * sig_nrows_ - 1);
  }
  double get_sig(arma::uword i, arma::uword ii) const {
    return sig_list_.rows(i * sig_nrows_, (i + 1) * sig_nrows_ - 1).operator()(ii,ii);
  }
  arma::mat get_A(arma::uword i) const {
    return A_list_.rows(i * A_nrows_, (i + 1) * A_nrows_ - 1);
  }
  arma::mat get_Asub(arma::uword i) const {
    return A_list_sub_.rows(i * A_nrows_, (i + 1) * A_nrows_ - 1);
  }
  arma::mat get_rm1A(arma::uword i) const {
    return rm1A_list_.rows(i * rm1A_nrows_, (i + 1) * rm1A_nrows_ - 1);
  }
  arma::mat get_M(arma::uword i) const {
    return M_list_.rows(i * M_nrows_, (i + 1) * M_nrows_ - 1);
  }
  arma::vec get_u(arma::uword i) const {
    return u_list_.subvec(i * u_nrows_, (i + 1) * u_nrows_ - 1);
  }
  arma::vec get_u(arma::uword i, const arma::uvec &idx_vec) const {
    return u_list_.rows(i * u_nrows_ + idx_vec);
  }
  arma::uvec join_idx(const arma::uvec &idx, arma::uword elem) {
    return arma::join_cols(idx, arma::uvec({elem}));
  }

  double comp_c_obj_fun() {
    for (arma::uword j = 0; j < nlist_; ++j) {
      new_val_vec_(j) = c_obj_fun(get_M(j), get_C(j));
    }
    return rd_mode_ == 1 ? arma::dot(new_val_vec_, weights_) : new_val_vec_.max();
  }
  double check_all_neighbours(double &diff) {
    // if reordering doesn't find a better solution then check all the neighbours
    // check optimality and see if there are any neighbours that improve the solution
    // I have got it to go through in order of the ordered observations, and then break the loop if
    // it reaches one with a positive gradient because we know the ones after that will be worse

    if (trace_) Rcpp::Rcout << "Checking optimality..." << std::endl;
    for (arma::uword obs = nfix_; obs < n_; ++obs)  {
      if (trace_) Rcpp::Rcout << "\rChecking neighbour block: " << obs+1 << " of " << n_;
      Update_rm1A_list(obs);

      arma::uvec inx_in_no_obs = uvec_minus(idx_in_, obs);
      bool idx_in_updated = false;
      for (arma::uword obs_j = 0; obs_j < N_ - n_; ++obs_j) {
        if (idx_in_updated) {
          inx_in_no_obs = uvec_minus(idx_in_, obs);
          idx_in_updated = false;
        }
        const arma::uvec inx_out_no_obsj = uvec_minus(idx_out_, obs_j);

        arma::vec val_in_mat(nlist_, arma::fill::zeros);
        for (arma::uword idx = 0; idx < nlist_; ++idx) {
          val_in_mat(idx) = add_one_helper(get_rm1A(idx),
                     idx, inx_in_no_obs, idx_out_(obs_j));
        } // for loop
        double val_in = rd_mode_==1 ? arma::sum(val_in_mat * weights_) : val_in_mat.max();

        if (val_ - val_in < 0) {
          Update_A_list_sub(inx_in_no_obs, idx_out_(obs_j));
          bool flag = Update_M_list(join_idx(inx_in_no_obs, idx_out_(obs_j)), true, 1);
          if (!flag) continue;
          new_val_ = comp_c_obj_fun();

          if (new_val_ - val_ < 0) {
            diff = new_val_ - val_;
            if (trace_) Rcpp::Rcout << "\nImprovement found: " << new_val_ << std::endl;
            UpdateResult(inx_in_no_obs, inx_out_no_obsj, obs, obs_j);
            idx_in_updated = true;
            break;
          } // if (new_val - val < 0)
        } else break; // if (val - val_in < 0)
      } // for loop obs_j

      if (diff < 0) break;
    } // for loop obs
    return diff;
  }

  bool easy_swap(int &i, double &diff, bool reorder = false) {
    //++i;

    // if no easy swaps can be made, reorder the list and try the top one again
    if (reorder) reorder_obs();

    /*
    Update_rm1A_list(0);
    Update_A_list_sub(idx_in_.tail(n_-1), idx_out_(0));
    Update_M_list(join_idx(idx_in_.tail(n_-1), idx_out_(0)), true, 2);
    new_val_ = comp_c_obj_fun();
    diff = new_val_ - val_;
    */
    arma::uvec idx_in_m1 = uvec_minus(idx_in_, nfix_);
    Update_rm1A_list(nfix_);
    Update_A_list_sub(idx_in_m1, idx_out_(0));
    Update_M_list(join_idx(idx_in_m1, idx_out_(0)), true, 2);
    new_val_ = comp_c_obj_fun();
    diff = new_val_ - val_;

    // we are now looking for the smallest value rather than largest so diff < 0
    std::string str = reorder ? "Iter (reorder)" : "Iter ";
    if (trace_) Rcpp::Rcout << str << i << ": " << val_ << std::endl;

    bool flag = diff < 0;
    //if (flag) UpdateResult(idx_in_.rows(1, n_-1), idx_out_.rows(1, n_-1));
    if (flag) UpdateResult(idx_in_m1, idx_out_.rows(1, n_-1));

    return flag;
  }

  void grad_robust_step() {
    new_val_ = comp_c_obj_fun();

    idx_in_.head(20).t().print("DEBUG1: ");

    reorder_obs();

    idx_in_.head(20).t().print("DEBUG2: ");

    double diff = -1.0;
    int i = 0;
    // we now need diff to be negative
    while (diff < 0) {
      ++i;

      val_ = new_val_;
      best_val_vec_ = new_val_vec_;

      bool diff_l0 = easy_swap(i, diff);
      idx_in_.head(20).t().print("DEBUG3: ");

      if (!diff_l0) {
        // if no easy swaps can be made, reorder the list
        // and try the top one again
        diff_l0 = easy_swap(i, diff, true);
        // if reordering doesn't find a better solution
        // then check all the neighbours
        if (!diff_l0) check_all_neighbours(diff);
      }
    } // while loop
    idx_in_.head(20).t().print("DEBUG4: ");
  }

  double choose_swap_robust(arma::uvec &idx_in_tmp) {
    idx_in_tmp = idx_in_;
    init_idx_out();

    // find one index from idx_in to remove
    // which results in largest val of remove_one()
    const arma::vec val_out_vec = comp_val_out_vec();
    int idx_rm = nfix_ + val_out_vec.tail(n_ - nfix_).index_max();
    //int idx_rm = val_out_vec.index_max();

    // compute the new A without the index of idx_rm
    Update_rm1A_list(idx_rm);

    // remove index idx_rm from idx_in
    idx_in_.shed_row(idx_rm);

    // find one index from idx_out to add (swap)
    // which results in largest val of add_one()
    const arma::vec val_in_vec = comp_val_in_vec(true);
    arma::uword indexmax = val_in_vec.index_max();

    arma::uword idx_swap = idx_out_(indexmax);

    // compute the new A with the index of idx_swap
    Update_A_list_sub(idx_in_, idx_swap);
    idx_in_ = join_idx(idx_in_, idx_swap);
    idx_in_.swap(idx_in_tmp);

    return val_in_vec(indexmax);
  }

  void grad_robust_alg1() {
    new_val_ = comp_c_obj_fun();

    double diff = -1.0;
    int i = 0;
    // we now need diff to be negative
    while (diff < 0) {
      ++i;

      val_ = new_val_;
      arma::uvec out2;
      new_val_ = choose_swap_robust(out2);

      //we have to now recalculate all the lists of matrices for the new design proposed by the swap
      A_list_ = A_list_sub_;
      Update_M_list(out2, false, 2);
      Update_u_list();
      new_val_ = comp_c_obj_fun();

      diff = new_val_ - val_;
      // we are now looking for the smallest value rather than largest so diff<0
      if (diff < 0) idx_in_ = out2;
      if (trace_) Rcpp::Rcout << "\rIter " << i << ": " << diff << std::endl;
    }
  }

private:
  void Update_A_list() {
    for (arma::uword i = 0; i < nlist_; ++i) {
      const arma::mat tmp = sig_list_.submat(i * sig_nrows_ + idx_in_, idx_in_);
      A_list_.rows(i * A_nrows_, (i + 1) * A_nrows_ - 1) = tmp.i();
    }
  }

  void Update_A_list_sub(const arma::uvec &indices, arma::uword ii) {
    for (arma::uword i = 0; i < nlist_; ++i) {
      const arma::mat Asub = add_one_mat(get_rm1A(i), get_sig(i, ii),
                                         sig_list_.submat(i * sig_nrows_ + indices, arma::uvec({ii})));
      A_list_sub_.rows(i * A_nrows_, (i + 1) * A_nrows_ - 1) = Asub;


    }
  }

  void Update_rm1A_list(arma::uword obs) {
    for (arma::uword idx = 0; idx < nlist_; ++idx) {
      const arma::mat rm1A = remove_one_mat(get_A(idx), obs);
      rm1A_list_.rows(idx * rm1A_nrows_, (idx + 1) * rm1A_nrows_ - 1) = rm1A;
    }
  }

  // idx_vec ---- user defined indices to get X; if empty, idx_in_ will be used
  // useAsub ---- Should we get A from Asub_list_?
  // check_M == 0 ---- No positive definite check for M is needed
  //         == 1 ---- Do the check, return false if find not positive definite
  //         == 2 ---- Do the check, stop function if find not positive definite
  bool Update_M_list(const arma::uvec &idx_vec = {},
                     bool useAsub = false, int check_M = 0) {
    for (arma::uword i = 0; i < nlist_; ++i) {
      const arma::mat X = idx_vec.empty() ?
                          get_X(i, idx_in_) : get_X(i, idx_vec);
      const arma::mat A = useAsub? get_Asub(i) : get_A(i);
      const arma::mat M = X.t() * A * X;
      M_list_.rows(i * M_nrows_, (i + 1) * M_nrows_ - 1) = M;
      if (check_M && !M.is_sympd()) {
        if (check_M == 1) return false;
        else {
          const arma::vec colsum = arma::trans(arma::sum(M, 0));
          const arma::uvec colidx = arma::find(colsum == 0);
          Rcpp::Rcout << "ERROR MSG:\n"
                      << "Design " << i+1 << " has the following column(s) with zero colsum:\n"
                      << colidx.t() + 1 << std::endl;
          Rcpp::stop("M not positive semi-definite.");
        }
      }
    }
    return true;
  }

  void Update_u_list() {
    for (arma::uword i = 0; i < nlist_; ++i) {
      const arma::mat X = get_X(i);
      const arma::vec C = get_C(i);
      const arma::mat M_inv = arma::inv_sympd(get_M(i));
      u_list_.subvec(i * u_nrows_, (i + 1) * u_nrows_ - 1) = X * (M_inv * C);
    }
  }

  void UpdateResult(const arma::uvec &in_vec, const arma::uvec &out_vec,
                    arma::uword obs = 0, arma::uword obs_j = 0) {
    const arma::uword swap_idx = idx_in_(obs);
    idx_in_ = join_idx(in_vec, idx_out_(obs_j));
    idx_out_ = join_idx(out_vec, swap_idx);
    A_list_ = A_list_sub_;
    Update_u_list();
  }

  void init_idx_out() {
    // generate the complete index
    arma::vec idx = arma::linspace(0, u_nrows_ - 1, u_nrows_);
    arma::uvec uidx = arma::conv_to<arma::uvec>::from(idx);

    // get the index not included in complete index
    idx_out_ =  std_setdiff(uidx, idx_in_);
  }

  // evaluate remove_one
  arma::vec comp_val_out_vec() const {
    arma::mat val_out_mat(idx_in_.n_elem, nlist_, arma::fill::zeros);
#pragma omp parallel for
    for (std::size_t j = 0; j < nlist_; ++j) {
      const arma::mat A = get_A(j);
      const arma::vec u_idx_in = get_u(j, idx_in_);
      for (std::size_t i = 0; i < idx_in_.n_elem; ++i) {
        val_out_mat(i,j) = remove_one(A, i, u_idx_in);
      }
    }

    return rd_mode_==1 ? val_out_mat * weights_ : arma::vec(arma::max(val_out_mat, 1));
  }

  double add_one_helper(const arma::mat &A, arma::uword i,
                        const arma::uvec &idx_vec, arma::uword ii) {
    arma::uvec sig_ridx = i * sig_nrows_ + idx_vec;
    arma::uvec u_ridx = i * u_nrows_ + join_idx(idx_vec, ii);
    return add_one(A, get_sig(i, ii),
                   sig_list_.submat(sig_ridx, arma::uvec({ii})),
                   u_list_.rows(u_ridx));
  }

  // evaluate add_one
  arma::vec comp_val_in_vec(bool use_rm1A = false) {
    arma::mat val_in_mat(idx_out_.n_elem,nlist_,arma::fill::zeros);
#pragma omp parallel for
    for (arma::uword j = 0; j < nlist_; ++j) {
      const arma::mat A = use_rm1A ? get_rm1A(j) : get_A(j);
      for (arma::uword i = 0; i < idx_out_.n_elem; ++i) {
        val_in_mat(i,j) = add_one_helper(A, j, idx_in_, idx_out_(i));
      }
    }

    return rd_mode_==1 ? val_in_mat * weights_ : arma::vec(arma::max(val_in_mat, 1));
  }

  void reorder_obs() {
    init_idx_out();

    // find one index from idx_in to remove
    // which results in largest val of remove_one()
    arma::vec val_out_vec = comp_val_out_vec();
    arma::uvec indices1 = arma::zeros<uvec>(n_);
    if (nfix_ != 0) {
	for (arma::uword i = 0; i < nfix_; ++i) { indices1[i] = i; }
	indices1.tail(n_ - nfix_) = nfix_ + arma::sort_index(val_out_vec.tail(n_ - nfix_), "descend");
    } else {
	indices1 = arma::sort_index(val_out_vec, "descend");
    }

    // arma::uvec indices1 = arma::sort_index(val_out_vec, "descend");
    idx_in_ = idx_in_.rows(indices1);

    for (arma::uword j = 0; j < nlist_; ++j) {
      A_list_.rows(j*A_nrows_, (j+1)*A_nrows_-1)
        = A_list_.submat(j*A_nrows_ + indices1, indices1);
    }

    // find one index from idx_out to add (swap)
    // which results in largest val of add_one()
    arma::vec val_in_vec = comp_val_in_vec();
    arma::uvec indices2 = arma::sort_index(val_in_vec, "descend");
    idx_out_ = idx_out_.rows(indices2);
  }
};

// [[Rcpp::export]]
Rcpp::List GradRobustStep(arma::uvec idx_in, arma::vec C_list, arma::mat X_list, arma::mat sig_list, arma::vec weights,
			  arma::uword nfix = 0, arma::uword rd_mode = 1) {
  HillClimbing hc(idx_in, C_list, X_list, sig_list, weights, nfix, rd_mode, true);
  hc.grad_robust_step();
  return Rcpp::List::create(Named("idx_in") = hc.idx_in_,
                            Named("idx_out") = hc.idx_out_,
                            Named("best_val_vec") = hc.best_val_vec_);
}

// [[Rcpp::export]]
Rcpp::List GradRobustAlg1(arma::uvec idx_in, arma::vec C_list, arma::mat X_list, arma::mat sig_list, arma::vec weights,
			  arma::uword nfix = 0, arma::uword rd_mode = 1) {
  HillClimbing hc(idx_in, C_list, X_list, sig_list, weights, nfix, rd_mode, true);
  hc.grad_robust_alg1();
  return Rcpp::List::create(Named("idx_in") = hc.idx_in_,
                            Named("idx_out") = hc.idx_out_,
                            Named("best_val_vec") = hc.best_val_vec_);
}
