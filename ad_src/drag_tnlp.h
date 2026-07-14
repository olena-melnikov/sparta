/* ----------------------------------------------------------------------
   drag_tnlp: an IPOPT TNLP that minimizes drag(alpha) for the symmetric
   Bezier body over the 4 design variables alpha = [x1,y1,x2,y2], subject
   to box constraints only (no general constraints, m = 0):

       a <= x1, x2 <= b     (x-coords of control points P1, P2)
       c <= y1, y2 <= d     (y-coords of control points P1, P2)

   The objective and its gradient come from the existing drag_objective
   API (drag_avg / grad_fd). Which gradient is used is a compile-time
   choice matching the linked SPARTA library:
     - AD build   (-DSPARTA_AD): drag_avg() returns value + gradient in a
       single SPARTA run (forward-mode AD).
     - stock build             : value from drag_avg(), gradient from the
       common-random-number central finite difference grad_fd().

   Each evaluation is expensive (1 SPARTA run for AD value+grad, nseeds
   runs for a stock value, 9*nseeds for a stock FD gradient), so results
   are cached per design point to avoid recomputation between the paired
   eval_f / eval_grad_f calls IPOPT makes.
------------------------------------------------------------------------- */

#ifndef DRAG_TNLP_H
#define DRAG_TNLP_H

#include "IpTNLP.hpp"
#include "drag_objective.h"

#include <vector>

// One recorded IPOPT iterate (for the output trajectory + progress bar).
struct TrajPoint {
  int    iter;
  double alpha[4];
  double drag;
  double inf_pr;   // primal infeasibility (0 here: bounds only)
  double inf_du;   // dual infeasibility (optimality measure)
};

class DragTNLP : public Ipopt::TNLP {
 public:
  DragTNLP(const double alpha0[4],
           const double xlo[4], const double xhi[4],
           const int *seeds, int nseeds,
           const DragCase &c, double h,
           int max_iter, bool show_bar, double obj_scale);

  // --- IPOPT TNLP interface -------------------------------------------
  bool get_nlp_info(Ipopt::Index &n, Ipopt::Index &m,
                    Ipopt::Index &nnz_jac_g, Ipopt::Index &nnz_h_lag,
                    IndexStyleEnum &index_style) override;

  bool get_bounds_info(Ipopt::Index n, Ipopt::Number *x_l, Ipopt::Number *x_u,
                       Ipopt::Index m, Ipopt::Number *g_l,
                       Ipopt::Number *g_u) override;

  // Scale the (tiny, ~1e-21) drag objective to O(1) so it is not swamped
  // by IPOPT's log-barrier term. obj_scaling = 1/|drag(start)|.
  bool get_scaling_parameters(Ipopt::Number &obj_scaling,
                              bool &use_x_scaling, Ipopt::Index n,
                              Ipopt::Number *x_scaling, bool &use_g_scaling,
                              Ipopt::Index m, Ipopt::Number *g_scaling) override;

  bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number *x,
                          bool init_z, Ipopt::Number *z_L, Ipopt::Number *z_U,
                          Ipopt::Index m, bool init_lambda,
                          Ipopt::Number *lambda) override;

  bool eval_f(Ipopt::Index n, const Ipopt::Number *x, bool new_x,
              Ipopt::Number &obj_value) override;

  bool eval_grad_f(Ipopt::Index n, const Ipopt::Number *x, bool new_x,
                   Ipopt::Number *grad_f) override;

  bool eval_g(Ipopt::Index n, const Ipopt::Number *x, bool new_x,
              Ipopt::Index m, Ipopt::Number *g) override;

  bool eval_jac_g(Ipopt::Index n, const Ipopt::Number *x, bool new_x,
                  Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index *iRow,
                  Ipopt::Index *jCol, Ipopt::Number *values) override;

  bool intermediate_callback(Ipopt::AlgorithmMode mode, Ipopt::Index iter,
                             Ipopt::Number obj_value, Ipopt::Number inf_pr,
                             Ipopt::Number inf_du, Ipopt::Number mu,
                             Ipopt::Number d_norm,
                             Ipopt::Number regularization_size,
                             Ipopt::Number alpha_du, Ipopt::Number alpha_pr,
                             Ipopt::Index ls_trials,
                             const Ipopt::IpoptData *ip_data,
                             Ipopt::IpoptCalculatedQuantities *ip_cq) override;

  void finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n,
                         const Ipopt::Number *x, const Ipopt::Number *z_L,
                         const Ipopt::Number *z_U, Ipopt::Index m,
                         const Ipopt::Number *g, const Ipopt::Number *lambda,
                         Ipopt::Number obj_value, const Ipopt::IpoptData *ip_data,
                         Ipopt::IpoptCalculatedQuantities *ip_cq) override;

  // --- results (read after OptimizeTNLP returns) ----------------------
  double init_alpha[4];
  double init_drag;
  double final_alpha[4];
  double final_drag;
  int    solve_status;          // Ipopt::SolverReturn as int
  std::vector<TrajPoint> traj;

 private:
  // Objective + gradient at x, cached. Returns drag(x); fills grad if
  // want_grad (grad may be NULL). Reuses the cache when x is unchanged.
  double evaluate(const double x[4], bool want_grad, double grad[4]);

  double a0_[4];
  double xl_[4], xu_[4];
  std::vector<int> seeds_;
  DragCase c_;
  double h_;
  int  max_iter_;
  bool show_bar_;
  double obj_scale_;

  // value/gradient cache keyed on the 4-vector design point
  bool   fcache_valid_;
  double fcache_x_[4];
  double fcache_f_;
  bool   gcache_valid_;
  double gcache_x_[4];
  double gcache_g_[4];

  bool got_init_;               // captured init_drag on first eval
};

#endif
