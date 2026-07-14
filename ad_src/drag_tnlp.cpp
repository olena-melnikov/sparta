/* ----------------------------------------------------------------------
   drag_tnlp: IPOPT TNLP implementation. See drag_tnlp.h for the design.
------------------------------------------------------------------------- */

#include "drag_tnlp.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <unistd.h>   // isatty

using namespace Ipopt;

// ---- colorful progress bar (rendered on stderr) -----------------------

// 256-color code sampled along a cyan -> green -> yellow -> red ramp.
static int bar_color(double t)
{
  static const int ramp[] = {
    51, 50, 49, 48, 47, 46, 82, 118, 154, 190, 226, 220, 214, 208, 202, 196
  };
  const int n = (int)(sizeof(ramp) / sizeof(ramp[0]));
  int i = (int)(t * (n - 1) + 0.5);
  if (i < 0) i = 0;
  if (i >= n) i = n - 1;
  return ramp[i];
}

static void draw_bar(int iter, int max_iter, double drag, double inf_du)
{
  const int width = 30;
  double frac = (max_iter > 0) ? (double)iter / (double)max_iter : 0.0;
  if (frac > 1.0) frac = 1.0;
  int fill = (int)(frac * width + 0.5);

  fprintf(stderr, "\r\033[1m[\033[0m");
  for (int i = 0; i < width; i++) {
    if (i < fill) {
      double t = (width > 1) ? (double)i / (double)(width - 1) : 0.0;
      fprintf(stderr, "\033[38;5;%dm\xE2\x96\x88", bar_color(t));  // U+2588 block
    } else {
      fprintf(stderr, "\033[38;5;240m\xE2\x94\x80");               // U+2500 dash
    }
  }
  fprintf(stderr, "\033[0m\033[1m]\033[0m %3d%%  iter %d/%d  "
                  "\033[36mdrag=%.4e\033[0m  opt=%.1e   ",
          (int)(frac * 100.0), iter, max_iter, drag, inf_du);
  fflush(stderr);
}

// ---- helpers ----------------------------------------------------------

static bool same4(const double *a, const double *b)
{
  return a[0] == b[0] && a[1] == b[1] && a[2] == b[2] && a[3] == b[3];
}

// ---- construction -----------------------------------------------------

DragTNLP::DragTNLP(const double alpha0[4],
                   const double xlo[4], const double xhi[4],
                   const int *seeds, int nseeds,
                   const DragCase &c, double h,
                   int max_iter, bool show_bar, double obj_scale)
  : init_drag(0.0), final_drag(0.0), solve_status(-1),
    c_(c), h_(h), max_iter_(max_iter), show_bar_(show_bar),
    obj_scale_(obj_scale),
    fcache_valid_(false), fcache_f_(0.0), gcache_valid_(false),
    got_init_(false)
{
  for (int i = 0; i < 4; i++) {
    a0_[i]  = alpha0[i];
    xl_[i]  = xlo[i];
    xu_[i]  = xhi[i];
    init_alpha[i]  = alpha0[i];
    final_alpha[i] = alpha0[i];
  }
  for (int k = 0; k < nseeds; k++) seeds_.push_back(seeds[k]);
}

// ---- objective + gradient with caching --------------------------------

double DragTNLP::evaluate(const double x[4], bool want_grad, double grad[4])
{
  if (want_grad && gcache_valid_ && same4(x, gcache_x_)) {
    memcpy(grad, gcache_g_, 4 * sizeof(double));
    if (fcache_valid_ && same4(x, fcache_x_)) return fcache_f_;
  }
  if (!want_grad && fcache_valid_ && same4(x, fcache_x_)) return fcache_f_;

  const int ns = (int)seeds_.size();
  double v, g[4];

#ifdef SPARTA_AD
  // AD build: one SPARTA run yields value AND gradient (gradient is free).
  v = drag_avg(x, seeds_.data(), ns, c_, /*out=*/0, g);
  memcpy(fcache_x_, x, 4 * sizeof(double)); fcache_f_ = v; fcache_valid_ = true;
  memcpy(gcache_x_, x, 4 * sizeof(double));
  memcpy(gcache_g_, g, 4 * sizeof(double)); gcache_valid_ = true;
  if (want_grad) memcpy(grad, g, 4 * sizeof(double));
#else
  // Stock build: value from drag_avg (nseeds runs); gradient from the
  // CRN central finite difference (9*nseeds runs) only when requested.
  if (want_grad) {
    v = grad_fd(x, seeds_.data(), ns, c_, h_, g);   // returns unperturbed value
    memcpy(fcache_x_, x, 4 * sizeof(double)); fcache_f_ = v; fcache_valid_ = true;
    memcpy(gcache_x_, x, 4 * sizeof(double));
    memcpy(gcache_g_, g, 4 * sizeof(double)); gcache_valid_ = true;
    memcpy(grad, g, 4 * sizeof(double));
  } else {
    v = drag_avg(x, seeds_.data(), ns, c_, 0, 0);
    memcpy(fcache_x_, x, 4 * sizeof(double)); fcache_f_ = v; fcache_valid_ = true;
  }
#endif

  if (!got_init_) { init_drag = v; got_init_ = true; }
  return v;
}

// ---- TNLP interface ---------------------------------------------------

bool DragTNLP::get_nlp_info(Index &n, Index &m, Index &nnz_jac_g,
                            Index &nnz_h_lag, IndexStyleEnum &index_style)
{
  n = 4;
  m = 0;
  nnz_jac_g = 0;
  nnz_h_lag = 0;               // Hessian via limited-memory, not supplied
  index_style = TNLP::C_STYLE;
  return true;
}

bool DragTNLP::get_bounds_info(Index n, Number *x_l, Number *x_u,
                               Index m, Number *g_l, Number *g_u)
{
  (void)m; (void)g_l; (void)g_u;
  for (Index i = 0; i < n; i++) { x_l[i] = xl_[i]; x_u[i] = xu_[i]; }
  return true;
}

bool DragTNLP::get_scaling_parameters(Number &obj_scaling, bool &use_x_scaling,
                                      Index n, Number *x_scaling,
                                      bool &use_g_scaling, Index m,
                                      Number *g_scaling)
{
  (void)n; (void)x_scaling; (void)m; (void)g_scaling;
  obj_scaling   = obj_scale_;   // 1/|drag(start)| -> scaled objective ~ O(1)
  use_x_scaling = false;        // design vars are already O(1)
  use_g_scaling = false;        // no constraints
  return true;
}

bool DragTNLP::get_starting_point(Index n, bool init_x, Number *x,
                                  bool init_z, Number *z_L, Number *z_U,
                                  Index m, bool init_lambda, Number *lambda)
{
  (void)init_z; (void)z_L; (void)z_U; (void)m; (void)init_lambda; (void)lambda;
  if (init_x) for (Index i = 0; i < n; i++) x[i] = a0_[i];
  return true;
}

bool DragTNLP::eval_f(Index n, const Number *x, bool new_x, Number &obj_value)
{
  (void)n; (void)new_x;
  obj_value = evaluate(x, /*want_grad=*/false, 0);
  return true;
}

bool DragTNLP::eval_grad_f(Index n, const Number *x, bool new_x, Number *grad_f)
{
  (void)new_x;
  double g[4];
  evaluate(x, /*want_grad=*/true, g);
  for (Index i = 0; i < n; i++) grad_f[i] = g[i];
  return true;
}

bool DragTNLP::eval_g(Index n, const Number *x, bool new_x, Index m, Number *g)
{
  (void)n; (void)x; (void)new_x; (void)m; (void)g;
  return true;                 // m = 0, nothing to do
}

bool DragTNLP::eval_jac_g(Index n, const Number *x, bool new_x, Index m,
                          Index nele_jac, Index *iRow, Index *jCol,
                          Number *values)
{
  (void)n; (void)x; (void)new_x; (void)m; (void)nele_jac;
  (void)iRow; (void)jCol; (void)values;
  return true;                 // no constraint Jacobian
}

bool DragTNLP::intermediate_callback(AlgorithmMode mode, Index iter,
                                     Number obj_value, Number inf_pr,
                                     Number inf_du, Number mu, Number d_norm,
                                     Number regularization_size, Number alpha_du,
                                     Number alpha_pr, Index ls_trials,
                                     const IpoptData *ip_data,
                                     IpoptCalculatedQuantities *ip_cq)
{
  (void)mode; (void)mu; (void)d_norm; (void)regularization_size;
  (void)alpha_du; (void)alpha_pr; (void)ls_trials;

  TrajPoint p;
  p.iter   = (int)iter;
  p.drag   = obj_value;
  p.inf_pr = inf_pr;
  p.inf_du = inf_du;

  double x[4];
  bool ok = get_curr_iterate(ip_data, ip_cq, /*scaled=*/false, 4, x,
                             0, 0, 0, 0, 0);
  if (ok) memcpy(p.alpha, x, 4 * sizeof(double));
  else if (!traj.empty()) memcpy(p.alpha, traj.back().alpha, 4 * sizeof(double));
  else memcpy(p.alpha, a0_, 4 * sizeof(double));
  traj.push_back(p);

  if (show_bar_) draw_bar((int)iter, max_iter_, obj_value, inf_du);
  return true;
}

void DragTNLP::finalize_solution(SolverReturn status, Index n, const Number *x,
                                 const Number *z_L, const Number *z_U, Index m,
                                 const Number *g, const Number *lambda,
                                 Number obj_value, const IpoptData *ip_data,
                                 IpoptCalculatedQuantities *ip_cq)
{
  (void)z_L; (void)z_U; (void)m; (void)g; (void)lambda;
  (void)ip_data; (void)ip_cq;
  solve_status = (int)status;
  final_drag = obj_value;
  for (Index i = 0; i < n; i++) final_alpha[i] = x[i];
  if (show_bar_) { fprintf(stderr, "\n"); fflush(stderr); }
}
