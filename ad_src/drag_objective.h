/* ----------------------------------------------------------------------
   drag_objective: drag(alpha) for the symmetric Bezier body, evaluated
   by running SPARTA in-process through its library interface
   (src/library.h). No SPARTA sources are modified.

   Layering (matches the AD roadmap):
     alpha --(bezier_geom, differentiable seam)--> pts,norms
           --(SPARTA run, seeded)---------------> drag realization

   AD/SAA notes:
     - seed is an explicit argument: sample-average approximation is
       mean_k drag(alpha, seed_k); reusing the SAME seeds across alpha
       values gives common random numbers (smoother finite differences,
       consistent SAA gradients).
     - callers only see alpha[4] -> double, so internals can later swap
       the temp-file read_surf for an in-memory create_surf command and
       return d(drag)/d(pts) without touching call sites.
------------------------------------------------------------------------- */

#ifndef DRAG_OBJECTIVE_H
#define DRAG_OBJECTIVE_H

struct DragCase {
  // geometry (fixed, non-differentiated except alpha)
  double chord = 4.0;          // body length, nose (0,0) -> tail (chord,0)
  int nseg = 25;               // Bezier segments per half

  // placement / domain (defaults match examples/circle scale)
  double origin[2] = {3.0, 5.0};
  double boxhi[2]  = {10.0, 10.0};
  int grid[2]      = {20, 20};

  // gas + flow
  double vstream = 100.0;
  double nrho = 1.0, fnum = 0.001, tstep = 0.0001;
  const char *species_file = "air.species";   // must exist in cwd
  const char *vss_file     = "air.vss";

  // run protocol
  int nsettle = 500;   // steps to develop the flow (excluded from avg)
  int navg    = 500;   // steps drag is time-averaged over

  int verbose = 0;     // 1 = show SPARTA screen output
  int progress = 0;    // 1 = print one line per SPARTA run to stderr

  // physics-model switches (validation ladder: specular + nocoll gives
  // free-molecular flow with deterministic trajectories, where pathwise
  // AD and small-h CRN-FD must agree)
  int specular = 0;    // 1 = specular wall instead of diffuse
  int collisions = 1;  // 0 = no gas-gas collisions (free-molecular)
};

// One DSMC realization of drag (x-force on the body, summed over all
// surf elements, time-averaged over the navg window).
// If fy is non-NULL, also returns the side force (symmetry check).
// If grad is non-NULL (length 4) and this is the AD build, receives
// the forward-mode gradient d(drag)/d(alpha) of this realization;
// in the stock build grad is filled with NaN (use drag_has_ad()).
// Exits with an error message on failure.
double drag(const double alpha[4], int seed, const DragCase &c,
            double *fy = 0, double *grad = 0);

// 1 if compiled against the AD build (-DSPARTA_AD): drag() can fill
// grad. 0 in the stock build.
int drag_has_ad();

// Sample-average approximation over nseeds seeds.
// If out is non-NULL, out[k] receives the k-th realization.
// If grad is non-NULL (AD build), receives the seed-averaged gradient.
double drag_avg(const double alpha[4], const int *seeds, int nseeds,
                const DragCase &c, double *out = 0, double *grad = 0);

// Central finite-difference gradient of drag_avg with COMMON RANDOM
// NUMBERS: every perturbed evaluation reuses the same seed set, so the
// Monte Carlo noise largely cancels in the differences. This is the
// frozen-stream (pathwise) derivative -- the same quantity forward-mode
// AD will compute, hence the Phase 4 acceptance reference.
//   h    : absolute perturbation per alpha component
//   grad : out, length 4
// Returns drag_avg(alpha) (the unperturbed value). Cost: 9 x nseeds runs.
double grad_fd(const double alpha[4], const int *seeds, int nseeds,
               const DragCase &c, double h, double grad[4]);

#endif
