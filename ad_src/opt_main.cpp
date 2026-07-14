/* ----------------------------------------------------------------------
   opt_main: Phase 1 driver -- CRN-FD gradient of drag(alpha) and a
   simple projected-gradient-descent shape optimization.

   Usage (run from a directory with air.species/air.vss, e.g.
   bezier/, now nested under ad_src/):

     gradient only:
       opt_run --alpha x1,y1,x2,y2 [--seeds s1,s2,...] [--h H] [...]

     optimization:
       opt_run --alpha x1,y1,x2,y2 --iters N [--eta STEP] [...]

   Options:
     --alpha   start / evaluation point            [required]
     --seeds   CRN seed set (default 12345,777,31415)
     --h       FD step (default 0.05)
     --iters   descent iterations (default 0 = gradient only)
     --eta     descent step length in alpha-space (default 0.15);
               update is alpha -= eta * g/|g| (normalized: drag's
               absolute scale ~1e-21 makes raw gradient steps useless)
     --nsettle/--navg/--nseg/--chord  forwarded to DragCase
     --verbose SPARTA screen output

   Bounds (projected after each step, keeps the body valid & in box):
     x1,x2 in [0.2, chord-0.2],  y1,y2 in [0.05, 3.0]

   The gradients printed here are the acceptance reference for the AD
   build (Phase 4): AD must reproduce them at the same seeds within FD
   truncation error, since both compute the frozen-stream derivative.
------------------------------------------------------------------------- */

#include "drag_objective.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

static const char *USAGE =
  "usage: opt_run --alpha x1,y1,x2,y2 [--seeds s1,s2,...] [--h H]\n"
  "               [--iters N] [--eta STEP] [--nsettle N] [--navg N]\n"
  "               [--nseg N] [--chord L] [--verbose]\n";

static double clip(double v, double lo, double hi)
{
  return v < lo ? lo : (v > hi ? hi : v);
}

static void clamp(double a[4], double chord)
{
  const double xmin = 0.2, ymin = 0.05, ymax = 3.0;
  double xmax = chord - 0.2;
  a[0] = clip(a[0],xmin,xmax);
  a[2] = clip(a[2],xmin,xmax);
  a[1] = clip(a[1],ymin,ymax);
  a[3] = clip(a[3],ymin,ymax);
}

int main(int narg, char **arg)
{
  std::vector<double> alpha;
  std::vector<int> seeds;
  DragCase c;
  double h = 0.05, eta = 0.15;
  int iters = 0, quiet = 0;

  for (int i = 1; i < narg; i++) {
    if (!strcmp(arg[i],"--alpha") && i+1 < narg) {
      char *tok = strtok(arg[++i],",");
      while (tok) { alpha.push_back(atof(tok)); tok = strtok(NULL,","); }
    } else if (!strcmp(arg[i],"--seeds") && i+1 < narg) {
      char *tok = strtok(arg[++i],",");
      while (tok) { seeds.push_back(atoi(tok)); tok = strtok(NULL,","); }
    } else if (!strcmp(arg[i],"--h") && i+1 < narg) {
      h = atof(arg[++i]);
    } else if (!strcmp(arg[i],"--iters") && i+1 < narg) {
      iters = atoi(arg[++i]);
    } else if (!strcmp(arg[i],"--eta") && i+1 < narg) {
      eta = atof(arg[++i]);
    } else if (!strcmp(arg[i],"--nsettle") && i+1 < narg) {
      c.nsettle = atoi(arg[++i]);
    } else if (!strcmp(arg[i],"--navg") && i+1 < narg) {
      c.navg = atoi(arg[++i]);
    } else if (!strcmp(arg[i],"--nseg") && i+1 < narg) {
      c.nseg = atoi(arg[++i]);
    } else if (!strcmp(arg[i],"--chord") && i+1 < narg) {
      c.chord = atof(arg[++i]);
    } else if (!strcmp(arg[i],"--vstream") && i+1 < narg) {
      c.vstream = atof(arg[++i]);
    } else if (!strcmp(arg[i],"--specular")) {
      c.specular = 1;
    } else if (!strcmp(arg[i],"--nocoll")) {
      c.collisions = 0;
    } else if (!strcmp(arg[i],"--verbose")) {
      c.verbose = 1;
    } else if (!strcmp(arg[i],"--quiet")) {
      quiet = 1;
    } else { fprintf(stderr,"%s",USAGE); return 1; }
  }

  if (alpha.size() != 4) {
    fprintf(stderr,"ERROR: --alpha needs exactly 4 values\n%s",USAGE);
    return 1;
  }
  if (seeds.empty()) { seeds.push_back(12345); seeds.push_back(777);
                       seeds.push_back(31415); }

  if (!quiet) c.progress = 1;

  double a[4] = {alpha[0],alpha[1],alpha[2],alpha[3]};
  double g[4];

  if (iters == 0) {
    double d = grad_fd(a,seeds.data(),(int)seeds.size(),c,h,g);
    printf("alpha = %.6g %.6g %.6g %.6g\n",a[0],a[1],a[2],a[3]);
    printf("drag  = %.8e   (CRN mean over %d seeds)\n",d,(int)seeds.size());
    printf("grad  = %.6e %.6e %.6e %.6e   (h=%g)\n",g[0],g[1],g[2],g[3],h);
    return 0;
  }

  printf("# projected gradient descent: eta=%g h=%g seeds=%d "
         "nsettle=%d navg=%d\n",eta,h,(int)seeds.size(),c.nsettle,c.navg);
  printf("# iter    drag            |grad|        alpha\n");

  for (int it = 0; it <= iters; it++) {
    double d = grad_fd(a,seeds.data(),(int)seeds.size(),c,h,g);
    double gn = sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]+g[3]*g[3]);
    printf("%4d      %.8e  %.4e    %.5g %.5g %.5g %.5g\n",
           it,d,gn,a[0],a[1],a[2],a[3]);
    fflush(stdout);
    if (it == iters || gn == 0.0) break;
    for (int j = 0; j < 4; j++) a[j] -= eta*g[j]/gn;
    clamp(a,c.chord);
  }
  return 0;
}
