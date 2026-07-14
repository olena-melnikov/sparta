/* ----------------------------------------------------------------------
   bezier_gen: command-line driver for Path A of the AD roadmap.

   Generates a runnable SPARTA case (modeled on examples/circle) for a
   closed 2D body symmetric about its axis, built from two mirrored
   cubic Bezier halves (see BezierGeom::symmetric_body_to_lines):

     upper half: P0=(0,0), P1=alpha1, P2=alpha2, P3=(chord,0)
     lower half: mirror (y -> -y)

   Design vector alpha = x1,y1,x2,y2 in BODY coordinates (nose at the
   origin, tail at (chord,0), y1,y2 > 0). The body is translated by
   --origin into the simulation box.

   Outputs:
     data.<prefix> : SPARTA surf data file (clockwise, outward normals)
     in.<prefix>   : input deck (2d flow past the body + drag output)

   Usage:
     bezier_gen --alpha x1,y1,x2,y2 [options]

   Options (defaults match examples/circle where applicable):
     --alpha  <x1,y1,x2,y2>  free control points          [required]
     --chord  <L>      body length (default 4.0)
     --nseg   <int>    segments PER HALF, total 2*nseg (default 25)
     --origin <x> <y>  nose position in the box (default 3 5)
     --prefix <name>   -> in.<name>, data.<name> (default bezier_alpha)
     --boxhi  <x> <y>  domain upper corner, lower is 0 0 (default 10 10)
     --grid   <nx> <ny> grid cells (default 20 20)
     --vstream <vx>    free-stream x-velocity (default 100.0)
     --nsteps <int>    timesteps to run (default 1000)
     --seed   <int>    RNG seed (default 12345)

   This driver intentionally goes through a *file* interface: it is the
   validation baseline. The in-memory create_surf command (Path B) must
   reproduce these results bit-for-bit at the same seed.
------------------------------------------------------------------------- */

#include "bezier_geom.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

static void usage_and_exit()
{
  fprintf(stderr,
    "usage: bezier_gen --alpha x1,y1,x2,y2 [--chord L] [--nseg N]\n"
    "                  [--origin X Y] [--prefix name] [--boxhi X Y]\n"
    "                  [--grid NX NY] [--vstream VX] [--nsteps N] [--seed S]\n");
  exit(1);
}

int main(int narg, char **arg)
{
  std::vector<double> alpha;
  double chord = 4.0;
  int nseg = 25;
  double origin[2] = {3.0, 5.0};
  std::string prefix = "bezier_alpha";
  double boxhi[2] = {10.0, 10.0};
  int grid[2] = {20, 20};
  double vstream = 100.0;
  int nsteps = 1000;
  int seed = 12345;

  // ---- parse command line ----

  for (int i = 1; i < narg; i++) {
    if (!strcmp(arg[i],"--alpha") && i+1 < narg) {
      char *tok = strtok(arg[++i],",");
      while (tok) { alpha.push_back(atof(tok)); tok = strtok(NULL,","); }
    } else if (!strcmp(arg[i],"--chord") && i+1 < narg) {
      chord = atof(arg[++i]);
    } else if (!strcmp(arg[i],"--nseg") && i+1 < narg) {
      nseg = atoi(arg[++i]);
    } else if (!strcmp(arg[i],"--origin") && i+2 < narg) {
      origin[0] = atof(arg[++i]); origin[1] = atof(arg[++i]);
    } else if (!strcmp(arg[i],"--prefix") && i+1 < narg) {
      prefix = arg[++i];
    } else if (!strcmp(arg[i],"--boxhi") && i+2 < narg) {
      boxhi[0] = atof(arg[++i]); boxhi[1] = atof(arg[++i]);
    } else if (!strcmp(arg[i],"--grid") && i+2 < narg) {
      grid[0] = atoi(arg[++i]); grid[1] = atoi(arg[++i]);
    } else if (!strcmp(arg[i],"--vstream") && i+1 < narg) {
      vstream = atof(arg[++i]);
    } else if (!strcmp(arg[i],"--nsteps") && i+1 < narg) {
      nsteps = atoi(arg[++i]);
    } else if (!strcmp(arg[i],"--seed") && i+1 < narg) {
      seed = atoi(arg[++i]);
    } else usage_and_exit();
  }

  if (alpha.size() != 4) {
    fprintf(stderr,"ERROR: --alpha needs exactly 4 values: x1,y1,x2,y2\n");
    usage_and_exit();
  }
  if (nseg < 2) { fprintf(stderr,"ERROR: nseg must be >= 2\n"); exit(1); }
  if (chord <= 0.0) { fprintf(stderr,"ERROR: chord must be > 0\n"); exit(1); }
  if (alpha[1] <= 0.0 || alpha[3] <= 0.0)
    fprintf(stderr,"warning: y1,y2 should be > 0 for a valid "
                   "(non-self-intersecting) body\n");

  // ---- geometry: alpha -> closed clockwise loop in body coords ----

  int ntot = 2*nseg;                    // total segments
  std::vector<double> pts(2*(ntot+1)), norms(2*ntot);
  BezierGeom::symmetric_body_to_lines(alpha.data(),chord,nseg,
                                      pts.data(),norms.data());

  // sanity checks: closed, clockwise (outward normals), non-degenerate

  if (!BezierGeom::is_closed(ntot,pts.data(),1.0e-12)) {
    fprintf(stderr,"ERROR: internal: loop not closed\n"); exit(1);
  }
  if (BezierGeom::signed_area(ntot,pts.data()) >= 0.0) {
    fprintf(stderr,"ERROR: loop not clockwise; check alpha "
                   "(y1,y2 > 0 expected)\n"); exit(1);
  }
  double minlen = BezierGeom::min_segment_length(ntot,pts.data());
  if (minlen <= 0.0) {
    fprintf(stderr,"ERROR: degenerate (zero-length) segment\n"); exit(1);
  }

  // ---- translate into box, drop duplicated closing point ----

  int npt = ntot;                       // unique points
  for (int i = 0; i < npt; i++) {
    pts[2*i]   += origin[0];
    pts[2*i+1] += origin[1];
    if (pts[2*i] <= 0.0 || pts[2*i] >= boxhi[0] ||
        pts[2*i+1] <= 0.0 || pts[2*i+1] >= boxhi[1]) {
      fprintf(stderr,"ERROR: point (%g,%g) outside box (0,0)-(%g,%g); "
              "adjust --origin/--chord/--alpha\n",
              pts[2*i],pts[2*i+1],boxhi[0],boxhi[1]);
      exit(1);
    }
  }

  // ---- write data.<prefix> (SPARTA surf format, cf. data.circle) ----

  std::string datafile = "data." + prefix;
  FILE *fp = fopen(datafile.c_str(),"w");
  if (!fp) { fprintf(stderr,"ERROR: cannot write %s\n",datafile.c_str()); exit(1); }

  fprintf(fp,"surf file from bezier_gen: symmetric body, "
             "alpha %.15g %.15g %.15g %.15g chord %g nseg %d\n\n",
          alpha[0],alpha[1],alpha[2],alpha[3],chord,nseg);
  fprintf(fp,"%d points\n%d lines\n\nPoints\n\n",npt,npt);
  for (int i = 0; i < npt; i++)
    fprintf(fp,"%d %.15g %.15g\n",i+1,pts[2*i],pts[2*i+1]);
  fprintf(fp,"\nLines\n\n");
  for (int i = 0; i < npt; i++)
    fprintf(fp,"%d %d %d\n",i+1,i+1,(i+1)%npt+1);
  fclose(fp);

  // ---- write in.<prefix> (modeled on examples/circle/in.circle,
  //      plus drag tally: fx,fy summed over all surf elements) ----

  std::string infile = "in." + prefix;
  fp = fopen(infile.c_str(),"w");
  if (!fp) { fprintf(stderr,"ERROR: cannot write %s\n",infile.c_str()); exit(1); }

  fprintf(fp,
    "################################################################################\n"
    "# 2d flow around a symmetric Bezier body, generated by bezier_gen\n"
    "# alpha (x1,y1,x2,y2): %g %g %g %g\n"
    "# chord %g  nseg/half %d  origin (%g,%g)  seed %d\n"
    "################################################################################\n\n"
    "seed                %d\n"
    "dimension           2\n"
    "global              gridcut 0.0 comm/sort yes\n\n"
    "boundary            o r p\n\n"
    "create_box          0 %g 0 %g -0.5 0.5\n"
    "create_grid         %d %d 1\n"
    "balance_grid        rcb cell\n\n"
    "global              nrho 1.0 fnum 0.001\n\n"
    "species             air.species N O\n"
    "mixture             air N O vstream %g 0 0\n\n"
    "read_surf           data.%s\n"
    "surf_collide        1 diffuse 300.0 0.0\n"
    "surf_modify         all collide 1\n\n"
    "collide             vss air air.vss\n\n"
    "fix                 in emit/face air xlo twopass\n\n"
    "timestep            0.0001\n\n"
    "# drag: per-surf force tallies, time-averaged, then summed\n"
    "compute             forces surf all all fx fy\n"
    "fix                 favg ave/surf all 10 10 100 c_forces[*] ave running\n"
    "compute             drag reduce sum f_favg[1] f_favg[2]\n\n"
    "stats               100\n"
    "stats_style         step cpu np nscoll c_drag[1] c_drag[2]\n"
    "run                 %d\n",
    alpha[0],alpha[1],alpha[2],alpha[3],chord,nseg,origin[0],origin[1],seed,
    seed,boxhi[0],boxhi[1],grid[0],grid[1],vstream,prefix.c_str(),nsteps);
  fclose(fp);

  printf("wrote %s (%d points, %d lines) and %s\n",
         datafile.c_str(),npt,npt,infile.c_str());
  printf("body: nose (%g,%g) tail (%g,%g)\n",
         origin[0],origin[1],origin[0]+chord,origin[1]);
  printf("run with: /path/to/spa_ < %s   (needs air.species, air.vss "
         "in the same directory)\n",infile.c_str());
  return 0;
}
