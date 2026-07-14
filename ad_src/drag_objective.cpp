/* ----------------------------------------------------------------------
   drag_objective implementation. See drag_objective.h.

   Uses only SPARTA's public library interface (sparta_open_no_mpi,
   sparta_command, sparta_extract_compute) -- no SPARTA sources touched.

   Geometry currently reaches SPARTA through a temporary surf data
   file + read_surf. This is the Path A baseline; the derivative chain
   is glued at this seam with the analytic Jacobian from bezier_geom
   (symmetric_body_jacobian). Path B replaces write_surf_tmp/read_surf
   with an in-memory create_surf command; nothing else here changes.
------------------------------------------------------------------------- */

#include "drag_objective.h"
#include "bezier_geom.h"
#include "library.h"

#ifdef SPARTA_AD
// AD build only: direct access to Surf::lines for derivative seeding.
// These headers compile with sfloat == forward dual (same -DSPARTA_AD).
#include "sparta.h"
#include "surf.h"
#endif

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <unistd.h>

/* ---------------------------------------------------------------------- */

static void die(const char *msg)
{
  fprintf(stderr,"drag_objective ERROR: %s\n",msg);
  exit(1);
}

static void cmd(void *spa, const char *str)
{
  // sparta_command takes a mutable char*; copy to be safe
  char buf[1024];
  snprintf(buf,sizeof(buf),"%s",str);
  if (!sparta_command(spa,buf)) die(buf);
}

/* ----------------------------------------------------------------------
   write clockwise closed loop (npt unique points, already translated
   into the box) as a SPARTA surf data file
------------------------------------------------------------------------- */

static void write_surf_tmp(const char *path, const double *pts, int npt)
{
  FILE *fp = fopen(path,"w");
  if (!fp) die("cannot write temp surf file");
  fprintf(fp,"temp surf file from drag_objective\n\n");
  fprintf(fp,"%d points\n%d lines\n\nPoints\n\n",npt,npt);
  for (int i = 0; i < npt; i++)
    fprintf(fp,"%d %.15g %.15g\n",i+1,pts[2*i],pts[2*i+1]);
  fprintf(fp,"\nLines\n\n");
  for (int i = 0; i < npt; i++)
    fprintf(fp,"%d %d %d\n",i+1,i+1,(i+1)%npt+1);
  fclose(fp);
}

/* ---------------------------------------------------------------------- */

double drag(const double alpha[4], int seed, const DragCase &c,
            double *fyout, double *gradout)
{
  clock_t t0 = clock();

  // fail loudly up front: SPARTA's own error exit is silent when its
  // screen output is suppressed
  if (access(c.species_file,R_OK) || access(c.vss_file,R_OK)) {
    fprintf(stderr,"drag_objective ERROR: cannot find %s / %s in the "
            "current directory.\nRun from a case directory such as "
            "ad_src/bezier (or set DragCase::species_file/vss_file).\n",
            c.species_file,c.vss_file);
    exit(1);
  }
  // ---- geometry: alpha -> closed clockwise loop, translated into box

  int ntot = 2*c.nseg;
  double *pts   = new double[2*(ntot+1)];
  double *norms = new double[2*ntot];
  BezierGeom::symmetric_body_to_lines(alpha,c.chord,c.nseg,pts,norms);

  if (BezierGeom::signed_area(ntot,pts) >= 0.0)
    die("body not clockwise; need y1,y2 > 0");
  if (BezierGeom::min_segment_length(ntot,pts) <= 0.0)
    die("degenerate segment");

  for (int i = 0; i < ntot; i++) {           // translate, drop closing pt
    pts[2*i]   += c.origin[0];
    pts[2*i+1] += c.origin[1];
    if (pts[2*i] <= 0.0 || pts[2*i] >= c.boxhi[0] ||
        pts[2*i+1] <= 0.0 || pts[2*i+1] >= c.boxhi[1])
      die("body outside box; adjust origin/chord/alpha");
  }

  char surfpath[256];
  snprintf(surfpath,sizeof(surfpath),"tmp_surf_%d_%d.data",
           (int) getpid(),seed);
  write_surf_tmp(surfpath,pts,ntot);
  delete [] pts;
  delete [] norms;

  // ---- open SPARTA in-process

  void *spa = NULL;
  if (c.verbose) {
    char *argv[] = {(char*)"sparta",(char*)"-log",(char*)"none"};
    sparta_open_no_mpi(3,argv,&spa);
  } else {
    char *argv[] = {(char*)"sparta",(char*)"-log",(char*)"none",
                    (char*)"-screen",(char*)"none"};
    sparta_open_no_mpi(5,argv,&spa);
  }
  if (!spa) die("sparta_open failed");

  // ---- case setup (mirrors examples/circle)

  char line[1024];
  snprintf(line,sizeof(line),"seed %d",seed);                     cmd(spa,line);
  cmd(spa,"dimension 2");
  cmd(spa,"global gridcut 0.0 comm/sort yes");
  cmd(spa,"boundary o r p");
  snprintf(line,sizeof(line),"create_box 0 %.15g 0 %.15g -0.5 0.5",
           c.boxhi[0],c.boxhi[1]);                                cmd(spa,line);
  snprintf(line,sizeof(line),"create_grid %d %d 1",
           c.grid[0],c.grid[1]);                                  cmd(spa,line);
  cmd(spa,"balance_grid rcb cell");
  snprintf(line,sizeof(line),"global nrho %.15g fnum %.15g",
           c.nrho,c.fnum);                                        cmd(spa,line);
  snprintf(line,sizeof(line),"species %s N O",c.species_file);    cmd(spa,line);
  snprintf(line,sizeof(line),"mixture air N O vstream %.15g 0 0",
           c.vstream);                                            cmd(spa,line);
  snprintf(line,sizeof(line),"read_surf %s",surfpath);            cmd(spa,line);

#ifdef SPARTA_AD
  // ---- seed d(surf points)/d(alpha) from the analytic Bezier Jacobian.
  // read_surf already computed normals and cut-cell volumes from
  // UNSEEDED points, so afterwards we issue a null move_surf: it re-runs
  // compute_line_normal and the full grid->surf2grid remap, letting the
  // seeded derivatives flow into normals and cut-cell volumes.
  {
    SPARTA_NS::SPARTA *sp = (SPARTA_NS::SPARTA *) spa;
    SPARTA_NS::Surf::Line *slines = sp->surf->lines;
    int nsurf = (int) sp->surf->nlocal;
    if (nsurf != 2*c.nseg) die("AD seeding: unexpected surf count");

    int npt = 2*c.nseg;                     // unique loop points
    int nrow = 2*(2*c.nseg+1);
    double *jac = new double[nrow*4];
    BezierGeom::symmetric_body_jacobian(c.nseg,jac);

    // translation by c.origin is alpha-independent: Jacobian unchanged.
    // line with id k+1 connects loop points k and (k+1)%npt (the order
    // write_surf_tmp emitted them)
    for (int i = 0; i < nsurf; i++) {
      int k1 = (int) (slines[i].id - 1);
      int k2 = (k1+1) % npt;
      for (int j = 0; j < 4; j++) {
        slines[i].p1[0].d[j] = jac[(2*k1)  *4 + j];
        slines[i].p1[1].d[j] = jac[(2*k1+1)*4 + j];
        slines[i].p2[0].d[j] = jac[(2*k2)  *4 + j];
        slines[i].p2[1].d[j] = jac[(2*k2+1)*4 + j];
      }
    }
    delete [] jac;
  }
  cmd(spa,"move_surf all trans 0 0 0");     // null move: propagate seeds
#endif

  if (c.specular) cmd(spa,"surf_collide 1 specular");
  else cmd(spa,"surf_collide 1 diffuse 300.0 0.0");
  cmd(spa,"surf_modify all collide 1");
  if (c.collisions) {
    snprintf(line,sizeof(line),"collide vss air %s",c.vss_file);  cmd(spa,line);
  }
  cmd(spa,"fix in emit/face air xlo twopass");
  snprintf(line,sizeof(line),"timestep %.15g",c.tstep);           cmd(spa,line);
  if (c.verbose) cmd(spa,"stats 100");

  // ---- settle (flow development, excluded from the average)

  snprintf(line,sizeof(line),"run %d",c.nsettle);                 cmd(spa,line);

  // ---- averaging window: per-surf fx,fy sampled every step,
  //      averaged over navg steps, summed over all surf elements

  cmd(spa,"compute forces surf all all fx fy");
  snprintf(line,sizeof(line),
           "fix favg ave/surf all 1 %d %d c_forces[*]",
           c.navg,c.navg);                                        cmd(spa,line);
  cmd(spa,"compute drag reduce sum f_favg[1] f_favg[2]");
  snprintf(line,sizeof(line),"run %d",c.navg);                    cmd(spa,line);

  // ---- extract through memory (style 0 = global, type 1 = vector)

#ifdef SPARTA_AD
  sfloat *vec = (sfloat *) sparta_extract_compute(spa,(char*)"drag",0,1);
  if (!vec) die("extract of compute drag failed");
  double fx = vec[0].v;
  if (fyout) *fyout = vec[1].v;
  if (gradout)
    for (int j = 0; j < 4; j++) gradout[j] = vec[0].d[j];
#else
  double *vec = (double *) sparta_extract_compute(spa,(char*)"drag",0,1);
  if (!vec) die("extract of compute drag failed");
  double fx = vec[0];
  if (fyout) *fyout = vec[1];
  if (gradout) {
    static int warned = 0;
    if (!warned) {
      fprintf(stderr,"drag_objective warning: gradient requested from "
              "the stock build; returning NaN (build with -DSPARTA_AD / "
              "make drag_run_ad)\n");
      warned = 1;
    }
    for (int j = 0; j < 4; j++) gradout[j] = 0.0/0.0;
  }
#endif

  sparta_close(spa);
  unlink(surfpath);

  if (c.progress)
    fprintf(stderr,"  [run] alpha=(%g,%g,%g,%g) seed=%d  drag=%.6e  (%.1fs)\n",
            alpha[0],alpha[1],alpha[2],alpha[3],seed,fx,
            (double)(clock()-t0)/CLOCKS_PER_SEC);
  return fx;
}

/* ---------------------------------------------------------------------- */

int drag_has_ad()
{
#ifdef SPARTA_AD
  return 1;
#else
  return 0;
#endif
}

double drag_avg(const double alpha[4], const int *seeds, int nseeds,
                const DragCase &c, double *out, double *grad)
{
  double sum = 0.0, g[4], gsum[4] = {0.0,0.0,0.0,0.0};
  for (int k = 0; k < nseeds; k++) {
    double d = drag(alpha,seeds[k],c,0,grad ? g : 0);
    if (out) out[k] = d;
    if (grad) for (int j = 0; j < 4; j++) gsum[j] += g[j];
    sum += d;
  }
  if (grad) for (int j = 0; j < 4; j++) grad[j] = gsum[j]/nseeds;
  return sum/nseeds;
}

/* ----------------------------------------------------------------------
   CRN central-difference gradient; see header
------------------------------------------------------------------------- */

double grad_fd(const double alpha[4], const int *seeds, int nseeds,
               const DragCase &c, double h, double grad[4])
{
  double ap[4],am[4];

  if (c.progress)
    fprintf(stderr,"grad_fd: central differences, %d runs total\n",
            (2*4+1)*nseeds);

  for (int j = 0; j < 4; j++) {
    if (c.progress)
      fprintf(stderr," component %d/4 (alpha[%d] +/- %g):\n",j+1,j,h);
    for (int i = 0; i < 4; i++) { ap[i] = alpha[i]; am[i] = alpha[i]; }
    ap[j] += h;
    am[j] -= h;
    double dp = drag_avg(ap,seeds,nseeds,c);
    double dm = drag_avg(am,seeds,nseeds,c);
    grad[j] = (dp-dm)/(2.0*h);
  }

  return drag_avg(alpha,seeds,nseeds,c);
}
