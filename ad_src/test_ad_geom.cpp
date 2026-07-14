/* ----------------------------------------------------------------------
   test_ad_geom: Phase 4 validation ladder, rungs 1-2 (AD build only).

   Rung 1: d(line normals)/d(alpha)     -- seeded via drag_objective's
           mechanism (read_surf -> seed -> null move_surf), compared to
           central FD of the normals computed by bezier_geom.
   Rung 2: d(cut-cell volumes)/d(alpha) -- AD derivatives on
           grid->cinfo[].volume compared to central FD of volumes from
           two extra SPARTA instances at alpha +/- h (matched by cell id).

   Usage: test_ad_geom [--alpha x1,y1,x2,y2] [--nseg N] [--h H]
   (standalone: needs no gas files; geometry only, no timesteps)
------------------------------------------------------------------------- */

#include "bezier_geom.h"
#include "library.h"
#include "sparta.h"
#include "surf.h"
#include "grid.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <vector>

#ifndef SPARTA_AD
#error "test_ad_geom must be compiled with -DSPARTA_AD"
#endif

static const double CHORD = 4.0;
static const double ORG[2] = {3.0,5.0};

static void cmd(void *spa, const char *str)
{
  char buf[512];
  snprintf(buf,sizeof(buf),"%s",str);
  if (!sparta_command(spa,buf)) {
    fprintf(stderr,"command failed: %s\n",str); exit(1);
  }
}

static void write_surf_file(const char *path, const double *pts, int npt)
{
  FILE *fp = fopen(path,"w");
  fprintf(fp,"test surf\n\n%d points\n%d lines\n\nPoints\n\n",npt,npt);
  for (int i = 0; i < npt; i++)
    fprintf(fp,"%d %.15g %.15g\n",i+1,pts[2*i],pts[2*i+1]);
  fprintf(fp,"\nLines\n\n");
  for (int i = 0; i < npt; i++)
    fprintf(fp,"%d %d %d\n",i+1,i+1,(i+1)%npt+1);
  fclose(fp);
}

// build a SPARTA instance for the given alpha; if seed_derivs, apply
// the Phase 3 seeding + null move_surf
static void *make_instance(const double alpha[4], int nseg, int seed_derivs)
{
  int ntot = 2*nseg;
  std::vector<double> pts(2*(ntot+1)), norms(2*ntot);
  BezierGeom::symmetric_body_to_lines(alpha,CHORD,nseg,
                                      pts.data(),norms.data());
  for (int i = 0; i < ntot; i++) {
    pts[2*i] += ORG[0]; pts[2*i+1] += ORG[1];
  }
  write_surf_file("tmp_test_surf.data",pts.data(),ntot);

  void *spa = NULL;
  char *argv[] = {(char*)"sparta",(char*)"-log",(char*)"none",
                  (char*)"-screen",(char*)"none"};
  sparta_open_no_mpi(5,argv,&spa);

  cmd(spa,"seed 12345");
  cmd(spa,"dimension 2");
  cmd(spa,"global gridcut 0.0 comm/sort yes");
  cmd(spa,"boundary o r p");
  cmd(spa,"create_box 0 10 0 10 -0.5 0.5");
  cmd(spa,"create_grid 20 20 1");
  cmd(spa,"balance_grid rcb cell");
  cmd(spa,"read_surf tmp_test_surf.data");

  if (seed_derivs) {
    SPARTA_NS::SPARTA *sp = (SPARTA_NS::SPARTA *) spa;
    SPARTA_NS::Surf::Line *slines = sp->surf->lines;
    int nsurf = (int) sp->surf->nlocal;
    std::vector<double> jac(2*(2*nseg+1)*4);
    BezierGeom::symmetric_body_jacobian(nseg,jac.data());
    for (int i = 0; i < nsurf; i++) {
      int k1 = (int)(slines[i].id-1), k2 = (k1+1)%ntot;
      for (int j = 0; j < 4; j++) {
        slines[i].p1[0].d[j] = jac[(2*k1)  *4+j];
        slines[i].p1[1].d[j] = jac[(2*k1+1)*4+j];
        slines[i].p2[0].d[j] = jac[(2*k2)  *4+j];
        slines[i].p2[1].d[j] = jac[(2*k2+1)*4+j];
      }
    }
    cmd(spa,"move_surf all trans 0 0 0");
  }
  return spa;
}

// collect per-line normals (values) keyed by line id, and per-cell
// volumes (values) keyed by cell id
static void collect_values(void *spa, std::map<int,double> nrm[2],
                           std::map<long,double> &vol)
{
  SPARTA_NS::SPARTA *sp = (SPARTA_NS::SPARTA *) spa;
  SPARTA_NS::Surf::Line *slines = sp->surf->lines;
  for (int i = 0; i < sp->surf->nlocal; i++) {
    nrm[0][(int)slines[i].id] = slines[i].norm[0].v;
    nrm[1][(int)slines[i].id] = slines[i].norm[1].v;
  }
  SPARTA_NS::Grid *grid = sp->grid;
  for (int icell = 0; icell < grid->nlocal; icell++) {
    if (grid->cells[icell].nsplit != 1) continue;   // skip split cells
    vol[(long)grid->cells[icell].id] = grid->cinfo[icell].volume.v;
  }
}

int main(int narg, char **arg)
{
  double alpha[4] = {1.3,1.0,2.7,0.8};
  int nseg = 25;
  double h = 1.0e-5;

  for (int i = 1; i < narg; i++) {
    if (!strcmp(arg[i],"--alpha") && i+1 < narg) {
      char *tok = strtok(arg[++i],",");
      for (int k = 0; k < 4 && tok; k++) { alpha[k]=atof(tok); tok=strtok(NULL,","); }
    } else if (!strcmp(arg[i],"--nseg") && i+1 < narg) nseg = atoi(arg[++i]);
    else if (!strcmp(arg[i],"--h") && i+1 < narg) h = atof(arg[++i]);
  }

  int ntot = 2*nseg;

  // ---- AD instance at alpha, with seeded derivatives

  void *spa = make_instance(alpha,nseg,1);
  SPARTA_NS::SPARTA *sp = (SPARTA_NS::SPARTA *) spa;

  // ---- rung 1: normals. FD reference from bezier_geom directly

  double maxerr_n = 0.0;
  {
    SPARTA_NS::Surf::Line *slines = sp->surf->lines;
    std::vector<double> np(2*(ntot+1)), nn(2*ntot);
    std::vector<double> pp(2*(ntot+1)), pn(2*ntot);
    for (int j = 0; j < 4; j++) {
      double ap[4],am[4];
      for (int k = 0; k < 4; k++) { ap[k]=alpha[k]; am[k]=alpha[k]; }
      ap[j]+=h; am[j]-=h;
      BezierGeom::symmetric_body_to_lines(ap,CHORD,nseg,np.data(),nn.data());
      BezierGeom::symmetric_body_to_lines(am,CHORD,nseg,pp.data(),pn.data());
      for (int i = 0; i < sp->surf->nlocal; i++) {
        int k1 = (int)(slines[i].id-1);
        for (int cxy = 0; cxy < 2; cxy++) {
          double fd = (nn[2*k1+cxy]-pn[2*k1+cxy])/(2.0*h);
          double ad = slines[i].norm[cxy].d[j];
          maxerr_n = fmax(maxerr_n,fabs(ad-fd));
        }
      }
    }
  }
  printf("rung 1  normals   max |AD - FD| = %.3e\n",maxerr_n);

  // ---- rung 2: cut-cell volumes. FD from two extra instances

  std::map<int,double> nrm0[2]; std::map<long,double> vol0;
  collect_values(spa,nrm0,vol0);

  double maxerr_v = 0.0, maxadv = 0.0;
  for (int j = 0; j < 4; j++) {
    double ap[4],am[4];
    for (int k = 0; k < 4; k++) { ap[k]=alpha[k]; am[k]=alpha[k]; }
    ap[j]+=h; am[j]-=h;
    void *sp_p = make_instance(ap,nseg,0);
    void *sp_m = make_instance(am,nseg,0);
    std::map<int,double> nrmp[2],nrmm[2]; std::map<long,double> volp,volm;
    collect_values(sp_p,nrmp,volp);
    collect_values(sp_m,nrmm,volm);
    sparta_close(sp_p); sparta_close(sp_m);

    SPARTA_NS::Grid *grid = sp->grid;
    for (int icell = 0; icell < grid->nlocal; icell++) {
      if (grid->cells[icell].nsplit != 1) continue;
      long id = (long) grid->cells[icell].id;
      if (!volp.count(id) || !volm.count(id)) continue;  // topology change
      double fd = (volp[id]-volm[id])/(2.0*h);
      double ad = grid->cinfo[icell].volume.d[j];
      maxadv = fmax(maxadv,fabs(ad));
      maxerr_v = fmax(maxerr_v,fabs(ad-fd));
    }
  }
  printf("rung 2  volumes   max |AD - FD| = %.3e   (max |dV/da| = %.3e)\n",
         maxerr_v,maxadv);

  sparta_close(spa);
  remove("tmp_test_surf.data");

  int ok = (maxerr_n < 1e-6) && (maxerr_v < 1e-4*fmax(maxadv,1.0));
  printf("%s\n", ok ? "PASS" : "FAIL");
  return !ok;
}
