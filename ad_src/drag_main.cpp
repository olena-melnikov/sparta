/* ----------------------------------------------------------------------
   drag_main: terminal driver for the drag(alpha) objective.

   Usage:
     drag_run [--input file] [--alpha x1,y1,x2,y2] [--seeds s1,s2,...]
              [--nsettle N] [--navg N] [--nseg N] [--chord L] [--verbose]

   --input reads alpha/seeds/case settings from a "key = value" text
   file (see drag.in.example); any CLI flag given in addition to
   --input overrides the corresponding value from the file.

   Run from a directory containing air.species and air.vss
   (e.g. bezier/, now nested under ad_src/).

   Prints one drag realization per seed, then mean and std error --
   the SAA estimate of drag(alpha).
------------------------------------------------------------------------- */

#include "drag_objective.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

static const char *USAGE =
    "usage: drag_run [--input file] [--alpha x1,y1,x2,y2] "
    "[--seeds s1,s2,...] [--nsettle N] [--navg N] "
    "[--nseg N] [--chord L] [--verbose]\n";

/* ----------------------------------------------------------------------
   parse "key = value" input file; recognizes the same keys as the CLI
   flags (alpha, seeds, nsettle, navg, nseg, chord, verbose). Unknown
   keys and blank/'#' comment lines are ignored.
------------------------------------------------------------------------- */

static void parse_input_file(const char *path, std::vector<double> &alpha,
                              std::vector<int> &seeds, DragCase &c)
{
  FILE *fp = fopen(path,"r");
  if (!fp) {
    fprintf(stderr,"ERROR: cannot open input file %s\n",path);
    exit(1);
  }

  char line[1024];
  while (fgets(line,sizeof(line),fp)) {
    char *hash = strchr(line,'#');
    if (hash) *hash = '\0';

    char *eq = strchr(line,'=');
    if (!eq) continue;
    *eq = '\0';
    char *key = strtok(line," \t\n");
    char *val = eq+1;
    if (!key) continue;

    if (!strcmp(key,"alpha")) {
      alpha.clear();
      char *tok = strtok(val,", \t\n");
      while (tok) { alpha.push_back(atof(tok)); tok = strtok(NULL,", \t\n"); }
    } else if (!strcmp(key,"seeds")) {
      seeds.clear();
      char *tok = strtok(val,", \t\n");
      while (tok) { seeds.push_back(atoi(tok)); tok = strtok(NULL,", \t\n"); }
    } else if (!strcmp(key,"nsettle")) {
      c.nsettle = atoi(val);
    } else if (!strcmp(key,"navg")) {
      c.navg = atoi(val);
    } else if (!strcmp(key,"nseg")) {
      c.nseg = atoi(val);
    } else if (!strcmp(key,"chord")) {
      c.chord = atof(val);
    } else if (!strcmp(key,"verbose")) {
      char *tok = strtok(val," \t\n");
      c.verbose = tok ? atoi(tok) : 1;
    }
  }
  fclose(fp);
}

int main(int narg, char **arg)
{
  std::vector<double> alpha;
  std::vector<int> seeds;
  DragCase c;

  // first pass: --input, so later CLI flags on the command line can
  // still override values loaded from the file
  for (int i = 1; i < narg; i++) {
    if (!strcmp(arg[i],"--input") && i+1 < narg) {
      parse_input_file(arg[i+1],alpha,seeds,c);
      break;
    }
  }

  for (int i = 1; i < narg; i++) {
    if (!strcmp(arg[i],"--input") && i+1 < narg) {
      i++;  // already handled above
    } else if (!strcmp(arg[i],"--alpha") && i+1 < narg) {
      alpha.clear();
      char *tok = strtok(arg[++i],",");
      while (tok) { alpha.push_back(atof(tok)); tok = strtok(NULL,","); }
    } else if (!strcmp(arg[i],"--seeds") && i+1 < narg) {
      seeds.clear();
      char *tok = strtok(arg[++i],",");
      while (tok) { seeds.push_back(atoi(tok)); tok = strtok(NULL,","); }
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
    } else {
      fprintf(stderr,"%s",USAGE);
      return 1;
    }
  }

  if (alpha.size() != 4) {
    fprintf(stderr,"ERROR: need exactly 4 alpha values: x1,y1,x2,y2 "
            "(via --alpha or --input)\n");
    return 1;
  }
  if (seeds.empty()) seeds.push_back(12345);

  int n = (int) seeds.size();
  std::vector<double> vals(n);

  int ad = drag_has_ad();
  double gsum[4] = {0,0,0,0};

  // n is the number of seeds
  // if ad is not turned on, just a zero gradient will bre returned
  for (int k = 0; k < n; k++) {
    double fy, g[4];
    vals[k] = drag(alpha.data(),seeds[k],c,&fy,ad ? g : 0);
    printf("seed %-8d drag = %.8e   (fy = %.2e)\n",seeds[k],vals[k],fy);
    if (ad) {
      printf("              grad = %.6e %.6e %.6e %.6e\n",
             g[0],g[1],g[2],g[3]);
      for (int j = 0; j < 4; j++) gsum[j] += g[j];
    }
    else {
      printf("AD not available in this build; no gradient computed\n");
    }
  }

  double mean = 0.0;
  for (int k = 0; k < n; k++) mean += vals[k];
  mean /= n;

  if (n > 1) {
    double var = 0.0;
    for (int k = 0; k < n; k++) var += (vals[k]-mean)*(vals[k]-mean);
    var /= (n-1);
    printf("drag(alpha) SAA mean = %.8e   stderr = %.2e   (n=%d)\n",
           mean,sqrt(var/n),n);
  } else {
    printf("drag(alpha) = %.8e   (single realization)\n",mean);
  }

  if (ad)
    printf("d(drag)/d(alpha) AD mean = %.6e %.6e %.6e %.6e\n",
           gsum[0]/n,gsum[1]/n,gsum[2]/n,gsum[3]/n);
  return 0;
}
