/* ----------------------------------------------------------------------
   svg_shape: write a self-contained SVG overlaying the initial and the
   optimized symmetric-Bezier body outlines. No external dependencies --
   opens in any browser. Geometry comes from BezierGeom.
------------------------------------------------------------------------- */

#ifndef SVG_SHAPE_H
#define SVG_SHAPE_H

#include <string>

// Draw both bodies into one SVG at `path`:
//   - initial shape: faint dashed grey outline + control points
//   - optimized shape: bold colored fill + outline + control points
// alpha_init / alpha_final are length-4 design vectors [x1,y1,x2,y2].
// Returns 0 on success, nonzero if the file could not be written.
int write_shapes_svg(const std::string &path,
                     const double alpha_init[4], double drag_init,
                     const double alpha_final[4], double drag_final,
                     double chord, int nseg);

#endif
