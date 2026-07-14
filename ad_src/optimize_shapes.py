#!/usr/bin/env python3
"""
optimize_shapes.py -- drag-minimization driver comparing FD and AD
gradients on the symmetric Bezier body (Phase 5 of AD_PLAN.md).

Runs the SAME projected gradient descent as opt_run (update
alpha -= eta * g/|g|, clamped to x1,x2 in [0.2, chord-0.2],
y1,y2 in [0.05, 3.0]); the only difference between the two variants is
the gradient estimator:

  fd : CRN central differences via  opt_run --quiet   (9 x nseeds runs)
  ad : forward-mode pathwise AD via drag_run_ad       (nseeds runs)

NOTE (see AD_PLAN.md, Phase 4): the pathwise AD gradient misses the
capture-area term, so the two variants are EXPECTED to differ -- that
comparison is the point of the final figure.

Usage (from a case directory such as ad_src/bezier):

  # run one optimization cycle, save the trajectory
  python3 ../optimize_shapes.py run --method fd --save fd.json
  python3 ../optimize_shapes.py run --method ad --save ad.json

  # side-by-side final shapes, equal scaling
  python3 ../optimize_shapes.py plot fd.json ad.json -o shapes_ad_vs_fd.png

Common run options:
  --alpha0 x1,y1,x2,y2   start point        (default 1.3,1.0,2.7,0.8)
  --iters N              descent iterations (default 5)
  --eta STEP             step length        (default 0.15)
  --h H                  FD step, fd method (default 0.05)
  --seeds s1,s2,...      CRN seed set       (default 12345)
  --nsettle/--navg N     SPARTA run protocol (default 300/300)
  --chord L              body length        (default 4.0)
  --bindir PATH          where opt_run / drag_run_ad live
                         (default: directory of this script)
"""

import argparse, json, math, os, re, subprocess, sys

# ----------------------------------------------------------------------
# geometry: mirror of BezierGeom::symmetric_body_to_lines (body frame)
# ----------------------------------------------------------------------

def bezier_point(ctrl, t):
    """de Casteljau on control points [(x,y), ...]"""
    px = [c[0] for c in ctrl]
    py = [c[1] for c in ctrl]
    n = len(ctrl)
    for lev in range(1, n):
        for i in range(n - lev):
            px[i] = (1.0 - t) * px[i] + t * px[i + 1]
            py[i] = (1.0 - t) * py[i] + t * py[i + 1]
    return px[0], py[0]

def body_outline(alpha, chord, nseg=100):
    """closed outline (x,y arrays), upper half + mirrored lower half"""
    ctrl = [(0.0, 0.0), (alpha[0], alpha[1]), (alpha[2], alpha[3]),
            (chord, 0.0)]
    upper = [bezier_point(ctrl, i / nseg) for i in range(nseg + 1)]
    lower = [(x, -y) for (x, y) in reversed(upper[:-1])]
    pts = upper + lower
    return [p[0] for p in pts], [p[1] for p in pts]

def body_area(alpha, chord, n=400):
    """cross-sectional area (2D 'volume') of the closed body"""
    xs, ys = body_outline(alpha, chord, n)
    a, m = 0.0, len(xs)
    for i in range(m):
        j = (i + 1) % m
        a += xs[i] * ys[j] - xs[j] * ys[i]
    return abs(a) / 2.0

def area_grad(alpha, chord, h=1e-7):
    """d(area)/d(alpha): area is smooth and analytic, so tiny-h central
    FD is exact to ~1e-9 -- no Monte Carlo involved"""
    g = []
    for j in range(4):
        ap = alpha[:]; am = alpha[:]
        ap[j] += h; am[j] -= h
        g.append((body_area(ap, chord) - body_area(am, chord)) / (2 * h))
    return g

def project_area(alpha, target, chord, tol=1e-10, maxit=8):
    """Newton projection onto area(alpha) = target"""
    for _ in range(maxit):
        err = body_area(alpha, chord) - target
        if abs(err) < tol:
            break
        gA = area_grad(alpha, chord)
        n2 = sum(x * x for x in gA)
        alpha = [alpha[j] - err * gA[j] / n2 for j in range(4)]
    return alpha

# ----------------------------------------------------------------------
# gradient evaluations via the compiled drivers
# ----------------------------------------------------------------------

def run_cmd(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"ERROR: {' '.join(cmd)} failed:\n{r.stdout}\n{r.stderr}")
    return r.stdout

def eval_fd(bindir, alpha, seeds, h, nsettle, navg, chord):
    out = run_cmd([os.path.join(bindir, "opt_run"),
                   "--alpha", ",".join(map(str, alpha)),
                   "--seeds", ",".join(map(str, seeds)),
                   "--h", str(h), "--nsettle", str(nsettle),
                   "--navg", str(navg), "--chord", str(chord), "--quiet"])
    drag = float(re.search(r"drag\s*=\s*([-\d.eE+]+)", out).group(1))
    g = [float(x) for x in
         re.search(r"grad\s*=\s*(.*?)\s*\(h=", out).group(1).split()]
    return drag, g

def eval_ad(bindir, alpha, seeds, nsettle, navg, chord):
    out = run_cmd([os.path.join(bindir, "drag_run_ad"),
                   "--alpha", ",".join(map(str, alpha)),
                   "--seeds", ",".join(map(str, seeds)),
                   "--nsettle", str(nsettle), "--navg", str(navg),
                   "--chord", str(chord)])
    m = re.search(r"SAA mean\s*=\s*([-\d.eE+]+)", out)
    if not m:
        m = re.search(r"drag\(alpha\)\s*=\s*([-\d.eE+]+)", out)
    drag = float(m.group(1))
    g = [float(x) for x in
         re.search(r"AD mean\s*=\s*(.*)", out).group(1).split()]
    return drag, g

# ----------------------------------------------------------------------
# projected gradient descent (identical to opt_main.cpp)
# ----------------------------------------------------------------------

def clamp(a, chord):
    xmin, ymin, ymax = 0.2, 0.05, 3.0
    xmax = chord - 0.2
    a[0] = min(max(a[0], xmin), xmax)
    a[2] = min(max(a[2], xmin), xmax)
    a[1] = min(max(a[1], ymin), ymax)
    a[3] = min(max(a[3], ymin), ymax)
    return a

def optimize(args):
    alpha = [float(x) for x in args.alpha0.split(",")]
    seeds = [int(s) for s in args.seeds.split(",")]
    assert len(alpha) == 4, "--alpha0 needs 4 values"

    traj = {"method": args.method, "eta": args.eta, "h": args.h,
            "seeds": seeds, "chord": args.chord,
            "nsettle": args.nsettle, "navg": args.navg,
            "fixarea": bool(args.fixarea), "lam": args.lam,
            "area0": args.A0 if args.A0 > 0.0
                     else body_area(alpha, args.chord),
            "alpha0": alpha[:], "iterates": []}

    D0 = args.D0 if args.D0 > 0.0 else None

    for it in range(args.iters + 1):
        if args.method == "fd":
            drag, g = eval_fd(args.bindir, alpha, seeds, args.h,
                              args.nsettle, args.navg, args.chord)
        else:
            drag, g = eval_ad(args.bindir, alpha, seeds,
                              args.nsettle, args.navg, args.chord)
        if D0 is None:
            D0 = drag
        A = body_area(alpha, args.chord)

        if args.lam > 0.0:
            # payload-fairing trade-off: J = D/D0 - lam * A/A0
            gA = area_grad(alpha, args.chord)
            g = [g[j] / D0 - args.lam * gA[j] / traj["area0"]
                 for j in range(4)]
            J = drag / D0 - args.lam * A / traj["area0"]
        else:
            J = None

        gn = math.sqrt(sum(x * x for x in g))
        rec = {"iter": it, "alpha": alpha[:], "drag": drag, "gnorm": gn,
               "area": A}
        if J is not None:
            rec["J"] = J
        traj["iterates"].append(rec)
        msg = (f"[{args.method}] iter {it:2d}  drag = {drag:.6e}  "
               f"|g| = {gn:.3e}  area = {A:.4f}")
        if J is not None:
            msg += f"  J = {J:+.4f}"
        print(msg + "  alpha = "
              + " ".join(f"{x:.4f}" for x in alpha), flush=True)
        if it == args.iters or gn == 0.0:
            break
        if args.fixarea:
            # tangential step: remove the area-changing component of g,
            # then Newton-project back onto area = A0 and clamp
            gA = area_grad(alpha, args.chord)
            nA = math.sqrt(sum(x * x for x in gA))
            dot = sum(g[j] * gA[j] for j in range(4)) / (nA * nA)
            gt = [g[j] - dot * gA[j] for j in range(4)]
            gtn = math.sqrt(sum(x * x for x in gt))
            if gtn == 0.0:
                break
            alpha = [alpha[j] - args.eta * gt[j] / gtn for j in range(4)]
            alpha = project_area(alpha, traj["area0"], args.chord)
            alpha = clamp(alpha, args.chord)
            alpha = project_area(alpha, traj["area0"], args.chord)
        else:
            alpha = clamp([alpha[j] - args.eta * g[j] / gn
                           for j in range(4)], args.chord)

    with open(args.save, "w") as f:
        json.dump(traj, f, indent=1)
    print(f"saved trajectory -> {args.save}")

# ----------------------------------------------------------------------
# side-by-side shape plot, equal scaling
# ----------------------------------------------------------------------

def plot(args):
    runs = []
    for path in args.traj:
        with open(path) as f:
            runs.append(json.load(f))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        for r in runs:  # fallback: dump outlines as CSV
            a = r["iterates"][-1]["alpha"]
            xs, ys = body_outline(a, r["chord"])
            fn = f"shape_{r['method']}.csv"
            with open(fn, "w") as f:
                for x, y in zip(xs, ys):
                    f.write(f"{x},{y}\n")
            print(f"matplotlib missing; wrote {fn}")
        return

    fig, axes = plt.subplots(1, len(runs), figsize=(6 * len(runs), 4.4),
                             sharex=True, sharey=True)
    if len(runs) == 1:
        axes = [axes]

    label = {"fd": "FD gradient (CRN, includes capture term)",
             "ad": "AD gradient (pathwise, misses capture term)"}

    for ax, r in zip(axes, runs):
        a0 = r["alpha0"]
        af = r["iterates"][-1]["alpha"]
        d0 = r["iterates"][0]["drag"]
        df = r["iterates"][-1]["drag"]

        xs, ys = body_outline(a0, r["chord"])
        ax.plot(xs + xs[:1], ys + ys[:1], "--", color="0.55", lw=1.4,
                label=f"initial  (drag {d0:.2e})")
        xs, ys = body_outline(af, r["chord"])
        afinal = r["iterates"][-1].get("area")
        lab = f"final    (drag {df:.2e}"
        if afinal is not None:
            lab += f", area {afinal:.3f}"
        lab += ")"
        ax.plot(xs + xs[:1], ys + ys[:1], "-", color="C0", lw=2.2,
                label=lab)
        ax.fill(xs, ys, color="C0", alpha=0.15)

        ax.set_title(label.get(r["method"], r["method"])
                     + f"\nfinal alpha = ("
                     + ", ".join(f"{x:.3f}" for x in af) + ")",
                     fontsize=10)
        ax.set_aspect("equal")          # identical scaling both panels
        ax.grid(alpha=0.3)
        ax.axhline(0.0, color="0.8", lw=0.8, zorder=0)
        ax.legend(fontsize=8, loc="upper right")
        ax.set_xlabel("x (body frame)")
    axes[0].set_ylabel("y")

    fig.suptitle("Drag minimization: FD vs AD gradient "
                 f"({runs[0]['iterates'][-1]['iter']} iterations, "
                 f"same start, step, bounds, seeds)", fontsize=11)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"wrote {args.out}")

def overlay(args):
    runs = []
    for path in args.traj:
        with open(path) as f:
            runs.append(json.load(f))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 4.6))

    # initial shape once (same for all runs)
    r0 = runs[0]
    xs, ys = body_outline(r0["alpha0"], r0["chord"])
    ax.plot(xs + xs[:1], ys + ys[:1], "--", color="0.6", lw=1.4,
            label=f"initial (drag {r0['iterates'][0]['drag']:.2e})")

    colors = {"fd": "C0", "ad": "C3"}
    names  = {"fd": "FD final", "ad": "AD final"}
    for r in runs:
        af = r["iterates"][-1]["alpha"]
        df = r["iterates"][-1]["drag"]
        A  = r["iterates"][-1].get("area")
        xs, ys = body_outline(af, r["chord"])
        c = colors.get(r["method"], None)
        lab = f"{names.get(r['method'], r['method'])} (drag {df:.2e}"
        if A is not None:
            lab += f", area {A:.3f}"
        lab += ")"
        ax.plot(xs + xs[:1], ys + ys[:1], "-", color=c, lw=2.2, label=lab)
        ax.fill(xs, ys, color=c, alpha=0.10)

    ax.set_aspect("equal")
    ax.grid(alpha=0.3)
    ax.axhline(0.0, color="0.85", lw=0.8, zorder=0)
    ax.set_xlabel("x (body frame)")
    ax.set_ylabel("y")
    ax.legend(fontsize=9, loc="upper right")
    ax.set_title("Final shapes overlaid: FD vs AD gradient "
                 "(same start, step, bounds, seeds)", fontsize=11)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"wrote {args.out}")

# ----------------------------------------------------------------------

def main():
    here = os.path.dirname(os.path.abspath(__file__))
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="mode", required=True)

    pr = sub.add_parser("run", help="run one optimization cycle")
    pr.add_argument("--method", choices=["fd", "ad"], required=True)
    pr.add_argument("--alpha0", default="1.3,1.0,2.7,0.8")
    pr.add_argument("--iters", type=int, default=5)
    pr.add_argument("--eta", type=float, default=0.15)
    pr.add_argument("--h", type=float, default=0.05)
    pr.add_argument("--seeds", default="12345")
    pr.add_argument("--nsettle", type=int, default=300)
    pr.add_argument("--navg", type=int, default=300)
    pr.add_argument("--chord", type=float, default=4.0)
    pr.add_argument("--fixarea", action="store_true",
                    help="constrain area(alpha) = area(alpha0)")
    pr.add_argument("--lam", type=float, default=0.0,
                    help="payload trade-off: minimize D/D0 - lam*A/A0")
    pr.add_argument("--D0", type=float, default=0.0,
                    help="fix drag scale (for continuation runs)")
    pr.add_argument("--A0", type=float, default=0.0,
                    help="fix area scale (for continuation runs)")
    pr.add_argument("--bindir", default=here)
    pr.add_argument("--save", required=True)
    pr.set_defaults(func=optimize)

    pp = sub.add_parser("plot", help="side-by-side final shapes")
    pp.add_argument("traj", nargs="+", help="trajectory json files")
    pp.add_argument("-o", "--out", default="shapes_ad_vs_fd.png")
    pp.set_defaults(func=plot)

    po = sub.add_parser("overlay",
                        help="all final shapes on one pair of axes")
    po.add_argument("traj", nargs="+", help="trajectory json files")
    po.add_argument("-o", "--out", default="shapes_overlay.png")
    po.set_defaults(func=overlay)

    args = p.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
