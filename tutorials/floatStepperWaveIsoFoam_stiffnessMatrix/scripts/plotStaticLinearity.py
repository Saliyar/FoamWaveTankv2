#!/usr/bin/env python3
"""Plot the static MoorDyn force against the linear stiffness-matrix prediction.

Usage:
    python3 plotStaticLinearity.py <dir with sweep_*.dat> [--out fig.png]

The sweep files are written by
    moorDynStiffnessMatrix ... -sweep '(dof amplitude n)'
with columns: offset, MoorDyn Fx Fy Fz Mx My Mz, linear Fx Fy Fz Mx My Mz.
The linear prediction is F = preload - C(:, dof) * offset, so its slope is the
(negative) stiffness column that was derived at offset 0.

For each sweep the force component of the same DoF (the diagonal term) and the
coupled component with the largest linear term (largest |C(i, dof)|, i != dof)
are drawn. The error of the diagonal term is printed per offset, relative to
the largest |F - preload| over the sweep (MoorDyn's line relaxation tolerance
threshIC limits the accuracy of very small force changes).
"""
import argparse
import glob
import os
import sys

import numpy as np

DOF = ["x", "y", "z", "rx", "ry", "rz"]
FORCE = ["Fx [N]", "Fy [N]", "Fz [N]", "Mx [N m]", "My [N m]", "Mz [N m]"]
UNIT = ["m", "m", "m", "rad", "rad", "rad"]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("dir")
    p.add_argument("--out", default="staticLinearity.png")
    a = p.parse_args()

    files = sorted(glob.glob(os.path.join(a.dir, "sweep_*.dat")),
                   key=lambda f: DOF.index(os.path.basename(f)[6:-4]))
    if not files:
        print("no sweep_*.dat in", a.dir)
        return 1

    panels = []
    for f in files:
        dof = DOF.index(os.path.basename(f)[6:-4])
        d = np.loadtxt(f)
        off, md, lin = d[:, 0], d[:, 1:7], d[:, 7:13]
        # diagonal term, then the coupled term with the largest |C(i, dof)|
        slope = np.abs(lin[-1] - lin[0])
        coupled = [i for i in range(6) if i != dof]
        ic = max(coupled, key=lambda i: slope[i])
        panels.append((dof, dof, off, md, lin))
        if slope[ic] > 1e-6:
            panels.append((dof, ic, off, md, lin))

        print(f"sweep {DOF[dof]}: relative error of {FORCE[dof].split()[0]}"
              " (linear - MoorDyn) / max|MoorDyn - preload|")
        ref = np.abs(md[:, dof] - md[len(md) // 2, dof]).max()
        for k in range(len(off)):
            if off[k] > 0:
                err = (lin[k, dof] - md[k, dof]) / ref
                print(f"   offset {off[k]:+8.4f} {UNIT[dof]:3s}  {err:+7.2%}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not found, no figure written")
        return 0

    n = len(panels)
    ncol = 2
    nrow = (n + 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(10, 3.2 * nrow), squeeze=False)
    for ax, (dof, i, off, md, lin) in zip(axes.flat, panels):
        ax.plot(off, md[:, i], "-", color="#257de0", lw=2, label="MoorDyn (static)")
        ax.plot(off, lin[:, i], "--", color="#f47230", lw=1.6,
                label="C matrix: preload - C" + f"[{FORCE[i].split()[0]},{DOF[dof]}]" + r"$\cdot\Delta$")
        ax.set_xlabel(f"offset in {DOF[dof]} [{UNIT[dof]}]", color="#333333")
        ax.set_ylabel(FORCE[i], color="#333333")
        ax.grid(True, color="#e5e5e5", lw=0.8)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color("#bbbbbb")
        ax.tick_params(colors="#555555")
        ax.xaxis.set_major_locator(plt.MaxNLocator(5))
        ax.legend(frameon=False, fontsize=8)
    for ax in axes.flat[n:]:
        ax.set_visible(False)
    fig.suptitle("Static mooring force: MoorDyn against the derived stiffness matrix",
                 color="#111111")
    fig.tight_layout()
    fig.savefig(a.out, dpi=130)
    print("Wrote", a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
