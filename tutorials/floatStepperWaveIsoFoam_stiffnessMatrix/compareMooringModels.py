#!/usr/bin/env python3
"""Compare body motion and mooring force of two FloatStepper runs.

Usage:
    python3 compareMooringModels.py <caseA> <caseB> [--labels A B] [--out fig.png]

Each case must have
    postProcessing/floaterMotion/0/floaterStateWriter.dat   (body state)
    log.floatStepperWaveIsoFoam                              (mooring force)

The mooring force is read from the log lines written by moorDynR2
("<t>: force (...), moment (...)") or mooringStiffnessMatrix
("mooringStiffnessMatrix <name> t = <t> xi (...) force (...) moment (...)").
With several outer correctors the last value per time step is kept.

Prints amplitude and mean of surge, heave and pitch over the last two wave
periods and the RMS difference between the cases.
"""
import argparse
import os
import re
import sys

import numpy as np

VEC = r"\(([^)]*)\)"


def read_state(case):
    path = os.path.join(case, "postProcessing", "floaterMotion", "0",
                        "floaterStateWriter.dat")
    t, cor, pitch = [], [], []
    cur = None
    with open(path) as f:
        for line in f:
            if line.startswith("Time ="):
                cur = float(line.split("=")[1])
            elif line.startswith("Centre of rotation"):
                c = [float(v) for v in re.search(VEC, line).group(1).split()]
            elif line.startswith("Orientation"):
                q = [float(v) for v in re.search(VEC, line).group(1).split()]
                t.append(cur)
                cor.append(c)
                # rotation about y: Q = [[c 0 s], [0 1 0], [-s 0 c]]
                pitch.append(np.arctan2(q[2], q[0]))
    t = np.array(t)
    cor = np.array(cor)
    _, idx = np.unique(t[::-1], return_index=True)
    idx = len(t) - 1 - idx
    return t[idx], cor[idx], np.array(pitch)[idx]


def read_force(case):
    path = os.path.join(case, "log.floatStepperWaveIsoFoam")
    re_md = re.compile(r"([-+0-9.eE]+): force " + VEC + r", moment " + VEC)
    re_cm = re.compile(r"mooringStiffnessMatrix \S+ t = (\S+) xi .*force "
                       + VEC + r" moment " + VEC)
    data = {}
    with open(path) as f:
        for line in f:
            m = re_cm.search(line) or re_md.search(line)
            if m:
                tt = float(m.group(1))
                F = [float(v) for v in m.group(2).split()]
                M = [float(v) for v in m.group(3).split()]
                data[tt] = F + M
    t = np.array(sorted(data))
    return t, np.array([data[k] for k in t])


def stats(t, y, t0):
    s = y[t >= t0]
    return 0.5 * (s.max() - s.min()), s.mean()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("caseA")
    p.add_argument("caseB")
    p.add_argument("--labels", nargs=2, default=["MoorDyn", "C matrix"])
    p.add_argument("--period", type=float, default=2.0)
    p.add_argument("--out", default="mooringComparison.png")
    a = p.parse_args()

    runs = []
    for case in (a.caseA, a.caseB):
        ts, cor, pitch = read_state(case)
        tf, F = read_force(case)
        runs.append(dict(ts=ts, cor=cor, pitch=pitch, tf=tf, F=F))

    tend = min(r["ts"][-1] for r in runs)
    t0 = tend - 2 * a.period
    tc = np.linspace(0, tend, 2001)

    rows = [
        ("Surge [m]", "ts", lambda r: r["cor"][:, 0] - r["cor"][0, 0]),
        ("Heave [m]", "ts", lambda r: r["cor"][:, 2] - r["cor"][0, 2]),
        ("Pitch [deg]", "ts", lambda r: np.degrees(r["pitch"])),
        ("Mooring Fx [N]", "tf", lambda r: r["F"][:, 0]),
        ("Mooring Fz [N]", "tf", lambda r: r["F"][:, 2]),
        ("Mooring My [N m]", "tf", lambda r: r["F"][:, 4]),
    ]

    print(f"Statistics over t = {t0:.2f} .. {tend:.2f} s")
    print(f"{'quantity':18s} {'amp ' + a.labels[0]:>14s} {'amp ' + a.labels[1]:>14s}"
          f" {'mean ' + a.labels[0]:>14s} {'mean ' + a.labels[1]:>14s}"
          f" {'RMS diff/amp':>12s}")
    series = []
    for name, tk, fn in rows:
        ys = [fn(r) for r in runs]
        tsr = [r[tk] for r in runs]
        st = [stats(t, y, t0) for t, y in zip(tsr, ys)]
        ia = np.interp(tc, tsr[0], ys[0])
        ib = np.interp(tc, tsr[1], ys[1])
        rms = np.sqrt(np.mean((ia - ib) ** 2)) / max(st[0][0], 1e-30)
        print(f"{name:18s} {st[0][0]:14.5g} {st[1][0]:14.5g} {st[0][1]:14.5g}"
              f" {st[1][1]:14.5g} {rms:12.3%}")
        series.append((name, tsr, ys))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not found, no figure written")
        return 0

    colors = ["#257de0", "#f47230"]
    styles = ["-", "--"]
    fig, axes = plt.subplots(3, 2, figsize=(11, 8), sharex=True)
    for ax, (name, tsr, ys) in zip(axes.T.flat, series):
        for k in range(2):
            ax.plot(tsr[k], ys[k], styles[k], color=colors[k], lw=2 if k == 0 else 1.6,
                    label=a.labels[k])
        ax.set_ylabel(name, color="#333333")
        ax.grid(True, color="#e5e5e5", lw=0.8)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color("#bbbbbb")
        ax.tick_params(colors="#555555")
    for ax in axes[-1]:
        ax.set_xlabel("Time [s]", color="#333333")
    axes[0, 0].legend(frameon=False, loc="upper left")
    fig.suptitle("Moored box in regular waves: MoorDyn vs linear 6x6 stiffness matrix",
                 color="#111111")
    fig.tight_layout()
    fig.savefig(a.out, dpi=130)
    print(f"Wrote {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
