# etaHOS - etaI per column next to the box (no solver run), written fields.
# python3 eta_wall.py S6a|B05
import sys, numpy as np
exec(open(__file__.replace('eta_wall.py', 'rI_terms.py')).read().split("if __name__")[0])
case = sys.argv[1]
if case == 'S6a':
    C = '/gpfs/hpc-home/projects/dtu_00075/people/sitali/SWENSE_testing/diag/momentum_rev1/floating/S6a_rev1'
    times = ['14.9978', '18.998', '22.9982']
    sys.path.insert(0, '/home/projects/dtu_00075/people/sitali/SWENSE_testing/postprocessingCodes')
    import plot_cfd_vs_hos as P
    eH = lambda xs, t: P.hos_eta(C + '/constant/modes_HOS_SWENSE_clean.dat', P.DEFAULT_READER, np.array(xs), np.array([t]))[0]
else:
    C = Z + case; times = ['15', '19', '22.9']
    eH = lambda xs, t: np.array([etaHOS(x, t) for x in xs])
for ts in times:
    t = float(ts); M = Mesh(C, ts)
    aI, _ = readField('%s/%s/alphaI' % (C, ts), 1, M.nC, M.bnd)
    zUp = np.full(M.nC, -1e9); zDn = np.full(M.nC, 1e9)
    vert = np.where(abs(M.Sf[:, 2]) > 0.7*M.magSf)[0]
    for side, cells in ((0, M.own), (1, M.nei)):
        fs = vert if side == 0 else vert[vert < M.nI]
        np.maximum.at(zUp, cells[fs], M.Cf[fs, 2]); np.minimum.at(zDn, cells[fs], M.Cf[fs, 2])
    Hc = zUp - zDn
    cols = columns(M); nRow = max(len(c) for c in cols.values())
    cb = [M.Cf[fs] for n, fs, _ in M.bfaces if n == 'box'][0]
    bx0, bx1 = cb[:, 0].min(), cb[:, 0].max()
    rows = []
    for x, cs in cols.items():
        if len(cs) < nRow: continue
        xc = M.C[cs, 0].mean(); dist = max(bx0 - xc, xc - bx1, 0)
        tag = 'box' if dist <= 0.3 else ('ctl' if 0.9 <= abs(xc) <= 1.1 else None)
        if tag is None: continue
        cs = cs[np.argsort(M.C[cs, 2])]
        rows.append((tag, xc, zDn[cs[0]] + (np.clip(aI[cs], 0, 1)*Hc[cs]).sum(), Hc[cs].sum()))
    xs = [r[1] for r in rows]; e = eH(xs, t)
    d = np.array([1e3*(e[i] - r[2]) for i, r in enumerate(rows)]); tg = np.array([r[0] for r in rows])
    for k in ('box', 'ctl'):
        m = tg == k
        if m.any():
            print('%s t %s %s: n %d  max|.| %.1f  rms %.1f  mean %+.1f mm  (depth %.3f..%.3f, box x %.3f..%.3f)' % (case, ts, k, m.sum(), abs(d[m]).max(), np.sqrt((d[m]**2).mean()), d[m].mean(), min(r[3] for r in rows), max(r[3] for r in rows), bx0, bx1))
    print('   per column x: ' + ' '.join('%.2f:%+.0f' % (xs[i], d[i]) for i in np.argsort(xs) if tg[i] == 'box'))
