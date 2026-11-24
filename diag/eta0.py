# Still-water level error from a written alpha field: per full-height column
# eta = z_bottom + sum alpha*h (cell heights from the face points).
# python3 eta0.py case time [field]
import sys, numpy as np
exec(open(__file__.replace('eta0.py', 'rI_terms.py')).read().split("if __name__")[0])
C, ts = sys.argv[1], sys.argv[2]; fld = sys.argv[3] if len(sys.argv) > 3 else 'alpha.water'
M = Mesh(C, ts); a, _ = readField('%s/%s/%s' % (C, ts, fld), 1, M.nC, M.bnd)
z0 = np.full(M.nC, 1e9); z1 = np.full(M.nC, -1e9)
for cells, fs in ((M.own, np.arange(len(M.F))), (M.nei, np.arange(M.nI))):
    for f in fs:
        zz = M.p[M.F[f]][:, 2]; c = cells[f]
        z0[c] = min(z0[c], zz.min()); z1[c] = max(z1[c], zz.max())
cols = columns(M); nRow = max(len(c) for c in cols.values())
E = []
for x, cs in cols.items():
    if len(cs) < nRow: continue
    cs = np.asarray(cs); E.append((M.C[cs, 0].mean(), z0[cs].min() + (a[cs]*(z1[cs] - z0[cs])).sum()))
E = np.array(E); i = np.argmax(abs(E[:, 1]))
zc = M.C[:, 2]; lo = zc < zc.min() + 0.05; hi = zc > zc.max() - 0.05
print('sidedness: mean alpha in the bottom 5 cm %.3f (must be 1), top 5 cm %.3f (must be 0)' % (a[lo].mean(), a[hi].mean()))
print('%s t %s %s: %d full columns, max|eta| %.3e m at x %.3f, rms %.3e m' % (C.split('/')[-1], ts, fld, len(E), abs(E[i, 1]), E[i, 0], np.sqrt((E[:, 1]**2).mean())))
