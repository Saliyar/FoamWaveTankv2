# etaHOS - etaI at sample x next to the box: etaI(x) = zbottom + sum alphaI*h
# over the cells whose x-extent contains x (works on split-hex meshes).
# python3 eta_wall_x.py S6a|B05
import sys, numpy as np
exec(open(__file__.replace('eta_wall_x.py', 'rI_terms.py')).read().split("if __name__")[0])
case = sys.argv[1]
if case == 'S6a':
    C = '/gpfs/hpc-home/projects/dtu_00075/people/sitali/SWENSE_testing/diag/momentum_rev1/floating/S6a_rev1'
    times = ['2.0001', '8.0004', '14.9978', '18.998', '22.9982']
    sys.path.insert(0, '/home/projects/dtu_00075/people/sitali/SWENSE_testing/postprocessingCodes')
    import plot_cfd_vs_hos as P
    eH = lambda xs, t: P.hos_eta(C + '/constant/modes_HOS_SWENSE_clean.dat', P.DEFAULT_READER, np.array(xs), np.array([t]))[0]
else:
    C = Z + case; times = ['15', '19', '22.9']
    eH = lambda xs, t: np.array([etaHOS(x, t) for x in xs])
for ts in times:
    t = float(ts); M = Mesh(C, ts)
    aI, _ = readField('%s/%s/alphaI' % (C, ts), 1, M.nC, M.bnd)
    fx0 = np.array([M.p[fv][:, 0].min() for fv in M.F]); fx1 = np.array([M.p[fv][:, 0].max() for fv in M.F])
    fz0 = np.array([M.p[fv][:, 2].min() for fv in M.F]); fz1 = np.array([M.p[fv][:, 2].max() for fv in M.F])
    cx0 = np.full(M.nC, 1e9); cx1 = np.full(M.nC, -1e9); cz0 = np.full(M.nC, 1e9); cz1 = np.full(M.nC, -1e9)
    for cells, fs in ((M.own, np.arange(len(M.F))), (M.nei, np.arange(M.nI))):
        np.minimum.at(cx0, cells[fs], fx0[fs]); np.maximum.at(cx1, cells[fs], fx1[fs])
        np.minimum.at(cz0, cells[fs], fz0[fs]); np.maximum.at(cz1, cells[fs], fz1[fs])
    hc = M.V/((cx1 - cx0)*(M.p[:, 1].max() - M.p[:, 1].min()))   # height of a cell (2D, one cell in y)
    cb = [M.Cf[fs] for n, fs, _ in M.bfaces if n == 'box'][0]; bx0, bx1 = cb[:, 0].min(), cb[:, 0].max()
    zb0 = M.p[:, 2].min()
    dist = np.arange(0.02, 0.301, 0.04) + 0.0013
    xs = list(bx0 - dist[::-1]) + list(bx1 + dist) + [-1.0013, 1.0013]
    eI = []
    for x in xs:
        m = (cx0 < x) & (cx1 > x)
        eI.append(zb0 + (np.clip(aI[m], 0, 1)*hc[m]).sum())
    d = 1e3*(eH(xs, t) - np.array(eI))
    nb = 2*len(dist)
    print('%s t %s box x %.3f..%.3f | next to box (0.02-0.30 m): max|.| %.1f rms %.1f mean %+.1f mm | upstream side mean %+.1f, downstream side mean %+.1f | control x -1, +1: %+.1f %+.1f mm'
          % (case, ts, bx0, bx1, abs(d[:nb]).max(), np.sqrt((d[:nb]**2).mean()), d[:nb].mean(), d[:len(dist)].mean(), d[len(dist):nb].mean(), d[-2], d[-1]))
