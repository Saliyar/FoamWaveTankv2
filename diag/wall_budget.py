# B05: alphaI volume change per column next to the box over one period,
# split into the box-wall face flux and the internal faces (= total - wall).
# Wall flux rate per cell: alphaI_cell (zeroGradient on box) * phiIRel_b,
# phiIRel_b = UInc_b . Sf - meshPhi_b (the projection does not change the box
# flux: psiProj zeroGradient there). Trapezoid over the 0.1 s writes.
import sys, numpy as np
exec(open(__file__.replace('wall_budget.py', 'rI_terms.py')).read().split("if __name__")[0])
CASE = sys.argv[1] if len(sys.argv) > 1 else 'B05'
C = Z + CASE
def sample(ts):
    M = Mesh(C, ts)
    aI, aIb = readField('%s/%s/alphaI' % (C, ts), 1, M.nC, M.bnd)
    UI, UIb = readField('%s/%s/UInc' % (C, ts), 3, M.nC, M.bnd)
    name, fs, cells = [b for b in M.bfaces if b[0] == 'box'][0]
    s = open('%s/%s/meshPhi' % (C, ts)).read(); k = s.index('box', s.index('boundaryField')); blk = s[k:s.index('}', k)]
    m = re.search(r'value\s+nonuniform\s+List<scalar>\s*\d+\s*\(([^)]*)\)', blk)
    mphi = np.fromstring(m.group(1), sep=' ') if m else np.zeros(len(fs))
    phiRel = np.einsum('ij,ij->i', UIb['box'], M.Sf[fs]) - mphi      # >0: out of the fluid into the wall
    wall = np.zeros(M.nC); ab = aIb.get('box') if isinstance(aIb, dict) else None
    ab = np.asarray(ab, float) if ab is not None and np.size(ab) == len(fs) else aI[cells]   # patch value (cell value for zeroGradient)
    np.add.at(wall, cells, -ab*phiRel)   # alphaI volume rate into the cell
    return M, aI, wall, cells
def budget(t0, t1):
    ts = ['%g' % round(t, 1) for t in np.arange(t0, t1 + 1e-9, 0.1)]
    M0, a0, _, cells = sample(ts[0]); cols = columns(M0)
    sel = {x: cs for x, cs in cols.items() if np.intersect1d(cs, cells).size}
    W = {x: 0.0 for x in sel}; prev = None
    for k, t in enumerate(ts):
        M, a, wall, _ = sample(t)
        r = {x: wall[cs].sum() for x, cs in sel.items()}
        if prev is not None:
            for x in sel: W[x] += 0.5*(r[x] + prev[x])*0.1
        prev = r
        if k == 0: V0 = {x: (a[cs]*M.V[cs]).sum() for x, cs in sel.items()}
        V1 = {x: (a[cs]*M.V[cs]).sum() for x, cs in sel.items()}
    A = 0.05*0.2
    print(CASE + ' %s-%s s (one period ~ %.3f s); per column: alphaI volume change as level (mm over 0.05 x 0.2 m)' % (ts[0], ts[-1], float(ts[-1]) - float(ts[0])))
    print('  %8s %6s %9s %9s %9s' % ('x0 (m)', 'where', 'total', 'wall', 'internal'))
    for x in sorted(sel):
        cs = sel[x]; where = 'side' if len(cs) == max(len(c) for c in cols.values()) else 'over/under'
        tot = V1[x] - V0[x]; print('  %8.3f %6s %9.2f %9.2f %9.2f' % (x, where, 1e3*tot/A, 1e3*W[x]/A, 1e3*(tot - W[x])/A))
for a, b in ((15.0, 17.0), (19.0, 21.0)):
    budget(a, b)
