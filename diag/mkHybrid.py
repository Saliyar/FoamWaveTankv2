#!/usr/bin/env python3
"""Hybrid conformal blockMeshDict for the moored-box tank (2D, one cell in y).

x: fine dx 0.0125 for |x| < 0.6, geometric grading (ratio <= 1.1 per cell,
   also across block joints) to the coarse dx 0.05 outside.
z: plain-fine rows from below the box bottom to above the max crest
   (5 x 0.01214 below the box, the fine box rows 10 x 0.01316, 5 x 0.0147
   above the box), graded at <= 1.1 to coarse-size rows, and 3 uniform rows
   at the top whose points are all > outerDistance (0.4 m) from the box.
Same rows across the tank; the box is the hole x +-0.1, z -0.0786..0.053.
"""
import math, sys
R = 1.1
FX, CX = 0.0125, 0.05            # fine and coarse dx
BOXZ0, BOXZ1 = -0.0786, 0.053

def geo(a, q, n):                # n cells a q, a q^2, ...
    return [a*q**k for k in range(1, n + 1)]

def grade_x(L):
    """cells from fine (next to 0.0125) to coarse 0.05 over length L, then the
    rest uniform 0.05: find n, q with q <= R, last*q ~ 0.05"""
    best = None
    for n in range(5, 40):
        q = (CX/FX)**(1.0/(n + 1))
        if q > R: continue
        s = sum(geo(FX, q, n))
        best = (n, q, s); break
    return best

def band(L, a, b, c, r=R):
    """cell sizes over length L between neighbour sizes a (start) and b (end):
    s_k = min(a q^k, b q^(N+1-k), c), the fewest cells N for which some
    q <= r fills L; q by bisection. Every ratio, also to the neighbours, is
    <= q <= r."""
    def prof(q, N):
        return [min(a*q**k, b*q**(N + 1 - k), c) for k in range(1, N + 1)]
    for N in range(1, 400):
        if sum(prof(r, N)) >= L:
            lo, hi = 1.0, r
            if sum(prof(lo, N)) > L:
                continue
            for _ in range(200):
                mid = 0.5*(lo + hi)
                if sum(prof(mid, N)) < L: lo = mid
                else: hi = mid
            cells = prof(hi, N)
            k = L/sum(cells)
            return [x*k for x in cells]
    raise SystemExit('band failed')

# x: -2.5 | coarse | -xc | graded | -0.6 | fine | -0.1 | box | 0.1 | fine | 0.6 | graded | xc | coarse | 2.5
n, q, Lg = grade_x(None)
xc = 0.6 + Lg
nCoarse = round((2.5 - xc)/CX)
xg = geo(FX, q, n)                                   # from x = 0.6 outward
xseg = {
    'coarse': [(2.5 - xc)/nCoarse]*nCoarse,
    'graded': xg,
    'fine': [FX]*round(0.5/FX),
    'box': [FX]*round(0.2/FX),
}
X = [-2.5, -xc, -0.6, -0.1, 0.1, 0.6, xc, 2.5]
xcells = [xseg['coarse'], xg[::-1], xseg['fine'], xseg['box'], xseg['fine'], xg, xseg['coarse']]

# z bands
dzB, dzBox, dzA = 0.0121428571, 0.01316, 0.0147
z_f0 = BOXZ0 - 5*dzB
z_f1 = BOXZ1 + 5*dzA
z_t = 0.455                                           # top 3 rows above this
zbot = band(z_f0 + 0.5, dzB, 1e9, 0.0405)[::-1]      # -0.5 .. z_f0, grows downward to coarse dz
ztop_mid = band(z_t - z_f1, dzA, (0.5 - z_t)/3, 0.049)
Z = [-0.5, z_f0, BOXZ0, BOXZ1, z_f1, z_t, 0.5]
zcells = [zbot, [dzB]*5, [dzBox]*10, [dzA]*5, ztop_mid, [(0.5 - z_t)/3]*3]

def sections(cells):
    """blockMesh multi-grading sections ((relLen n expansion) ...)"""
    out = []; i = 0; N = len(cells)
    while i < N:
        j = i + 1
        if j < N:
            rr = cells[j]/cells[i]
            while j + 1 < N and abs(cells[j + 1]/cells[j] - rr) < 1e-6: j += 1
        seg = cells[i:j + 1]
        out.append((sum(seg), len(seg), seg[-1]/seg[0]))
        i = j + 1
    L = sum(cells)
    return '(' + ' '.join('(%.10g %d %.10g)' % (l/L, k, e) for l, k, e in out) + ')'

def check(name, segs):
    allc = [c for s in segs for c in s]
    rmax = max(max(a/b, b/a) for a, b in zip(allc, allc[1:]))
    return '%s: %d cells, max ratio %.4f, min %.5f max %.5f' % (name, len(allc), rmax, min(allc), max(allc))

nx = len([c for s in xcells for c in s]); nz = len([c for s in zcells for c in s])
nbox = len(xseg['box'])*10
info = [check('x', xcells), check('z', zcells),
        'cells %d (nx %d x nz %d minus box %d)' % (nx*nz - nbox, nx, nz, nbox),
        'x graded band |x| 0.6..%.4f: %d cells q %.4f; coarse dx %.5f x %d' % (xc, n, q, xseg['coarse'][0], nCoarse),
        'z fine band %.4f..%.4f; top 3 rows from z %.3f (dz %.4f)' % (z_f0, z_f1, z_t, (0.5 - z_t)/3)]
for s in info: print('//', s, file=sys.stderr)

Y = [-0.1, 0.1]
def v(i, j, k): return (k*2 + j)*len(X) + i
verts = []
for k, z in enumerate(Z):
    for j, y in enumerate(Y):
        for i, x in enumerate(X):
            verts.append('    (%.10g %.10g %.10g)' % (x, y, z))
blocks = []; inlet = []; outlet = []; bottom = []; atm = []; w1 = []; w2 = []; box = []
for k in range(len(Z) - 1):
    for i in range(len(X) - 1):
        if i == 3 and k == 2:
            continue
        h = [v(i, 0, k), v(i + 1, 0, k), v(i + 1, 1, k), v(i, 1, k),
             v(i, 0, k + 1), v(i + 1, 0, k + 1), v(i + 1, 1, k + 1), v(i, 1, k + 1)]
        blocks.append('    hex (%s) (%d 1 %d) simpleGrading (%s 1 %s)' % (
            ' '.join(map(str, h)), len(xcells[i]), len(zcells[k]), sections(xcells[i]), sections(zcells[k])))
        if i == 0: inlet.append((h[0], h[3], h[7], h[4]))
        if i == len(X) - 2: outlet.append((h[1], h[2], h[6], h[5]))
        if k == 0: bottom.append((h[0], h[1], h[2], h[3]))
        if k == len(Z) - 2: atm.append((h[4], h[5], h[6], h[7]))
        w1.append((h[0], h[1], h[5], h[4])); w2.append((h[3], h[2], h[6], h[7]))
# box faces (hole i=3, k=2), normals out of the fluid into the box
box = [(v(3, 0, 2), v(4, 0, 2), v(4, 1, 2), v(3, 1, 2)),      # under the box: z = BOXZ0
       (v(3, 0, 3), v(3, 1, 3), v(4, 1, 3), v(4, 0, 3)),      # top: z = BOXZ1
       (v(3, 0, 2), v(3, 1, 2), v(3, 1, 3), v(3, 0, 3)),      # x = -0.1
       (v(4, 0, 2), v(4, 0, 3), v(4, 1, 3), v(4, 1, 2))]      # x = +0.1
def patch(name, typ, faces):
    return '    %s\n    {\n        type %s;\n        faces\n        (\n%s\n        );\n    }' % (
        name, typ, '\n'.join('            (%d %d %d %d)' % f for f in faces))
print('''FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// Hybrid conformal mesh, generated by mkHybrid.py
%s
scale 1;
vertices
(
%s
);
blocks
(
%s
);
edges ();
boundary
(
%s
);
mergePatchPairs ();''' % ('\n'.join('// ' + s for s in info), '\n'.join(verts), '\n'.join(blocks), '\n'.join([
    patch('inlet', 'patch', inlet), patch('outlet', 'patch', outlet), patch('bottom', 'patch', bottom),
    patch('wall1', 'empty', w1), patch('wall2', 'empty', w2), patch('box', 'patch', box),
    patch('atmosphere', 'patch', atm)])))
