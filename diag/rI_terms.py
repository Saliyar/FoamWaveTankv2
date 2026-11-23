# Post-processing of written fields (no solver run): terms of the incident
# residual rI in band cells, column level error etaHOS - etaI, and the
# correlation of band G_f with rho_w g (etaHOS - etaI)/dz.
# Usage: python3 rI_terms.py  (cases ZD and B05, times 15 19 22.9)
import numpy as np, re, os, sys
Z = os.path.dirname(os.path.abspath(__file__)) + '/'
RW, RA, G = 1000.0, 1.0, 9.81
NUW, NUA = 1e-6, 1.48e-5
H = .1195; T = 1.993; h = .5; w = 2*np.pi/T; L = L0 = G*T*T/(2*np.pi)
for i in range(100): L = L0*np.tanh(2*np.pi*h/L)
K = 2*np.pi/L; s = np.tanh(K*h); a2 = K*H*H/4*(3 - s*s)/(4*s**3)
def etaHOS(x, t, r=12.):
    th = K*x - w*t + 1.5*np.pi; f = np.sin(np.pi/2*min(t, r)/r)
    return f*(H/2*np.cos(th) + a2*np.cos(2*th))

# ---------- reading ----------
def _list(s, start):
    m = re.compile(r'(\d+)\s*\(').search(s, start); n = int(m.group(1))
    b = s.index('\n)', m.end()) if n > 0 else m.end()
    return s[m.end():b], n
def readList(f, nc):
    s = open(f).read(); i = s.find('FoamFile'); i = s.index('}', i)
    body, n = _list(s, i)
    v = np.fromstring(body.replace('(', ' ').replace(')', ' '), sep=' ')
    return v.reshape(n, nc) if nc > 1 else v
def readFaces(f):
    s = open(f).read(); i = s.index('}', s.find('FoamFile'))
    body, n = _list(s, i)
    return [list(map(int, m.split())) for m in re.findall(r'\d+\(([^)]*)\)', body)]
def readBoundary(f):
    s = open(f).read(); out = []
    for m in re.finditer(r'(\w+)\s*\{[^}]*?type\s+(\w+);[^}]*?nFaces\s+(\d+);[^}]*?startFace\s+(\d+);', s):
        out.append((m.group(1), m.group(2), int(m.group(3)), int(m.group(4))))
    return out
def readField(f, nc, nCells, patches):
    s = open(f).read(); i = s.index('internalField')
    line = s[i:s.index('\n', i)]
    if 'nonuniform' in line:
        a = s.index('(', s.index('List', i)); b = s.index('\n)', a)
        v = np.fromstring(s[a+1:b].replace('(', ' ').replace(')', ' '), sep=' ')
        iv = v.reshape(-1, nc) if nc > 1 else v
    else:
        u = np.fromstring(line.split('uniform')[1].replace('(', ' ').replace(')', ' ').replace(';', ' '), sep=' ')
        iv = np.tile(u, (nCells, 1)) if nc > 1 else np.full(nCells, u[0])
    bf = {}
    j = s.index('boundaryField', i)
    for name, typ, n, st in patches:
        if typ == 'empty' or n == 0: continue
        k = s.find(name, j); k2 = s.index('}', k); blk = s[k:k2]
        m = re.search(r'value\s+(nonuniform|uniform)', blk)
        if m is None: bf[name] = None; continue
        if m.group(1) == 'uniform':
            u = np.fromstring(blk[m.end():blk.index(';', m.end())].replace('(', ' ').replace(')', ' '), sep=' ')
            bf[name] = np.tile(u, (n, 1)) if nc > 1 else np.full(n, u[0])
        else:
            a = s.index('(', k + m.end()); b = s.index('\n)', a)
            v = np.fromstring(s[a+1:b].replace('(', ' ').replace(')', ' '), sep=' ')
            bf[name] = v.reshape(-1, nc) if nc > 1 else v
    return iv, bf

class Mesh:
    def __init__(s, case, t):
        P = case + '/constant/polyMesh/'
        s.F = readFaces(P + 'faces'); s.own = readList(P + 'owner', 1).astype(int)
        s.nei = readList(P + 'neighbour', 1).astype(int); s.bnd = readBoundary(P + 'boundary')
        s.nI = len(s.nei); s.nC = s.own.max() + 1
        s.p0 = readList(P + 'points', 3)
        pt = '%s/%s/polyMesh/points' % (case, t)
        s.p = readList(pt, 3) if os.path.exists(pt) else s.p0
        s.geom()
    def geom(s):
        p = s.p; nf = len(s.F)
        Sf = np.zeros((nf, 3)); Cf = np.zeros((nf, 3))
        for f, fv in enumerate(s.F):
            q = p[fv]; c0 = q.mean(0); A = np.zeros(3); C = np.zeros(3); At = 0
            for k in range(len(q)):
                a = np.cross(q[k] - c0, q[(k+1) % len(q)] - c0)/2; m = np.linalg.norm(a)
                A += a; C += m*(q[k] + q[(k+1) % len(q)] + c0)/3; At += m
            Sf[f] = A; Cf[f] = C/At
        s.Sf, s.Cf = Sf, Cf; s.magSf = np.linalg.norm(Sf, axis=1)
        cest = np.zeros((s.nC, 3)); cnt = np.zeros(s.nC)
        np.add.at(cest, s.own, Cf); np.add.at(cnt, s.own, 1)
        np.add.at(cest, s.nei, Cf[:s.nI]); np.add.at(cnt, s.nei, 1); cest /= cnt[:, None]
        V = np.zeros(s.nC); C = np.zeros((s.nC, 3))
        for side, cells, fs in ((1, s.own, np.arange(nf)), (-1, s.nei, np.arange(s.nI))):
            pv = side*np.einsum('ij,ij->i', Sf[fs], Cf[fs] - cest[cells])/3
            np.add.at(V, cells, pv); np.add.at(C, cells, pv[:, None]*(0.75*Cf[fs] + 0.25*cest[cells]))
        s.V, s.C = V, C/V[:, None]
        o, n = s.own[:s.nI], s.nei
        d = s.C[n] - s.C[o]
        s.wl = np.einsum('ij,ij->i', Sf[:s.nI], s.C[n] - Cf[:s.nI])/np.einsum('ij,ij->i', Sf[:s.nI], d)
        s.dI = d
        s.bfaces = []   # (patch name, face ids, owner cells) of non-empty patches
        for name, typ, nn, st in s.bnd:
            if typ == 'empty' or nn == 0: continue
            fs = np.arange(st, st + nn); s.bfaces.append((name, fs, s.own[fs]))
        # 2D reconstruct tensor (x, z) over internal and non-empty boundary faces
        M = np.zeros((s.nC, 2, 2))
        def addM(fs, cells):
            n2 = Sf[fs][:, [0, 2]]; np.add.at(M, cells, np.einsum('ij,ik->ijk', n2, n2)/s.magSf[fs][:, None, None])
        addM(np.arange(s.nI), o); addM(np.arange(s.nI), n)
        for _, fs, cells in s.bfaces: addM(fs, cells)
        s.Minv = np.linalg.inv(M)
    def interp(s, phi):   # linear, internal faces
        return s.wl.reshape((-1,) + (1,)*(phi.ndim - 1))*phi[s.own[:s.nI]] + (1 - s.wl.reshape((-1,) + (1,)*(phi.ndim - 1)))*phi[s.nei]
    def snGradI(s, phi):  # internal faces, orthogonal part (mesh is orthogonal or nearly)
        nhat = s.Sf[:s.nI]/s.magSf[:s.nI, None]
        return (phi[s.nei] - phi[s.own[:s.nI]])/np.einsum('ij,ij->i', nhat, s.dI)
    def snGradB(s, cells, fs, vb, vc):
        nhat = s.Sf[fs]/s.magSf[fs, None]
        return (vb - vc)/np.einsum('ij,ij->i', nhat, s.Cf[fs] - s.C[cells])
    def reconstruct(s, fluxI, fluxB):  # flux = value*|Sf| on faces; returns cell (x, z)
        r = np.zeros((s.nC, 2)); nI = s.nI
        n2 = s.Sf[:nI][:, [0, 2]]/s.magSf[:nI, None]
        np.add.at(r, s.own[:nI], n2*fluxI[:, None]); np.add.at(r, s.nei, n2*fluxI[:, None])
        for (name, fs, cells) in s.bfaces:
            if name in fluxB:
                np.add.at(r, cells, (s.Sf[fs][:, [0, 2]]/s.magSf[fs, None])*fluxB[name][:, None])
        return np.einsum('ijk,ik->ij', s.Minv, r)

def rhoOf(a): a = np.clip(a, 0, 1); return a*RW + (1 - a)*RA
def muOf(a): a = np.clip(a, 0, 1); return a*RW*NUW + (1 - a)*RA*NUA

def zoneMask(case, nC):
    s = open(case + '/constant/polyMesh/cellZones').read()
    m = re.search(r'bodyZone\s*\{.*?cellLabels\s+List<label>\s+(\d+)\s*\(([^)]*)\)', s, re.S)
    z = np.zeros(nC, bool); z[np.fromstring(m.group(2), sep=' ', dtype=int)] = True; return z

def tstr(t): return '%g' % t

def analyse(case, t, dt=0.1):
    Cp = Z + case
    M = Mesh(Cp, tstr(t)); nC = M.nC; pat = M.bnd
    rd = lambda f, nc, tt=t: readField('%s/%s/%s' % (Cp, tstr(tt), f), nc, nC, pat)
    aI, aIb = rd('alphaI', 1); UI, UIb = rd('UInc', 3)
    pH, pHb = rd('pHydIrgh', 1); pD, pDb = rd('pDynD', 1)
    rI = rhoOf(aI); rIb = {k: rhoOf(v) for k, v in aIb.items() if v is not None}
    # ddt(rhoI UI): central difference of the cell content (rhoI UI V)
    def content(tt):
        Mt = Mesh(Cp, tstr(tt)); a, _ = rd('alphaI', 1, tt); u, _ = rd('UInc', 3, tt)
        return (rhoOf(a)[:, None]*u*Mt.V[:, None])[:, [0, 2]]
    ddt = (content(t + dt) - content(t - dt))/(2*dt*M.V[:, None])
    # div(mI uI): flux F = rhoI_f (UI_f . Sf - meshPhi), UI_f linear
    UIf = M.interp(UI); mphi = np.zeros(M.nI)
    mp = '%s/%s/meshPhi' % (Cp, tstr(t))
    if os.path.exists(mp): mphi = readList_surface(mp, M.nI)
    F = M.interp(rI)*(np.einsum('ij,ij->i', UIf, M.Sf[:M.nI]) - mphi)
    conv = np.zeros((nC, 3))
    np.add.at(conv, M.own[:M.nI], F[:, None]*UIf); np.add.at(conv, M.nei, -F[:, None]*UIf)
    for name, fs, cells in M.bfaces:
        ub = UIb.get(name); rb = rIb.get(name)
        if ub is None or rb is None: continue
        Fb = rb*np.einsum('ij,ij->i', ub, M.Sf[fs]); np.add.at(conv, cells, Fb[:, None]*ub)
    conv = (conv/M.V[:, None])[:, [0, 2]]
    # face pressure part G_f = snGrad(pHrgh + pDyn) + ghf snGrad(rhoI), split
    gh = -G*M.C[:, 2]; ghf = -G*M.Cf[:, 2]
    pHyd = pH + rI*gh
    def faceSplit(fsI=True):
        out = {}
        sg = lambda c, cb: (M.snGradI(c), {n: M.snGradB(cells, fs, cb[n], c[cells]) for n, fs, cells in M.bfaces if cb.get(n) is not None})
        pHydb = {n: pHb[n] + rIb[n]*ghf[fs] for n, fs, _ in M.bfaces if pHb.get(n) is not None and n in rIb}
        ghrb = {n: rIb[n]*ghf[fs] for n, fs, _ in M.bfaces if n in rIb}
        sH, sHb = sg(pHyd, pHydb); sD, sDb = sg(pD, pDb)
        sR, sRb = sg(rI, rIb); sRG, sRGb = sg(rI*gh, ghrb); sP, sPb = sg(pH, pHb)
        Ghyd = sP + ghf[:M.nI]*sR
        Ghydb = {n: sPb[n] + ghf[fs]*sRb[n] for n, fs, _ in M.bfaces if n in sPb and n in sRb}
        grav = Ghyd - sH; gravb = {n: Ghydb[n] - sHb[n] for n in Ghydb if n in sHb}
        out['Gf'] = Ghyd + sD
        mI = M.magSf[:M.nI]; mB = {n: M.magSf[fs] for n, fs, _ in M.bfaces}
        out['pHyd'] = M.reconstruct(sH*mI, {n: sHb[n]*mB[n] for n in sHb})
        out['grav'] = M.reconstruct(grav*mI, {n: gravb[n]*mB[n] for n in gravb})
        out['pDyn'] = M.reconstruct(sD*mI, {n: sDb[n]*mB[n] for n in sDb})
        out['Ghyd_f'] = Ghyd; out['GDyn_f'] = sD
        return out
    fsp = faceSplit()
    # viscous: laplacian(mu_f, UI), Gauss linear, internal faces only
    muf = M.interp(muOf(aI)); vis = np.zeros((nC, 3))
    flx = (muf*M.magSf[:M.nI]/np.einsum('ij,ij->i', M.Sf[:M.nI]/M.magSf[:M.nI, None], M.dI))[:, None]*(UI[M.nei] - UI[M.own[:M.nI]])
    np.add.at(vis, M.own[:M.nI], flx); np.add.at(vis, M.nei, -flx)
    vis = (vis/M.V[:, None])[:, [0, 2]]
    terms = dict(ddt=ddt, conv=conv, pHyd=fsp['pHyd'], grav=fsp['grav'], pDyn=fsp['pDyn'], visc=-vis)
    terms['pHyd+grav'] = fsp['pHyd'] + fsp['grav']
    terms['total'] = ddt + conv + fsp['pHyd'] + fsp['grav'] + fsp['pDyn'] - vis
    terms['ddt+conv'] = ddt + conv
    terms['inert+hyd'] = ddt + conv + fsp['pHyd'] + fsp['grav']
    terms['total_x'] = np.c_[terms['total'][:, 0], 0*ddt[:, 0]]
    terms['total_z'] = np.c_[0*ddt[:, 0], terms['total'][:, 1]]
    return M, aI, terms, fsp

def readList_surface(f, nI):
    s = open(f).read(); i = s.index('internalField'); line = s[i:s.index('\n', i)]
    if 'nonuniform' not in line:
        return np.full(nI, float(line.split('uniform')[1].strip(' ;')))
    a = s.index('(', s.index('List', i)); b = s.index('\n)', a)
    return np.fromstring(s[a+1:b], sep=' ')[:nI]

def columns(M):
    # topological columns: group cells by set-up x of their centre (blockMesh grid)
    c0 = np.zeros((M.nC, 3)); cnt = np.zeros(M.nC)
    q = np.array([M.p0[fv].mean(0) for fv in M.F])
    np.add.at(c0, M.own, q); np.add.at(cnt, M.own, 1); np.add.at(c0, M.nei, q[:M.nI]); np.add.at(cnt, M.nei, 1)
    x0 = np.round((c0/cnt[:, None])[:, 0], 6)
    return {x: np.where(x0 == x)[0] for x in np.unique(x0)}

def run(case, times, region):
    print('\n=== %s ===' % case)
    acc = {}; colrows = []; pairs = []
    for t in times:
        M, aI, terms, fsp = analyse(case, t)
        band = (aI > 0.01) & (aI < 0.99)
        sel = band & region(M, case)
        for k, v in terms.items():
            acc.setdefault(k, []).append(np.linalg.norm(v[sel], axis=1))
        # columns: eta error and band G_f
        cols = columns(M); nRow = max(len(c) for c in cols.values())
        zUp = np.full(M.nC, -1e9); zDn = np.full(M.nC, 1e9)
        vert = abs(M.Sf[:, 2]) > 0.7*M.magSf
        for f in np.where(vert)[0]:
            for c in ((M.own[f], M.nei[f]) if f < M.nI else (M.own[f],)):
                zUp[c] = max(zUp[c], M.Cf[f, 2]); zDn[c] = min(zDn[c], M.Cf[f, 2])
        cellH = zUp - zDn; zBot = zDn
        hz = np.zeros(M.nC)
        for f, fv in enumerate(M.F): pass
        for x, cs in cols.items():
            if len(cs) < nRow: continue           # cut by the body
            if not region(M, case)[cs].any(): continue
            zc = M.C[cs, 2]; o = np.argsort(zc); cs = cs[o]
            # cell heights from the moved mesh: z of the top face centre minus
            # z of the bottom face centre (faces with |n_z| > 0.7)
            dz = cellH[cs]
            zb = zBot[cs[0]]
            etaI = zb + (np.clip(aI[cs], 0, 1)*dz).sum()
            xc = M.C[cs, 0].mean(); eH = etaHOS(xc, t)
            k = np.argmin(abs(zc[o] - etaI)); dzb = dz[k]
            # band G_f: vertical internal faces touching band cells of the column, signed upward, max |.|
            bc = set(cs[band[cs]]) | {cs[k]}
            best = 0.0
            for f in range(M.nI):
                if abs(M.Sf[f, 2]) < 0.9*M.magSf[f]: continue
                if M.own[f] in bc or M.nei[f] in bc:
                    if M.own[f] in cs and M.nei[f] in cs:
                        g = fsp['Gf'][f]*np.sign(M.Sf[f, 2])
                        if abs(g) > abs(best): best = g
            colrows.append((t, xc, eH - etaI))
            pairs.append((RW*G*(eH - etaI)/dzb, best))
    print('rI terms in band cells (N/m3), %d cell samples over %s s:' % (len(acc['total'][0]) + len(acc['total'][1]) + len(acc['total'][2]), times))
    print('  %-10s %10s %10s' % ('term', 'max', 'rms'))
    for k in ('ddt', 'conv', 'pHyd', 'grav', 'pHyd+grav', 'pDyn', 'visc', 'ddt+conv', 'inert+hyd', 'total', 'total_x', 'total_z'):
        v = np.concatenate(acc[k]); print('  %-10s %10.1f %10.1f' % (k, v.max(), np.sqrt((v*v).mean())))
    e = np.array([r[2] for r in colrows])*1e3
    print('etaHOS - etaI per column: n %d, max |.| %.2f mm, rms %.2f mm, mean %+.2f mm' % (len(e), abs(e).max(), np.sqrt((e*e).mean()), e.mean()))
    P = np.array(pairs); x, y = P[:, 0], P[:, 1]
    A = np.vstack([x, np.ones_like(x)]).T; (sl, ic), res, *_ = np.linalg.lstsq(A, y, rcond=None)
    r2 = 1 - ((y - (sl*x + ic))**2).sum()/((y - y.mean())**2).sum()
    print('band G_f vs rho_w g (etaHOS - etaI)/dz: slope %.3f, intercept %.1f, R2 %.3f (n %d); |x| rms %.1f, |G| rms %.1f N/m3'
          % (sl, ic, r2, len(x), np.sqrt((x*x).mean()), np.sqrt((y*y).mean())))

def regZD(M, case):
    z = zoneMask(Z + case, M.nC); ax = abs(M.C[:, 0])
    return (~z) & (ax > 0.5) & (ax < 1.5)

def regB05(M, case):
    for name, fs, cells in M.bfaces:
        if name == 'box':
            cb = M.Cf[fs]
    d = np.min(np.linalg.norm(M.C[:, None, [0, 2]] - cb[None, :, [0, 2]], axis=2), axis=1)
    return d <= 0.3

if __name__ == '__main__':
    times = [15.0, 19.0, 22.9]
    run('ZD', times, regZD)
    run('B05', times, regB05)
