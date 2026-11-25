"""Reader for HOS-NWT modes_HOS_SWENSE.dat (ASCII, 'velocity' layout).

Layout (from HOS-NWT sources/IO/output.f90 and Grid2Grid hosNWT.inc):
  formatted direct-access file, fixed records of 18*n1 characters, no newlines,
  each value written as ES17.10 plus one space.
  Block 0 (8*n2 records): header in record 1:
      n1, n3_add, n2, dt_out*, T_stop*, xlen*, ylen*, h[m]   (* = non-dimensional)
  Block j >= 1: state at t* = (j-1)*dt_out*, 8 records per y-mode iy:
      modeX, modeY, modeZ, modet, modeFS, modeFSt   (n1 values)
      modeAdd, modeAddt                              (n3_add values, rest is padding)

All values are non-dimensional: lengths / h, time / sqrt(h/g), g = 1, bottom at z = -1.

Usage: python3 hos_nwt_reader.py modes_HOS_SWENSE.dat
"""
import sys
import numpy as np

NAMES = ["modeX", "modeY", "modeZ", "modet", "modeFS", "modeFSt", "modeAdd", "modeAddt"]


class HOSNWTModes:
    def __init__(self, path, l_add=2):
        self.raw = open(path).read()
        v = [float(self.raw[18 * i:18 * i + 17]) for i in range(8)]
        self.n1, self.n3_add, self.n2 = int(v[0]), int(v[1]), int(v[2])
        self.dt, self.T_stop, self.xlen, self.ylen, self.h = v[3:8]
        self.reclen = 18 * self.n1
        self.nblocks = len(self.raw) // (self.reclen * 8 * self.n2)
        self.ntimes = self.nblocks - 1
        self.kx = np.arange(self.n1) * np.pi / self.xlen
        self.ky = np.arange(self.n2) * np.pi / self.ylen if self.n2 > 1 else np.zeros(1)
        # wavemaker (additional) vertical wavenumbers; l_add is a compile-time
        # parameter in common_vars.f90 (default 2), not stored in the file
        self.kx_add = (2 * np.arange(1, self.n3_add + 1) - 1) * np.pi / (2 * l_add)
        self.Lscale, self.Tscale = self.h, np.sqrt(self.h / 9.81)

    def _record(self, irec, n):
        s = self.raw[irec * self.reclen: irec * self.reclen + 18 * n]
        return np.array([float(s[18 * j:18 * j + 17]) for j in range(n)])

    def modes(self, itime):
        """Modes at t* = itime*dt (itime = 0 .. ntimes-1). Arrays shaped (n, n2)."""
        base = (itime + 1) * 8 * self.n2
        out = {nm: [] for nm in NAMES}
        for iy in range(self.n2):
            for j, nm in enumerate(NAMES):
                n = self.n3_add if j >= 6 else self.n1
                out[nm].append(self._record(base + 8 * iy + j, n))
        return {nm: np.array(v).T for nm, v in out.items()}

    def fields_2d(self, itime, x, z, include_add=True):
        """eta, eta_t on x, and u, w, phi_t at (x, z) for a 2D run (n2 = 1). Non-dimensional."""
        m = self.modes(itime)
        k = self.kx
        C = np.cos(np.outer(x, k)); S = np.sin(np.outer(x, k))
        eta = C @ m["modeFS"][:, 0]
        eta_t = C @ m["modeFSt"][:, 0]
        ch = np.cosh(k * (z + 1)) / np.cosh(k)
        sh = np.ones_like(k); sh[1:] = np.sinh(k[1:] * (z + 1)) / np.sinh(k[1:])
        u = S @ (m["modeX"][:, 0] * ch)
        w = C @ (m["modeZ"][:, 0] * sh)
        phit = C @ (m["modet"][:, 0] * ch)
        if include_add:
            ka = self.kx_add[None, :]; X = np.asarray(x)[:, None]; L = self.xlen
            den = 1 + np.exp(-2 * ka * L)
            csh = (np.exp(-ka * X) + np.exp(-ka * (2 * L - X))) / den   # cosh(ka(L-x))/cosh(ka L)
            snh = (np.exp(-ka * X) - np.exp(-ka * (2 * L - X))) / den   # sinh(ka(L-x))/cosh(ka L)
            a, at = m["modeAdd"][:, 0], m["modeAddt"][:, 0]
            cz, sz = np.cos(self.kx_add * (z + 1)), np.sin(self.kx_add * (z + 1))
            u += -(snh * (a * self.kx_add * cz)).sum(1)
            w += -(csh * (a * self.kx_add * sz)).sum(1)
            phit += (csh * (at * cz)).sum(1)
        return dict(eta=eta, eta_t=eta_t, u=u, w=w, phi_t=phit)


if __name__ == "__main__":
    f = HOSNWTModes(sys.argv[1])
    print(f"n1={f.n1} n2={f.n2} n3_add={f.n3_add} dt*={f.dt:.5g} T_stop*={f.T_stop:.5g} "
          f"xlen*={f.xlen} ylen*={f.ylen} h={f.h} m, {f.ntimes} saved times")
    m = f.modes(f.ntimes // 2)
    i = np.argmax(np.abs(m["modeFS"][1:, 0])) + 1
    print(f"mid-run dominant mode i={i}: modeX/modeZ = {m['modeX'][i,0]/m['modeZ'][i,0]:.6f}, "
          f"-1/tanh(k) = {-1/np.tanh(f.kx[i]):.6f}")
