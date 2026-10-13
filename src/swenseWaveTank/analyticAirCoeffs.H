/*---------------------------------------------------------------------------*\
  analyticAirCoeffs.H  (header-only, standard C++11, no OpenFOAM types)

  Air-side modal coefficients from the HOS modes, orders 1, 2, 3.
  Derivation: air_coefficients_third_order.pdf, Section 15.7.
  Same algorithm as air_coeffs.m (methods 'linear', 'second', 'third').

  Everything is NON-DIMENSIONAL (HOS units: length h, time sqrt(h/g)).

  Basis (tank with walls):
      eta   = sum_mn e_mn  cos(kx_m x) cos(ky_n y)      e  = modeFS
      eta_t = sum_mn et_mn cos(kx_m x) cos(ky_n y)      et = modeFSt
      psi   = w0 z + sum_(m,n)!=(0,0) A_mn cos(kx_m x) cos(ky_n y) exp(-k_mn z)
      kx_m = m pi / Lx,  ky_n = n pi / Ly,  k_mn = sqrt(kx^2 + ky^2)
  Storage: (m,n) -> m*My + n,  m = 0..Mx-1, n = 0..My-1  (My = 1 in 2D).
  Pass in the modes AFTER the spectral cut-off (zero above k_cut, or
  truncate Mx, My to the retained block).

  Algorithm:
      order 1:  w0 = et_00,  A = -et / k
      one pass: on z = 0 from the current A
                  psi_x  = -sum kx A   sin cos      psi_y  = -sum ky A   cos sin
                  psi_xz = +sum kx k A sin cos      psi_yz = +sum ky k A cos sin
                Gx = eta psi_x  [+ 1/2 eta^2 psi_xz  for order 3]
                Gy = eta psi_y  [+ 1/2 eta^2 psi_yz  for order 3]
                gx = projection of Gx on sin(kx x) cos(ky y)
                gy = projection of Gy on cos(kx x) sin(ky y)
                R  = kx gx + ky gy                      ( = div G as cos cos )
                A  = -(et + R) / k
      order 2:  one pass (G without the eta^2 term)
      order 3:  two passes (full G in both)
  No derivatives of eta are needed.

  Projection is exact: midpoint grid with Nq = 4 M points per direction
  (DCT-II / DST-II orthogonality, no aliasing for the cubic eta^2 psi term).
  Cost per pass: O(Nqx Nqy (Mx + My)) with separable tables.
\*---------------------------------------------------------------------------*/

#ifndef analyticAirCoeffs_H
#define analyticAirCoeffs_H

#include <vector>
#include <cmath>
#include <stdexcept>

namespace airAnalytic
{

typedef std::vector<double> vec;

// f(i,j) = sum_mn c[m*My+n] X[i*Mx+m] Y[j*My+n]   (separable synthesis)
inline void synth
(
    int Mx, int My, int Nx, int Ny,
    const vec& c, const vec& X, const vec& Y, vec& f
)
{
    vec T(Nx*My, 0.0);
    for (int i = 0; i < Nx; ++i)
        for (int m = 0; m < Mx; ++m)
        {
            const double x = X[i*Mx + m];
            if (x == 0.0) continue;
            for (int n = 0; n < My; ++n) T[i*My + n] += x*c[m*My + n];
        }
    f.assign(Nx*Ny, 0.0);
    for (int i = 0; i < Nx; ++i)
        for (int j = 0; j < Ny; ++j)
        {
            double s = 0.0;
            for (int n = 0; n < My; ++n) s += T[i*My + n]*Y[j*My + n];
            f[i*Ny + j] = s;
        }
}

// c[m*My+n] = wx[m] wy[n] sum_ij f(i,j) X[i*Mx+m] Y[j*My+n]   (projection)
inline void project
(
    int Mx, int My, int Nx, int Ny,
    const vec& f, const vec& X, const vec& Y,
    const vec& wx, const vec& wy, vec& c
)
{
    vec T(Nx*My, 0.0);
    for (int i = 0; i < Nx; ++i)
        for (int n = 0; n < My; ++n)
        {
            double s = 0.0;
            for (int j = 0; j < Ny; ++j) s += f[i*Ny + j]*Y[j*My + n];
            T[i*My + n] = s;
        }
    c.assign(Mx*My, 0.0);
    for (int m = 0; m < Mx; ++m)
        for (int n = 0; n < My; ++n)
        {
            double s = 0.0;
            for (int i = 0; i < Nx; ++i) s += X[i*Mx + m]*T[i*My + n];
            c[m*My + n] = wx[m]*wy[n]*s;
        }
}

// Midpoint-grid tables for one direction: C = cos(k_m x_i), S = sin(k_m x_i),
// cosine weights (1/N for m = 0, 2/N else), sine weights (2/N).
// N = 1 (2D, y direction) gives C = 1, S = 0, weights 1.
inline void tables
(
    int M, int N, double L, vec& C, vec& S, vec& wc, vec& ws
)
{
    C.assign(N*M, 0.0); S.assign(N*M, 0.0); wc.assign(M, 0.0); ws.assign(M, 0.0);
    for (int i = 0; i < N; ++i)
    {
        const double x = (i + 0.5)*L/N;
        for (int m = 0; m < M; ++m)
        {
            const double k = m*M_PI/L;
            C[i*M + m] = std::cos(k*x);
            S[i*M + m] = std::sin(k*x);
        }
    }
    for (int m = 0; m < M; ++m)
    {
        wc[m] = (N == 1) ? 1.0 : ((m == 0) ? 1.0/N : 2.0/N);
        ws[m] = (N == 1) ? 1.0 : 2.0/N;
    }
}

// Main entry. e, et: Mx*My modes. order 1, 2 or 3. Returns A (A[0] = 0), w0.
inline void analyticAirCoeffs
(
    int Mx, int My, double Lx, double Ly,
    const vec& e, const vec& et, int order,
    vec& A, double& w0
)
{
    if (order < 1 || order > 3)
        throw std::invalid_argument("analyticAirCoeffs: order must be 1, 2 or 3");

    const int NM = Mx*My;
    vec kx(NM), ky(NM), k(NM);
    for (int m = 0; m < Mx; ++m)
        for (int n = 0; n < My; ++n)
        {
            const int id = m*My + n;
            kx[id] = m*M_PI/Lx;
            ky[id] = (My > 1) ? n*M_PI/Ly : 0.0;
            k[id] = std::sqrt(kx[id]*kx[id] + ky[id]*ky[id]);
        }

    // Order 1 (eq. Aet1): w0 = et_00, A = -et/k
    w0 = et[0];
    A.assign(NM, 0.0);
    for (int id = 1; id < NM; ++id) A[id] = -et[id]/k[id];
    if (order == 1) return;

    // Midpoint grids and tables
    const int Nx = 4*Mx;
    const int Ny = (My > 1) ? 4*My : 1;
    vec Cx, Sx, wcx, wsx, Cy, Sy, wcy, wsy;
    tables(Mx, Nx, Lx, Cx, Sx, wcx, wsx);
    tables(My, Ny, (My > 1 ? Ly : 1.0), Cy, Sy, wcy, wsy);

    vec eta;
    synth(Mx, My, Nx, Ny, e, Cx, Cy, eta);           // eta on the grid

    const int npass = (order == 2) ? 1 : 2;
    const bool quad = (order == 3);                    // 1/2 eta^2 term
    vec c(NM), px, py, pxz, pyz, Gx(Nx*Ny), Gy(Nx*Ny), gx, gy;

    for (int pass = 0; pass < npass; ++pass)
    {
        // psi_x = -sum kx A sin cos ; psi_xz = +sum kx k A sin cos
        for (int id = 0; id < NM; ++id) c[id] = -kx[id]*A[id];
        synth(Mx, My, Nx, Ny, c, Sx, Cy, px);
        for (int id = 0; id < NM; ++id) c[id] = kx[id]*k[id]*A[id];
        synth(Mx, My, Nx, Ny, c, Sx, Cy, pxz);
        if (My > 1)
        {
            // psi_y = -sum ky A cos sin ; psi_yz = +sum ky k A cos sin
            for (int id = 0; id < NM; ++id) c[id] = -ky[id]*A[id];
            synth(Mx, My, Nx, Ny, c, Cx, Sy, py);
            for (int id = 0; id < NM; ++id) c[id] = ky[id]*k[id]*A[id];
            synth(Mx, My, Nx, Ny, c, Cx, Sy, pyz);
        }

        // G = eta grad psi0 (+ 1/2 eta^2 grad d_z psi0)
        for (int q = 0; q < Nx*Ny; ++q)
        {
            const double h2 = quad ? 0.5*eta[q]*eta[q] : 0.0;
            Gx[q] = eta[q]*px[q] + h2*pxz[q];
            if (My > 1) Gy[q] = eta[q]*py[q] + h2*pyz[q];
        }

        // R = kx [Gx]_{sin cos} + ky [Gy]_{cos sin}
        project(Mx, My, Nx, Ny, Gx, Sx, Cy, wsx, wcy, gx);
        if (My > 1) project(Mx, My, Nx, Ny, Gy, Cx, Sy, wcx, wsy, gy);

        // A = -(et + R)/k  (eqs. Aet2, Aet3); w0 unchanged
        for (int id = 1; id < NM; ++id)
        {
            double R = kx[id]*gx[id];
            if (My > 1) R += ky[id]*gy[id];
            A[id] = -(et[id] + R)/k[id];
        }
    }
}

} // namespace airAnalytic

#endif
