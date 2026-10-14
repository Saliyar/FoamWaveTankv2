/*---------------------------------------------------------------------------*\
  Air2AirEnhanced implementation (see Air2AirEnhanced.H)
\*---------------------------------------------------------------------------*/

#include "Air2AirEnhanced.H"
#include "analyticAirCoeffs.H"
#include "error.H"
#include "messageStream.H"
#include "Switch.H"

#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <limits>

namespace
{
const int nSaveRecords = 6;
const int nBandRecords = 3;
const double dtRecord = 0.05;          // record spacing of the HOS data [s]
}

namespace Foam
{

// * * * * * * * * * * * * * * File reading  * * * * * * * * * * * * * * * * //

// 0-based record rec: n fixed-width values of 18 characters
void Air2AirEnhanced::readRecord(long rec, int n, std::vector<scalar>& v) const
{
    std::string buf(18*n, ' ');
    in_.clear();
    in_.seekg(rec*reclen_);
    in_.read(&buf[0], 18*n);
    if (!in_)
    {
        FatalErrorInFunction
            << "Cannot read record " << rec << " of " << file_
            << exit(FatalError);
    }
    v.resize(n);
    for (int i = 0; i < n; ++i)
    {
        std::string f = buf.substr(18*i, 18);
        for (char& c : f) if (c == 'D' || c == 'd') c = 'E';   // Fortran D exponent
        v[i] = std::strtod(f.c_str(), nullptr);
    }
}

// modeFS and modeFSt of time index it, block Mx x My, (m,n) -> m*My + n
void Air2AirEnhanced::readModes
(
    int it, std::vector<scalar>& e, std::vector<scalar>& et, int Mx, int My
) const
{
    e.assign(Mx*My, 0.0);
    et.assign(Mx*My, 0.0);
    std::vector<scalar> r;
    // block it+1 (block 0 is the header), 8 records per y-mode:
    // modeX, modeY, modeZ, modet, modeFS, modeFSt, modeAdd, modeAddt
    const long base = long(it + 1)*8*n2_;
    for (int n = 0; n < My; ++n)
    {
        readRecord(base + 8*n + 4, Mx, r);
        for (int m = 0; m < Mx; ++m) e[m*My + n] = r[m];
        readRecord(base + 8*n + 5, Mx, r);
        for (int m = 0; m < Mx; ++m) et[m*My + n] = r[m];
    }
}

// * * * * * * * * * * * * * * * Cut-off * * * * * * * * * * * * * * * * * * //

void Air2AirEnhanced::setupCutoff(const dictionary& dict)
{
    cutoff_ = dict.lookupOrDefault<word>("cutoff", word("auto"));
    const scalar C = dict.lookupOrDefault<scalar>("cutoffC", 0.5);
    const scalar Kp = dict.lookupOrDefault<scalar>("cutoffKp", 4.0);
    const scalar Lc = dict.lookupOrDefault<scalar>("cutoffLength", 1.0);
    const scalar etaOverride = dict.lookupOrDefault<scalar>("cutoffEtaMax", -1.0);
    const int stride = std::max(1, int(dict.lookupOrDefault<label>("cutoffScanStride", 1)));

    // Scan the whole run once: max|eta| in physical space and the
    // time-averaged mode energy (uncut modeFS)
    const int Gx = 2*n1_, Gy = (n2_ > 1) ? 2*n2_ : 1;
    std::vector<scalar> Cx(Gx*n1_), Cy(Gy*n2_);
    for (int i = 0; i < Gx; ++i)
        for (int m = 0; m < n1_; ++m)
            Cx[i*n1_ + m] = std::cos(m*M_PI*(i + 0.5)/Gx);
    for (int j = 0; j < Gy; ++j)
        for (int n = 0; n < n2_; ++n)
            Cy[j*n2_ + n] = (n2_ > 1) ? std::cos(n*M_PI*(j + 0.5)/Gy) : 1.0;

    etaMax_ = 0;
    std::vector<scalar> E(n1_*n2_, 0.0), e, et, T(Gx*n2_);
    int nscan = 0;
    for (int it = 0; it < nt_; it += stride)
    {
        readModes(it, e, et, n1_, n2_);
        for (int q = 0; q < n1_*n2_; ++q) E[q] += e[q]*e[q];
        std::fill(T.begin(), T.end(), 0.0);
        for (int i = 0; i < Gx; ++i)
            for (int m = 0; m < n1_; ++m)
            {
                const scalar c = Cx[i*n1_ + m];
                for (int n = 0; n < n2_; ++n) T[i*n2_ + n] += c*e[m*n2_ + n];
            }
        for (int i = 0; i < Gx; ++i)
            for (int j = 0; j < Gy; ++j)
            {
                scalar s = 0;
                for (int n = 0; n < n2_; ++n) s += T[i*n2_ + n]*Cy[j*n2_ + n];
                etaMax_ = std::max(etaMax_, std::fabs(s));
            }
        ++nscan;
    }
    if (etaOverride > 0) etaMax_ = etaOverride/h_;

    auto kmn = [&](int m, int n)
    {
        const scalar kx = m*M_PI/xlen_, ky = (n2_ > 1) ? n*M_PI/ylen_ : 0.0;
        return std::sqrt(kx*kx + ky*ky);
    };

    if (cutoff_ == "none")
    {
        kcut_ = std::numeric_limits<scalar>::max();
    }
    else if (cutoff_ == "auto")
    {
        if (etaMax_ <= 0)
        {
            FatalErrorInFunction << "cutoff auto: max|eta| is zero" << exit(FatalError);
        }
        kcut_ = C/etaMax_;
    }
    else if (cutoff_ == "kp")
    {
        // peak of the time-averaged spectrum, binned by |k| in steps of pi/xlen
        const scalar dk = M_PI/xlen_;
        std::vector<scalar> bin(n1_ + n2_ + 2, 0.0);
        for (int m = 0; m < n1_; ++m)
            for (int n = 0; n < n2_; ++n)
            {
                if (m == 0 && n == 0) continue;
                const int b = std::min(int(bin.size()) - 1, int(std::floor(kmn(m, n)/dk + 0.5)));
                bin[b] += E[m*n2_ + n];
            }
        const int bp = int(std::max_element(bin.begin() + 1, bin.end()) - bin.begin());
        kcut_ = Kp*bp*dk;
    }
    else if (cutoff_ == "wavelength")
    {
        kcut_ = 2*M_PI/(Lc/h_);
    }
    else
    {
        FatalErrorInFunction
            << "Air2Air cutoff '" << cutoff_ << "' unknown; use auto, kp, wavelength or none"
            << exit(FatalError);
    }

    // smallest block holding every retained mode
    Mx_ = (kcut_ >= kmn(n1_ - 1, 0)) ? n1_
        : std::min(n1_, int(std::floor(kcut_*xlen_/M_PI)) + 1);
    My_ = (n2_ == 1) ? 1
        : ((kcut_ >= kmn(0, n2_ - 1)) ? n2_
          : std::min(n2_, int(std::floor(kcut_*ylen_/M_PI)) + 1));
    Mx_ = std::max(Mx_, 2);
    keep_.assign(Mx_*My_, 0);
    int nkeep = 0;
    scalar Ekeep = 0, Etot = 0, kmax = 0;
    for (int m = 0; m < n1_; ++m)
        for (int n = 0; n < n2_; ++n)
        {
            if (m == 0 && n == 0) continue;
            Etot += E[m*n2_ + n];
            kmax = std::max(kmax, kmn(m, n));
            if (m < Mx_ && n < My_ && kmn(m, n) <= kcut_)
            {
                keep_[m*My_ + n] = 1;
                Ekeep += E[m*n2_ + n];
                ++nkeep;
            }
        }
    keep_[0] = 1;

    Info<< "Air2Air spectral cut-off: " << cutoff_ << nl
        << "    max|eta| over " << nscan << " records = " << etaMax_*h_ << " m" << nl;
    if (cutoff_ != "none")
    {
        Info<< "    k_cut = " << kcut_/h_ << " 1/m,  L_cut = " << 2*M_PI*h_/kcut_ << " m" << nl
            << "    k*max|eta|: " << kmax*etaMax_ << " before, " << kcut_*etaMax_ << " after" << nl;
    }
    else
    {
        Info<< "    k_max*max|eta| = " << kmax*etaMax_
            << (kmax*etaMax_ > 1 && order_ >= 2
                ? "  WARNING: > 1, the expansion is not valid for the short modes" : "")
            << nl;
    }
    Info<< "    retained " << nkeep << " of " << n1_*n2_ - 1 << " modes, energy removed "
        << 100*(Etot > 0 ? 1 - Ekeep/Etot : 0) << " %" << endl;
}

// * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * //

Air2AirEnhanced::Air2AirEnhanced(const string& modesFile, const dictionary& dict)
:
    file_(modesFile),
    in_(modesFile.c_str(), std::ios::in | std::ios::binary),
    n1_(0), n2_(0), n3add_(0), nt_(0), reclen_(0),
    dt_(0), xlen_(0), ylen_(0), h_(0),
    g_(dict.lookupOrDefault<scalar>("g", 9.81)),
    rhoAir_(dict.lookupOrDefault<scalar>("rhoAir", 1.225)),
    order_(dict.lookupOrDefault<label>("order", 2)),
    bernoulli_(false),
    verbose_(dict.lookupOrDefault<Switch>("verbose", Switch(false))),
    kcut_(0), etaMax_(0), Mx_(0), My_(0),
    printed_(false),
    w0_(0), w0t_(0),
    updated_(false)
{
    if (!in_.good())
    {
        FatalErrorInFunction << "Cannot open modes file " << file_ << exit(FatalError);
    }
    if (order_ < 0 || order_ > 3)
    {
        FatalErrorInFunction << "Air2Air order must be 0, 1, 2 or 3" << exit(FatalError);
    }
    const word pm = dict.lookupOrDefault<word>("pressureModel", word("linear"));
    if (pm == "bernoulli") bernoulli_ = true;
    else if (pm != "linear")
    {
        FatalErrorInFunction << "pressureModel must be linear or bernoulli" << exit(FatalError);
    }

    // Header: record 0 holds n1, n3_add, n2, dt*, T_stop*, xlen*, ylen*, h[m]
    reclen_ = 18*8;
    std::vector<scalar> hd;
    readRecord(0, 8, hd);
    n1_ = int(hd[0] + 0.5); n3add_ = int(hd[1] + 0.5); n2_ = int(hd[2] + 0.5);
    dt_ = hd[3]; xlen_ = hd[5]; ylen_ = hd[6]; h_ = hd[7];
    reclen_ = 18L*n1_;
    in_.seekg(0, std::ios::end);
    const long size = long(in_.tellg());
    nt_ = int(size/(reclen_*8*n2_)) - 1;
    if (n1_ < 2 || nt_ < 2 || h_ <= 0)
    {
        FatalErrorInFunction << "Bad modes file header in " << file_ << exit(FatalError);
    }

    setupCutoff(dict);

    cache_.assign(nSaveRecords, slot());
    for (slot& c : cache_) c.index = -1;

    Info<< "Air2AirEnhanced: order " << order_
        << ", pressureModel " << (bernoulli_ ? "bernoulli" : "linear") << nl
        << "    modes file " << file_ << nl
        << "    n1 = " << n1_ << ", n2 = " << n2_ << ", records = " << nt_
        << ", dt_HOS = " << getDt() << " s, h = " << h_ << " m" << nl
        << "    tank " << xlen_*h_ << " x " << ylen_*h_ << " m" << endl;
}

Air2AirEnhanced::~Air2AirEnhanced()
{}

// * * * * * * * * * * * * * * Time interpolation (Grid2Grid scheme) * * * * * //
//
// Grid2Grid (hosVol2Vol.inc: nInterp = interpolationOrder + 1 = 4, nSaveT =
// nInterp + 2 = 6, nBand = floor(nSaveT/2) = 3, correctIdx = ceiling(t/dt),
// window start = correctIdx - nBand clamped to the record range) interpolates
// the HOS data with a cubic B-spline over a window of six records, with
// automatic knots and not-a-knot end conditions (bspline_oo_module.f90, db1ink
// with iknot = 0). Air2Air uses the same window, nodes and knot vector for
// every mode coefficient A_k and for w0.

namespace
{


// Cubic B-spline basis B_j(x), j = 0..5, and derivatives, on a knot vector
// of 10 knots (order 4); x in [kv[0], kv[9]]
void bsplBasis(const scalar* kv, const scalar x, scalar* B, scalar* dB)
{
    const int nK = 10;
    scalar N[10][4];
    for (int i = 0; i < nK; ++i)
        for (int p = 0; p < 4; ++p) N[i][p] = 0;

    for (int i = 0; i < nK - 1; ++i)
    {
        if (kv[i] < kv[i+1])
        {
            if ((x >= kv[i] && x < kv[i+1]) || (x == kv[nK-1] && kv[i+1] == kv[nK-1]))
            {
                N[i][0] = 1;
            }
        }
    }

    for (int p = 1; p <= 3; ++p)
    {
        for (int i = 0; i < nK - 1 - p; ++i)
        {
            scalar a = 0, b = 0;
            const scalar d1 = kv[i+p] - kv[i];
            const scalar d2 = kv[i+p+1] - kv[i+1];
            if (d1 > 0) a = (x - kv[i])/d1*N[i][p-1];
            if (d2 > 0) b = (kv[i+p+1] - x)/d2*N[i+1][p-1];
            N[i][p] = a + b;
        }
    }

    for (int j = 0; j < 6; ++j)
    {
        B[j] = N[j][3];
        scalar d = 0;
        const scalar d1 = kv[j+3] - kv[j];
        const scalar d2 = kv[j+4] - kv[j+1];
        if (d1 > 0) d += N[j][2]/d1;
        if (d2 > 0) d -= N[j+1][2]/d2;
        dB[j] = 3*d;
    }
}

// Inverse of a 6x6 matrix by Gauss-Jordan elimination with partial pivoting
void invert6(scalar M[6][6], scalar Minv[6][6])
{
    scalar a[6][12];
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 12; ++j)
            a[i][j] = (j < 6) ? M[i][j] : ((j - 6) == i ? 1.0 : 0.0);

    for (int c = 0; c < 6; ++c)
    {
        int piv = c;
        for (int r = c + 1; r < 6; ++r)
            if (std::fabs(a[r][c]) > std::fabs(a[piv][c])) piv = r;
        for (int j = 0; j < 12; ++j) std::swap(a[c][j], a[piv][j]);
        const scalar d = a[c][c];
        for (int j = 0; j < 12; ++j) a[c][j] /= d;
        for (int r = 0; r < 6; ++r)
        {
            if (r == c) continue;
            const scalar f = a[r][c];
            for (int j = 0; j < 12; ++j) a[r][j] -= f*a[c][j];
        }
    }

    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            Minv[i][j] = a[i][j+6];
}

}


// * * * * * * * * * * * * * * Coefficients  * * * * * * * * * * * * * * * * //

// Coefficients of record it; 'other' is the record that must stay cached
const Air2AirEnhanced::slot& Air2AirEnhanced::coeffs(int it, int iLeft) const
{
    for (const slot& c : cache_)
        if (c.index == it) return c;

    // Free a slot outside the current window [iLeft, iLeft + nSaveRecords - 1]
    slot* freeSlot = nullptr;
    for (slot& c : cache_)
    {
        if (c.index < iLeft || c.index > iLeft + nSaveRecords - 1)
        {
            freeSlot = &c;
            break;
        }
    }
    if (!freeSlot)
    {
        FatalErrorInFunction << "no free record slot in the Air2Air window"
            << exit(FatalError);
    }
    slot& c = *freeSlot;

    std::vector<scalar> e, et;
    readModes(it, e, et, Mx_, My_);
    for (int q = 0; q < Mx_*My_; ++q)
        if (!keep_[q]) { e[q] = 0; et[q] = 0; }

    if (order_ == 0)
    {
        c.A.assign(Mx_*My_, 0.0);
        c.w0 = 0;
    }
    else
    {
        airAnalytic::analyticAirCoeffs
        (
            Mx_, My_, xlen_, (n2_ > 1 ? ylen_ : 1.0), e, et, order_, c.A, c.w0
        );
    }
    c.index = it;

    if (verbose_ && !printed_)
    {
        scalar amax = 0;
        for (scalar a : c.A) amax = std::max(amax, std::fabs(a));
        const int st = My_;     // stride of m in the (m,n) layout
        Info<< "Air2AirEnhanced: record " << it << " order " << order_
            << " (non-dimensional) A[1..3] = "
            << c.A[1*st] << " " << (Mx_ > 2 ? c.A[2*st] : 0) << " "
            << (Mx_ > 3 ? c.A[3*st] : 0) << "  max|A| = " << amax << endl;
        printed_ = true;
    }
    return c;
}

void Air2AirEnhanced::update(scalar t)
{
    const scalar Ts = std::sqrt(h_/g_);

    // Grid2Grid window: correctIdx = ceiling(t/dt), start = correctIdx - nBand,
    // clamped to the records (hosVol2Vol.inc)
    int iLeft = int(std::ceil(std::max(t, scalar(0))/dtRecord)) - nBandRecords;
    iLeft = std::max(iLeft, 0);
    iLeft = std::min(iLeft, nt_ - nSaveRecords);

    const int N = Mx_*My_;
    std::vector<const slot*> win(nSaveRecords);
    for (int j = 0; j < nSaveRecords; ++j)
    {
        win[j] = &coeffs(iLeft + j, iLeft);
    }

    // Not-a-knot cubic B-spline nodes tau_j = (iLeft + j) dt, knots 4 x tau0,
    // tau2, tau3, 4 x tau5
    scalar tau[6];
    for (int j = 0; j < 6; ++j) tau[j] = (iLeft + j)*dtRecord;
    scalar kv[10] =
    {
        tau[0], tau[0], tau[0], tau[0], tau[2], tau[3],
        tau[5], tau[5], tau[5], tau[5]
    };

    scalar M[6][6], Minv[6][6];
    for (int i = 0; i < 6; ++i)
    {
        scalar B[6], dB[6];
        bsplBasis(kv, tau[i], B, dB);
        for (int j = 0; j < 6; ++j) M[i][j] = B[j];
    }
    invert6(M, Minv);

    scalar Bt[6], dBt[6];
    bsplBasis(kv, t, Bt, dBt);

    A_.resize(N);
    At_.resize(N);
    for (int q = 0; q < N; ++q)
    {
        scalar A = 0, dA = 0;
        for (int j = 0; j < 6; ++j)
        {
            scalar cj = 0;
            for (int i = 0; i < 6; ++i) cj += Minv[j][i]*win[i]->A[q];
            A += cj*Bt[j];
            dA += cj*dBt[j];
        }
        A_[q] = A;
        At_[q] = Ts*dA;            // d/dtau, tau = t/Ts
    }

    {
        scalar A = 0, dA = 0;
        for (int j = 0; j < 6; ++j)
        {
            scalar cj = 0;
            for (int i = 0; i < 6; ++i) cj += Minv[j][i]*win[i]->w0;
            A += cj*Bt[j];
            dA += cj*dBt[j];
        }
        w0_ = A;
        w0t_ = Ts*dA;
    }

    updated_ = true;
}


void Air2AirEnhanced::recordCoefficients
(
    int it,
    std::vector<scalar>& A,
    scalar& w0
) const
{
    const slot& c = coeffs(it, it);
    A = c.A;
    w0 = c.w0;
}

// * * * * * * * * * * * * * * Field evaluation  * * * * * * * * * * * * * * //

void Air2AirEnhanced::evaluate
(
    scalar xs, scalar ys, scalar zs,
    scalar& u, scalar& v, scalar& w, scalar& phi, scalar& phit
) const
{
    if (!updated_)
    {
        FatalErrorInFunction << "Air2AirEnhanced: call update(t) first" << exit(FatalError);
    }
    // cos/sin of multiples by recurrence
    std::vector<scalar> cx(Mx_), sx(Mx_), cy(My_), sy(My_);
    const scalar ax = M_PI*xs/xlen_, ca = std::cos(ax), sa = std::sin(ax);
    cx[0] = 1; sx[0] = 0;
    for (int m = 1; m < Mx_; ++m)
    {
        cx[m] = cx[m-1]*ca - sx[m-1]*sa;
        sx[m] = sx[m-1]*ca + cx[m-1]*sa;
    }
    cy[0] = 1; sy[0] = 0;
    if (My_ > 1)
    {
        const scalar ay = M_PI*ys/ylen_, cb = std::cos(ay), sb = std::sin(ay);
        for (int n = 1; n < My_; ++n)
        {
            cy[n] = cy[n-1]*cb - sy[n-1]*sb;
            sy[n] = sy[n-1]*cb + cy[n-1]*sb;
        }
    }

    u = 0; v = 0; w = w0_; phi = w0_*zs; phit = w0t_*zs;
    for (int m = 0; m < Mx_; ++m)
    {
        const scalar kx = m*M_PI/xlen_;
        for (int n = 0; n < My_; ++n)
        {
            const int q = m*My_ + n;
            if (q == 0 || !keep_[q]) continue;
            const scalar ky = (My_ > 1) ? n*M_PI/ylen_ : 0.0;
            const scalar k = std::sqrt(kx*kx + ky*ky);
            const scalar kz = k*zs;
            if (kz > 700) continue;
            const scalar E = std::exp(-kz);
            const scalar cc = cx[m]*cy[n];
            u -= kx*A_[q]*sx[m]*cy[n]*E;
            v -= ky*A_[q]*cx[m]*sy[n]*E;
            w -= k*A_[q]*cc*E;
            phi += A_[q]*cc*E;
            phit += At_[q]*cc*E;
        }
    }
}

vector Air2AirEnhanced::velocity(scalar x, scalar y, scalar z) const
{
    scalar u, v, w, phi, phit;
    evaluate(x/h_, y/h_, z/h_, u, v, w, phi, phit);
    const scalar U = std::sqrt(g_*h_);
    return vector(u*U, v*U, w*U);
}

scalar Air2AirEnhanced::pressure(scalar x, scalar y, scalar z) const
{
    scalar u, v, w, phi, phit;
    evaluate(x/h_, y/h_, z/h_, u, v, w, phi, phit);
    scalar p = -rhoAir_*g_*h_*phit;
    if (bernoulli_) p -= 0.5*rhoAir_*g_*h_*(u*u + v*v + w*w);
    return p;
}

scalar Air2AirEnhanced::potential(scalar x, scalar y, scalar z) const
{
    scalar u, v, w, phi, phit;
    evaluate(x/h_, y/h_, z/h_, u, v, w, phi, phit);
    return phi*h_*std::sqrt(g_*h_);
}

void Air2AirEnhanced::flow(scalar x, scalar y, scalar z, vector& U, scalar& Pd) const
{
    scalar u, v, w, phi, phit;
    evaluate(x/h_, y/h_, z/h_, u, v, w, phi, phit);
    const scalar Us = std::sqrt(g_*h_);
    U = vector(u*Us, v*Us, w*Us);
    Pd = -rhoAir_*g_*h_*phit;
    if (bernoulli_) Pd -= 0.5*rhoAir_*g_*h_*(u*u + v*v + w*w);
}

} // End namespace Foam
