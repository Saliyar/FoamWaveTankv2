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

    cache_[0].index = -1;
    cache_[1].index = -1;

    Info<< "Air2AirEnhanced: order " << order_
        << ", pressureModel " << (bernoulli_ ? "bernoulli" : "linear") << nl
        << "    modes file " << file_ << nl
        << "    n1 = " << n1_ << ", n2 = " << n2_ << ", records = " << nt_
        << ", dt_HOS = " << getDt() << " s, h = " << h_ << " m" << nl
        << "    tank " << xlen_*h_ << " x " << ylen_*h_ << " m" << endl;
}

Air2AirEnhanced::~Air2AirEnhanced()
{}

// * * * * * * * * * * * * * * Coefficients  * * * * * * * * * * * * * * * * //

// Coefficients of record it; 'other' is the record that must stay cached
const Air2AirEnhanced::slot& Air2AirEnhanced::coeffs(int it, int other) const
{
    for (int s = 0; s < 2; ++s)
        if (cache_[s].index == it) return cache_[s];

    const int s = (cache_[0].index == other) ? 1 : 0;
    slot& c = cache_[s];

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
    scalar s = std::max(scalar(0), t/Ts)/dt_;          // position in records
    s = std::min(s, scalar(nt_ - 1));
    const int i0 = std::min(int(std::floor(s)), nt_ - 2);
    const scalar f = s - i0;

    const slot& c0 = coeffs(i0, i0 + 1);
    const slot& c1 = coeffs(i0 + 1, i0);

    const int N = Mx_*My_;
    A_.resize(N);
    At_.resize(N);
    for (int q = 0; q < N; ++q)
    {
        A_[q] = (1 - f)*c0.A[q] + f*c1.A[q];
        At_[q] = (c1.A[q] - c0.A[q])/dt_;
    }
    w0_ = (1 - f)*c0.w0 + f*c1.w0;
    w0t_ = (c1.w0 - c0.w0)/dt_;
    updated_ = true;
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
