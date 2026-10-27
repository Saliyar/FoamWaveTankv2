// Phase 3 air-side checks (standalone, no OpenFOAM)
// Uses the same algorithms as Air2AirEnhanced (analyticAirCoeffs.H).
#include "analyticAirCoeffs.H"
#include <cstdio>
#include <cmath>

using namespace airAnalytic;

// Stokes-like steady wave in the nondimensional frame (k = 1, c = omega):
// eta(x,t) = sum_m a_m cos(m (x - c t)) + b_m sin(m (x - c t)), period 2 pi.
// Wave speed c = omega = 1 (nondim), L = 2 pi, k_m = m.
static const int M = 32;
static const double L = 2.0*M_PI;

static void etaCoeffs(double t, double c, vec& ec, vec& es, vec& etc, vec& ets)
{
    // base profile with cos and sin harmonics (decaying amplitudes)
    vec a(M, 0.0), b(M, 0.0);
    a[1] = 0.1; b[1] = 0.0;
    a[2] = 0.004; b[2] = 0.0;
    a[3] = 0.0008; b[3] = 0.0003;
    ec.assign(M, 0.0); es.assign(M, 0.0);
    etc.assign(M, 0.0); ets.assign(M, 0.0);
    for (int m = 0; m < M; ++m)
    {
        const double th = m*c*t;
        // shifted profile at t: cos(m(x - c t)) etc. from base cos/sin
        ec[m] = a[m]*std::cos(th) - b[m]*std::sin(th);
        es[m] = a[m]*std::sin(th) + b[m]*std::cos(th);
    }
    // eta_t = -c eta_x exactly: coefficients -c*m*(sin part), +c*m*(cos part)
    for (int m = 1; m < M; ++m)
    {
        const double w = c*m;
        etc[m] = -w*es[m];
        ets[m] = w*ec[m];
    }
}

int main()
{
    const double c = 1.0;
    const int order = 3;

    // ---- Check 1: periodic vs cosine basis on an even eta ----
    // Cosine code on [0, Lx], Lx = L/2, with modes m pi/Lx = m (since Lx = pi)
    {
        const int Mc = 16;
        vec e(Mc, 0.0), et(Mc, 0.0);
        e[1] = 0.1; e[2] = 0.004; e[3] = 0.0008;
        et[1] = -c*0.0; // at t = 0, eta_t = -c eta_x; even eta gives odd eta_t -> use
        // even eta_t as a test input (cos series): et = c*k*(sin part); use cos-only test
        et[1] = 0.05; et[2] = 0.002; et[3] = -0.0004;
        vec Ac; double w0c;
        analyticAirCoeffs(Mc, 1, M_PI, 1.0, e, et, order, Ac, w0c);

        vec ec(M, 0.0), es(M, 0.0), etc(M, 0.0), ets(M, 0.0);
        for (int m = 0; m < Mc; ++m) { ec[m] = e[m]; etc[m] = et[m]; }
        vec Ap, Bp; double w0p;
        analyticAirCoeffsPeriodic(M, L, ec, es, etc, ets, order, Ap, Bp, w0p);

        double maxd = 0.0, maxA = 0.0, maxB = 0.0;
        for (int m = 0; m < Mc; ++m)
        {
            maxd = std::max(maxd, std::fabs(Ap[m] - Ac[m]));
            maxA = std::max(maxA, std::fabs(Ac[m]));
            maxB = std::max(maxB, std::fabs(Bp[m]));
        }
        std::printf("CHECK1 periodic vs cosine (even input): max|dA| = %.3e, max|A| = %.3e, max|B| = %.3e\n",
                    maxd, maxA, maxB);
    }

    // ---- Check 2: A_t exact (rotation) vs central difference ----
    // Solve A(t) from eta(t) with the periodic algorithm at each t;
    // compare the analytic rotation rule with a central difference
    {
        const double T = 2.0*M_PI/c;   // period in time units
        const double tau = T/1000.0;
        for (int order_ : {1, 2, 3})
        {
            const double t0 = 0.37;

            auto solveAt = [&](double t, vec& Ac, vec& As, double& w0)
            {
                vec ec, es, etc, ets;
                etaCoeffs(t, c, ec, es, etc, ets);
                analyticAirCoeffsPeriodic(M, L, ec, es, etc, ets, order_, Ac, As, w0);
            };

            vec Ac0, As0; double w00;
            solveAt(t0, Ac0, As0, w00);

            // rotation rule: omega_m = c m, At_c = -omega B(t), At_s = +omega A(t)
            vec Acp, Asp; double w0p;
            solveAt(t0 + tau, Acp, Asp, w0p);
            vec Acm, Asm; double w0m;
            solveAt(t0 - tau, Acm, Asm, w0m);

            double maxRel = 0.0, maxAbs = 0.0;
            for (int m = 1; m < M; ++m)
            {
                const double om = c*m;
                const double rotAt = -om*As0[m];
                const double rotBt = om*Ac0[m];
                const double cdAt = (Acp[m] - Acm[m])/(2*tau);
                const double cdBt = (Asp[m] - Asm[m])/(2*tau);
                const double scale = std::max(1e-30, std::fabs(Ac0[m]) + std::fabs(As0[m]));
                maxAbs = std::max(maxAbs, std::max(std::fabs(rotAt - cdAt), std::fabs(rotBt - cdBt)));
                maxRel = std::max(maxRel, std::max(std::fabs(rotAt - cdAt), std::fabs(rotBt - cdBt))/scale);
            }
            std::printf("CHECK2 order %d: max|At_rot - At_cd| = %.3e, relative to |A| = %.3e\n",
                        order_, maxAbs, maxRel);
        }
    }

    // ---- Check 3: air kinematic mismatch at z = eta, orders 1..3 ----
    // psi = sum (A cos + B sin) exp(-k z) (+ w0 z = 0): w = d psi/dz,
    // u = d psi/dx. At z = eta(x): w - (eta_t + u eta_x) should decrease with order.
    {
        double mism[4] = {0, 0, 0, 0};
        for (int order_ = 1; order_ <= 3; ++order_)
        {
            vec ec, es, etc, ets;
            etaCoeffs(0.0, c, ec, es, etc, ets);
            vec Ac, As; double w0;
            analyticAirCoeffsPeriodic(M, L, ec, es, etc, ets, order_, Ac, As, w0);

            double maxm = 0.0;
            const int Nx = 400;
            for (int i = 0; i < Nx; ++i)
            {
                const double x = (i + 0.5)*L/Nx;
                double eta = 0, etx = 0, etam = 0;
                for (int m = 0; m < M; ++m)
                {
                    const double k = m;
                    eta += ec[m]*std::cos(k*x) + es[m]*std::sin(k*x);
                    etx += -k*ec[m]*std::sin(k*x) + k*es[m]*std::cos(k*x);
                    etam += etc[m]*std::cos(k*x) + ets[m]*std::sin(k*x);
                }
                // air velocities at z = eta (nondimensional), all modes
                double w = w0, u = 0;
                for (int m = 1; m < M; ++m)
                {
                    const double k = m;
                    const double E = std::exp(-k*eta);
                    const double cs = std::cos(k*x), sn = std::sin(k*x);
                    // psi_z = -k (A cos + B sin) E ; psi_x = k (B cos - A sin) E
                    w += -k*(Ac[m]*cs + As[m]*sn)*E;
                    u += k*(As[m]*cs - Ac[m]*sn)*E;
                }
                const double kin = w - (etam + u*etx);
                maxm = std::max(maxm, std::fabs(kin));
            }
            mism[order_] = maxm;
            std::printf("CHECK3 order %d: max |w - (eta_t + u eta_x)| at z = eta = %.3e\n", order_, maxm);
        }
        std::printf("CHECK3 ratio order2/order3 = %.3f\n", mism[2]/mism[3]);
    }

    return 0;
}
