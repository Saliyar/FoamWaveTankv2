/*---------------------------------------------------------------------------*\
  Air2AirEnhanced: analytical air field (orders 0-3) from HOS-NWT modes

  Air potential (non-dimensional HOS units: length h, time sqrt(h/g)):

      psi = w0 z + sum_(m,n)!=(0,0) A_mn cos(kx_m x) cos(ky_n y) exp(-k_mn z)

  z is measured from the MEAN water level (z - seaLevel), as in the
  derivation. A_mn and w0 come from modeFS (eta) and modeFSt (eta_t) by the
  perturbation recipe of air_coefficients_third_order.pdf, Section 15.7,
  implemented in analyticAirCoeffs.H (orders 1, 2, 3; no least squares).
  order 0 gives zero air velocity (sensitivity test only).

  Spectral cut-off (applied to eta and eta_t before the solve):
      cutoff auto        k_cut = cutoffC / max|eta|, max over ALL records
      cutoff kp          k_cut = cutoffKp * k_p (peak of time-averaged |eta_k|^2)
      cutoff wavelength  k_cut = 2 pi / cutoffLength [m]
      cutoff none        all modes
  Coefficients are solved once per HOS record (two-slot cache), interpolated
  linearly in time; dA/dt is the record difference over dt_HOS.

  This class reads the modes file itself (ASCII, fixed 18-character fields,
  records of 18*n1 characters, no newlines). It does not use HOSModesReader.

  Queries take wave-frame x, y [m] and z - seaLevel [m]; they return
  velocity [m/s], pressure [Pa] and potential [m^2/s].
\*---------------------------------------------------------------------------*/

#ifndef Air2AirEnhanced_H
#define Air2AirEnhanced_H

#include "dictionary.H"
#include "vector.H"

#include <vector>
#include <string>
#include <fstream>

namespace Foam
{

class Air2AirEnhanced
{
    struct slot
    {
        int index;                  // HOS record, -1 = empty
        scalar w0;
        std::vector<scalar> A;      // Mx*My, (m,n) -> m*My + n
    };

    // File
    std::string file_;
    mutable std::ifstream in_;
    int n1_, n2_, n3add_, nt_;
    long reclen_;
    scalar dt_, xlen_, ylen_, h_;           // dt, xlen, ylen non-dimensional; h [m]

    // Options
    scalar g_, rhoAir_;
    int order_;
    bool bernoulli_, verbose_;
    std::string cutoff_;

    // Cut-off
    scalar kcut_, etaMax_;
    int Mx_, My_;                            // retained block
    std::vector<char> keep_;                 // Mx*My mask (k <= kcut)

    // Cache and current state (non-dimensional)
    mutable slot cache_[2];
    mutable bool printed_;
    scalar w0_, w0t_;
    std::vector<scalar> A_, At_;
    bool updated_;

    // Helpers
    void readRecord(long rec, int n, std::vector<scalar>& v) const;
    void readModes(int it, std::vector<scalar>& e, std::vector<scalar>& et, int Mx, int My) const;
    void setupCutoff(const dictionary& dict);
    const slot& coeffs(int it, int other) const;
    void evaluate(scalar xs, scalar ys, scalar zs,
                  scalar& u, scalar& v, scalar& w, scalar& phi, scalar& phit) const;

    Air2AirEnhanced(const Air2AirEnhanced&);
    void operator=(const Air2AirEnhanced&);

public:

    Air2AirEnhanced(const string& modesFile, const dictionary& dict);
    ~Air2AirEnhanced();

    //- Set the coefficients for CFD time t [s]
    void update(scalar t);

    //- Velocity [m/s]; z measured from the mean water level [m]
    vector velocity(scalar x, scalar y, scalar z) const;

    //- Air pressure [Pa]: -rho psi_t (linear) or -rho(psi_t + |u|^2/2) (bernoulli)
    scalar pressure(scalar x, scalar y, scalar z) const;

    //- Velocity potential [m^2/s]
    scalar potential(scalar x, scalar y, scalar z) const;

    void flow(scalar x, scalar y, scalar z, vector& U, scalar& Pd) const;

    int getNt() const { return nt_; }
    scalar getDt() const { return dt_*std::sqrt(h_/g_); }   // [s]
    int getOrder() const { return order_; }
    scalar kCut() const { return kcut_/h_; }                 // [1/m]
    scalar h() const { return h_; }
};

} // End namespace Foam

#endif
