// isoAdvector-specific fields and controls for the alpha step
#include "alphaControls.H"
isoAdvection advector(alpha1, phi, U);

// isoAdvector for alphaI with the incident flux (same settings as alpha1).
// phiIRel equals the projected incident flux phiInc on the static mesh
isoAdvection advectorI(alphaI, phiIRel, incident.UInc());
