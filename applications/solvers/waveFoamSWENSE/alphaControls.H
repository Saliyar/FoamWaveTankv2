// Alpha (VoF) controls, read from the alpha solver dictionary of the case
const dictionary& alphaControls = mesh.solverDict(alpha1.name());

const label nAlphaCorr
(
    readLabel(alphaControls.lookup("nAlphaCorr"))
);

const label nAlphaSubCycles
(
    alphaControls.lookupOrDefault<label>("nAlphaSubCycles", 1)
);

const bool MULESCorr
(
    alphaControls.lookupOrDefault<Switch>("MULESCorr", false)
);

const bool alphaApplyPrevCorr
(
    alphaControls.lookupOrDefault<Switch>("alphaApplyPrevCorr", false)
);

const scalar icAlpha
(
    alphaControls.lookupOrDefault<scalar>("icAlpha", 0.0)
);

const scalar scAlpha
(
    alphaControls.lookupOrDefault<scalar>("scAlpha", 0.0)
);

const bool LTS = false;
