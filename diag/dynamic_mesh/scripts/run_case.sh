#!/bin/bash
# usage: run_case.sh LABEL BINARY SRCCASE MODE NSTEPS SCRATCH
# MODE: prod   production solver, case unchanged
#       static _dyn solver, dynamicMeshDict staticFvMesh
#       zero   _dyn solver, dynamicMotionSolverFvMesh with zero displacement
# Runs NSTEPS time steps (endTime = NSTEPS deltaT) in SCRATCH/LABEL.

LABEL=$1
BIN=$2
SRC=$3
MODE=$4
N=$5
SCRATCH=$6

source /home/projects/dtu_00075/people/sitali/SWENSE_testing/diag/box_fixed/env.sh > /dev/null 2>&1

OUT=$SCRATCH/$LABEL
rm -rf $OUT
mkdir -p $OUT
cp -r $SRC/0 $SRC/constant $SRC/system $OUT/
cd $OUT || exit 1

DT=$(awk '$1=="deltaT"{print $2}' system/controlDict | tr -d ';')
END=$(awk -v n=$N -v d=$DT 'BEGIN{printf "%.10g", n*d}')
sed -i "s/^endTime .*/endTime         $END;/" system/controlDict
sed -i "s/^writeControl .*/writeControl    timeStep;/" system/controlDict
sed -i "s/^writeInterval .*/writeInterval   ${WRITE_EVERY:-$N};/" system/controlDict

# Tutorials without a mesh: blockMesh as their Allrun does; the HOS modes file
# (not in the repository) is linked from MODES when it is missing
if [ ! -f constant/polyMesh/points ] && [ -f system/blockMeshDict ]
then
    blockMesh > log.blockMesh 2>&1
fi
if [ ! -f constant/modes_HOS_SWENSE_clean.dat ] && [ -n "$MODES" ]
then
    ln -s $MODES constant/modes_HOS_SWENSE_clean.dat
fi
# Coarse tutorials predate the psiProj and cAlpha entries of the current
# solvers: their fvSolution and fvSchemes are replaced by the box case ones
if [ -n "$FVSOL" ]
then
    cp $FVSOL system/fvSolution
fi
if [ -n "$FVSCHEMES" ]
then
    cp $FVSCHEMES system/fvSchemes
fi

if [ "$MODE" = "sixdof" ]
then
    # Free floating body: OpenFOAM v2506 sixDoFRigidBodyMotion under
    # dynamicMotionSolverFvMesh (box geometry measured from the mesh, see report)
    cat > constant/dynamicMeshDict <<EOF
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      dynamicMeshDict;
}

dynamicFvMesh       dynamicMotionSolverFvMesh;

motionSolverLibs    (sixDoFRigidBodyMotion);

motionSolver        sixDoFRigidBodyMotion;

updateControl       timeStep;
updateInterval      1;

patches             (boxFront boxBack boxBottom);
innerDistance       0.3;
outerDistance       1.0;

mass                509.024;
centreOfMass        (0 0 0.9033);
momentOfInertia     (1549.36 213.33 1354.58);
orientation         (1 0 0 0 1 0 0 0 1);
angularMomentum     (0 0 0);
g                   (0 0 -9.81);
rho                 rho;
report              on;

solver
{
    type            symplectic;
}

constraints
{
    plane
    {
        sixDoFRigidBodyMotionConstraint plane;
        normal          (0 1 0);
        point           (0 0 0);
    }
}

EOF
fi

if [ "$MODE" = "translate" ]; then cp /gpfs/hpc-home/projects/dtu_00075/apps/FoamWaveTank/wt_SWENSE_dynamic/diag/dynamic_mesh/scripts/translate.dynamicMeshDict.template constant/dynamicMeshDict; fi
if [ "$MODE" = "sixdof" ]; then cp /gpfs/hpc-home/projects/dtu_00075/apps/FoamWaveTank/wt_SWENSE_dynamic/diag/dynamic_mesh/scripts/pointDisplacement.template 0/pointDisplacement; fi

if [ "$MODE" = "static" ] || [ "$MODE" = "zero" ]
then
    if [ "$MODE" = "static" ]
    then
        DYN=staticFvMesh
    else
        DYN=dynamicMotionSolverFvMesh
    fi
    cat > constant/dynamicMeshDict <<EOF
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      dynamicMeshDict;
}

dynamicFvMesh       $DYN;

EOF
    if [ "$MODE" = "zero" ]
    then
        cat >> constant/dynamicMeshDict <<EOF
motionSolverLibs    (fvMotionSolvers);
motionSolver        displacementLaplacian;
diffusivity         uniform;
EOF
        cat > 0/pointDisplacement <<EOF
FoamFile
{
    version     2.0;
    format      ascii;
    class       pointVectorField;
    object      pointDisplacement;
}

dimensions      [0 1 0 0 0 0 0];
internalField   uniform (0 0 0);

boundaryField
{
    ".*"
    {
        type            fixedValue;
        value           uniform (0 0 0);
    }
}
EOF
    fi
fi

{
    echo "label: $LABEL"
    echo "mode: $MODE  steps: $N  deltaT: $DT  endTime: $END"
    echo "binary realpath: $(realpath $BIN)"
    echo "binary md5: $(md5sum $BIN | cut -d' ' -f1)"
    echo "source case: $SRC"
} > log.run

START=$(date +%s)
$BIN >> log.run 2>&1
RC=$?
END_T=$(date +%s)
echo "exit code: $RC" >> log.run
echo "wall time s: $((END_T - START))" >> log.run
echo "RUN_DONE $LABEL rc=$RC wall=$((END_T - START))"
