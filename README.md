# Pre- and post-processing utilities for dissolFoam

1) genBlockMeshDict: Python utility to generate blockMeshDict files.
usage:
    python /PATH/TO/YOUR/genBlockMeshDict/runGen.py -l              # List possible generators
    python /PATH/TO/YOUR/genBlockMeshDict/runGen.py -g generator    # Create dictionary
    python /PATH/TO/YOUR/genBlockMeshDict/runGen.py -d genDict      # Create blockMeshDict
The case dissolCirc uses genBLockMeshDict in the script makeSTL

2) surfRoughGen: OpenFOAM utility to generate a rough surface by Fourier synthesis.
usage:
    surfRoughGen
It has a dictionary system/surfRoughGenDict - the case dissolFrac uses surfRoughGen

3) orderBoundaries: OpenFoam utility to reorder boundary patches for normalMotionSlip
usage:
    orderBoundaries -overwrite
The case dissolCirc uses orderBoundaries

4) dissolCalc: OpenFOAM utility to postprocess fracture fields by integrating over the aperture
usage:
    dissolCalc fieldMap2D all 1000 100
Calculates all the 2D fields using 1000 cells in flow direction and 100 cells in transverse direction. It uses a dictionary fieldMap2Ddict


