# Pre- and post-processing utilities for dissolFoam

1) genBlockMeshDict: Python utility to generate blockMeshDict files.
<br>usage:
<br>    python /PATH/TO/YOUR/genBlockMeshDict/runGen.py -l              # List possible generators
<br>    python /PATH/TO/YOUR/genBlockMeshDict/runGen.py -g generator    # Create dictionary
<br>    python /PATH/TO/YOUR/genBlockMeshDict/runGen.py -d genDict      # Create blockMeshDict
<br>The case dissolCirc uses genBLockMeshDict in the script makeSTL


2) surfRoughGen: OpenFOAM utility to generate a rough surface by Fourier synthesis.
<br>usage:
<br>    surfRoughGen
<br>It has a dictionary system/surfRoughGenDict - the case dissolFrac uses surfRoughGen<\br>


3) orderBoundaries: OpenFoam utility to reorder boundary patches for normalMotionSlip
<br>usage:
<br>    orderBoundaries -overwrite
<br>The case dissolCirc uses orderBoundaries<\br>


4) dissolCalc: OpenFOAM utility to postprocess fracture fields by integrating over the aperture
<br>usage:
<br>    dissolCalc fieldMap2D all 1000 100
<br>Calculates all the 2D fields using 1000 cells in flow direction and 100 cells in transverse direction. It uses a dictionary fieldMap2Ddict


