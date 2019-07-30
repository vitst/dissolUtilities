# Pre- and post-processing utilities for dissolFoam

1) surfRoughGen: OpenFOAM utility to generate a rough surface by Fourier synthesis.
<br>usage:
<br>    surfRoughGen
<br>It has a dictionary system/surfRoughGenDict - the case dissolFrac uses surfRoughGen<\br>

2) runMeshUpdateOnce: OpenFOAM utility for a single mesh update cycle (after surfRoughGen)
<br>usage:
<br>    runMeshUpdateOnce


3) dissolCalc: OpenFOAM utility to postprocess fracture fields by integrating over the aperture
<br>usage:
<br>    dissolCalc fieldMap2D all 1000 100
<br>Calculates all the 2D fields using 1000 cells in flow direction and 100 cells in transverse direction. It uses a dictionary fieldMap2Ddict


