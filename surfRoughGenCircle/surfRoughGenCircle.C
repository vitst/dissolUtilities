/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-201X OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
  surfRoughGenCircle

Description
  Preprocessing utility to modify the fracture surface represented by
  two parallel plates. It overwrites 0 directory with modified surface.

  Limitations now are the next:
    - a geometry should include one or two flat plates

Usage
  - surfRoughGenCircle

Needs dictionary
  system/surfRoughGenDict.circle

\*---------------------------------------------------------------------------*/

#include <cmath>
#include <complex>
#include <fftw3.h>
#include <time.h>
#include <sys/types.h>
#include <iostream>
#include <unistd.h>
//#include <vector>

#include "coupledPatchInterpolation.H"
#include "dynamicFvMesh.H"
#include "fvCFD.H"
#include "pointPatchField.H"

#include "fixedValuePointPatchField.H"
#include "normalMotionSlipPointPatchVectorField.H"

#include "velocityMotionSolver.H"

#include "motionSolver.H"

/*
 #######################################################################################
 *    Main program body
 #######################################################################################
*/

class RoughnessGenerator 
{
protected:
    // Protected data
    unsigned long seed;
    int totalNum;
    vector center;
    vector direction;
    scalar radius;
    double rgh;
    double cutoff;
    double dHurst; // Fractal dimension D = 3 - dHurst
    word wayToApply;

public:
    // Constructor
    RoughnessGenerator(unsigned long seed_, int totalNum_, point center_, vector direction_,
            scalar radius_, scalar rgh_, scalar cutoff_, scalar dHurst_, word wayToApply_)
        : seed(seed_), totalNum(totalNum_), center(center_), direction(direction_),
            radius(radius_), rgh(rgh_), cutoff(cutoff_), dHurst(dHurst_), 
            wayToApply(wayToApply_)
    {}

    static unsigned long mix(unsigned long a, unsigned long b, unsigned long c)
    {
        a=a-b;  a=a-c;  a=a^(c >> 13);
        b=b-c;  b=b-a;  b=b^(a << 8);
        c=c-a;  c=c-b;  c=c^(b >> 13);
        a=a-b;  a=a-c;  a=a^(c >> 12);
        b=b-c;  b=b-a;  b=b^(a << 16);
        c=c-a;  c=c-b;  c=c^(b >> 5);
        a=a-b;  a=a-c;  a=a^(c >> 3);
        b=b-c;  b=b-a;  b=b^(a << 10);
        c=c-a;  c=c-b;  c=c^(b >> 15);
        return c;
    }
  
  
    void getSurfaceDisplacement(dynamicFvMesh &mesh, scalarField &wd, label &patchID) 
    {
        scalarField sFp(totalNum, 0.0);
        fftDisp(sFp);
    
        pointField pointFace = mesh.boundaryMesh()[patchID].faceCentres();
        pointField normlFace = mesh.boundaryMesh()[patchID].faceNormals();
    
        // scalar cL = constant::mathematical::pi * ;
        scalarField phi(totalNum, 0.0);
        forAll(pointFace, i)
        {
            vector radVec = pointFace[i] - center;
      
            scalar phi = Foam::atan2(radVec[1], radVec[0]);
            scalar phiDeg = phi * 180.0 / constant::mathematical::pi;
            if (phiDeg < 0)
              phiDeg += 360.0;
      
            int ind0 = std::trunc(phiDeg/360.0 * totalNum);
            int ind1 = std::ceil(phiDeg/360.0 * totalNum);
            if(ind1>totalNum) ind1 = totalNum;
      
            scalar x0 = radius * Foam::cos(ind0*360.0/totalNum);
            scalar y0 = radius * Foam::sin(ind0*360.0/totalNum);
            scalar x1 = radius * Foam::cos(ind1*360.0/totalNum);
            scalar y1 = radius * Foam::sin(ind1*360.0/totalNum);
      
            point p0(x0, y0, 0.05), p1(x1, y1, 0.05);
      
            scalar w0 = 1.0 / mag( pointFace[i] - p0);
            scalar w1 = 1.0 / mag( pointFace[i] - p1);
      
      
            wd[i] = (sFp[ind0] * w0 + sFp[ind1] * w1) / (w0 + w1);
        }
    
        Info << "Displacement calculated" << nl;
    }
  
private:
    scalar power(double ksq) 
    {
      if (ksq == 0)
        return 0; // <rad^2> ~ 1/ksq^(1+H)
      if (ksq > 1)
        return 0; // cutoff wavelength = cutLen
      scalar p = Foam::pow(ksq, -(dHurst + 1));
      //    p *= Foam::exp(-ksq);
      return std::sqrt(p);
    }
  
    void fftDisp(scalarField &disp)
    {
        scalar TwoPi = constant::mathematical::twoPi;
        scalar L = TwoPi * radius; 
        scalar cutoffL = L / cutoff;
    
        std::vector<std::complex<double>> f, F;
        f.resize(totalNum);
        F.resize(totalNum);
    
        Random rnd(seed);
    
        Info << "Displacement calc starts...." << nl;
    
        for (int n = 0; n < totalNum / 2 + 1; ++n)
        {
            scalar p = TwoPi * rnd.sample01<scalar>();
            double k = n * cutoffL / L;
            double ksq = k * k;
            scalar rad = power(ksq) * rnd.GaussNormal<scalar>();
      
            f[n] = rad * std::complex<double>(Foam::cos(p), Foam::sin(p));
            f[((totalNum - n) % totalNum)] =
                rad * std::complex<double>(Foam::cos(p), -Foam::sin(p));
        }
    
        f[totalNum / 2].imag(0.0);
    
        fftw_plan plan;
        plan = fftw_plan_dft_1d( totalNum, 
                                 reinterpret_cast<fftw_complex *>(&f[0]),
                                 reinterpret_cast<fftw_complex *>(&F[0]),
                                 FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);
    
        scalarField sF(totalNum);
        forAll(sF, ii) { sF[ii] = F[ii].real(); }
    
        scalarField sF2 = sqr(sF);
        scalar avSF = average(sF);
        scalar avSF2 = average(sF2);
    
        scalar factor = rgh / Foam::sqrt(mag(avSF2 - sqr(avSF)));
    
        sF *= factor;
    
        disp = sF;
    }
};

int main(int argc, char *argv[]) 
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
  
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    // // reading dictionary surfRoughGenDict.circle
    IOdictionary surfRoughGenDict(IOobject("surfRoughGenDict.circle", runTime.system(),
                                           mesh, IOobject::MUST_READ,
                                           IOobject::NO_WRITE));
  
    word wayToApply;
    if (!surfRoughGenDict.readIfPresent<word>("apply", wayToApply)) {
      SeriousErrorIn("main")
          << "There is no `apply` parameter in dictionary"
          << exit(FatalError);
    }
    word patchName;
    if (!surfRoughGenDict.readIfPresent<word>("patchName", patchName)) {
      SeriousErrorIn("main") << "There is no `patchName` parameter in dictionary"
                             << exit(FatalError);
    }
  
    vector center;
    if (!surfRoughGenDict.readIfPresent<vector>("center", center)) {
      SeriousErrorIn("main") << "There is no `center` parameter in dictionary"
                             << exit(FatalError);
    }
    vector direction;
    if (!surfRoughGenDict.readIfPresent<vector>("direction", direction)) {
      SeriousErrorIn("main") << "There is no `direction` parameter in dictionary"
                             << exit(FatalError);
    }
    scalar radius;
    if (!surfRoughGenDict.readIfPresent<scalar>("radius", radius)) {
      SeriousErrorIn("main") << "There is no `radius` parameter in dictionary"
                             << exit(FatalError);
    }
    int totalNum;
    if (!surfRoughGenDict.readIfPresent<int>("totalNum", totalNum)) {
      SeriousErrorIn("main") << "There is no `totalNum` parameter in dictionary"
                             << exit(FatalError);
    }
  
    unsigned long seed;
    if (!surfRoughGenDict.readIfPresent<unsigned long>("seed", seed)) {
      SeriousErrorIn("main") << "There is no `seed` parameter in dictionary"
                             << exit(FatalError);
    }
    if(seed==0)
        seed = RoughnessGenerator::mix(std::clock(), time(NULL), getpid());
  
    scalar rgh;
    if (!surfRoughGenDict.readIfPresent<scalar>("roughness", rgh)) {
      SeriousErrorIn("main") << "There is no `roughness` parameter in dictionary"
                             << exit(FatalError);
    }
    double dHurst;
    if (!surfRoughGenDict.readIfPresent<double>("dHurst", dHurst)) {
      SeriousErrorIn("main") << "There is no `dHurst` parameter in dictionary"
                             << exit(FatalError);
    }
    double cutoff;
    if (!surfRoughGenDict.readIfPresent<double>("cutoff", cutoff)) {
      SeriousErrorIn("main") << "There is no `cutoff` parameter in dictionary"
                             << exit(FatalError);
    }
  
    Info << "patch:         " << patchName << endl;
    Info << "apply (method):" << wayToApply << endl;
    Info << "center:        " << center << endl;
    Info << "direction:     " << direction << endl;
    Info << "radius:        " << radius << endl;
    Info << "totalNum:      " << totalNum << endl;
    Info << "seed:          " << seed << endl;
    Info << "roughness:     " << rgh << endl;
    Info << "dHurst:        " << dHurst << endl;
    Info << "cutoff:        " << cutoff << endl;
    Info << "Setup RoughnessGenerator class" << endl;
  
    RoughnessGenerator rg(seed, totalNum, center, direction,
            radius, rgh, cutoff, dHurst, wayToApply);
  
    double cpuTime = runTime.elapsedCpuTime();
  
    // Get patch ID for moving boundaries
    label patchID = mesh.boundaryMesh().findPatchID(patchName);
    if (patchID == -1) {
        SeriousErrorIn("main") << "patch " << patchName << " is missing"
                               << exit(FatalError);
    }
  
    coupledPatchInterpolation patchInterpolator(mesh.boundaryMesh()[patchID],
                                                mesh);
  
    const pointField &boundaryPoints = mesh.boundaryMesh()[patchID].localPoints();
    vectorField pointDispWall(boundaryPoints.size(), vector::zero);
    vectorField pointNface = mesh.boundaryMesh()[patchID].faceNormals();
    vectorField motionN = patchInterpolator.faceToPointInterpolate(pointNface);
    forAll(motionN, ii) motionN[ii] /= mag(motionN[ii]);
  
    scalarField faceDisp(pointNface.size(), 0.0);
    rg.getSurfaceDisplacement(mesh, faceDisp, patchID);
    scalarField pointDisp = patchInterpolator.faceToPointInterpolate(faceDisp);
  
    /*
    forAll(pointDisp, i) {
      pointDisp[i] = std::min(pointDisp[i],  maxDisp);
      pointDisp[i] = std::max(pointDisp[i], -maxDisp);
    }
    */
  
    Info << "Maximum and minimum face displacements  " << max(faceDisp) << "  "
         << min(faceDisp) << endl;
    Info << "Maximum and minimum point displacements  " << max(pointDisp) << "  "
         << min(pointDisp) << endl;
  
    forAll(pointDispWall, i) pointDispWall[i] = pointDisp[i] * motionN[i];
  
    pointVectorField &pointVelocity = const_cast<pointVectorField &>(
        mesh.objectRegistry::lookupObject<pointVectorField>("pointMotionU"));
  
    pointVelocity.boundaryFieldRef()[patchID] == pointDispWall;
    mesh.update();
  
    cpuTime = runTime.elapsedCpuTime() - cpuTime;
  
    Info << nl << "Time statistics:" << nl;
  
    int wlNP = mesh.boundaryMesh()[patchID].nPoints();
    Info << "Total number of points:                   " << mesh.nPoints() << nl;
    Info << "Number of points on the walls:            " << wlNP << nl;
    Info << "Running time:                             " << cpuTime << nl << endl;
  
    Info << "Overwriting points in current time directory." << nl;
    runTime.writeNow();
  
    Info << "End" << nl;
    return 0;
}

// ************************** End of the solver ************************** //
