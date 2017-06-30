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
  surfRoughGen

Description
  Preprocessing utility to modify the fracture surface represented by
  two parallel plates. It overwrites 0 directory with modified surface.
  
  Limitations now are the next:
    - a geometry should include two parallel plates
    - the number of faces on each surface should be 
      a power of 2 (for Fourier transform)

Usage
  - surfRoughGen

Needs dictionary
  system/surfRoughGenDict
    
\*---------------------------------------------------------------------------*/

#include <complex>
#include <fftw3.h>

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "coupledPatchInterpolation.H"
#include "pointPatchField.H"
//#include "meshRelax.H"

/*
 #######################################################################################
 *    Main program body
 #######################################################################################
*/

class RoughnessGenerator
{
protected:
  // Protected data
  int seed;
  int M;
  int N;
  double rgh;
  word wayToApply;
  double fractalParam; // Set Hausdorff dimension H = 3 - D_f
  double smoothing;
  
public:
  
  // Constructor
  RoughnessGenerator
  (
    int seed_,
    int M_,
    int N_,
    double rgh_,
    word wayToApply_,
    double fractalParam_,
    double smoothing_
  )
  :
    seed(seed_),
    M(M_),
    N(N_),
    rgh(rgh_),
    wayToApply(wayToApply_),
    fractalParam(fractalParam_),
    smoothing(smoothing_)
  {
  }

  void getSurfaceDisplacement
  (
    dynamicFvMesh& mesh,
    scalarField& wd,
    label& wallID,
    int majDir,
    int latDir
  )
  {
    scalarField sFn(M*N, 0.0);
    fftDisp(sFn);
    scalarField sFp(M*N, 0.0);
    if(wayToApply=="asymmetric")
    {
      // shift the seed to get two different numbers
      seed += 125522;
      fftDisp(sFp);
    }
    else
    {
      sFp = sFn;
    }

    pointField pointFace = mesh.boundaryMesh()[wallID].faceCentres();
    pointField normlFace = mesh.boundaryMesh()[wallID].faceNormals();
    scalar maxMaj = max( pointFace.component(majDir) );
    scalar maxLat = max( pointFace.component(latDir) );
    scalar minMaj = min( pointFace.component(majDir) );
    scalar minLat = min( pointFace.component(latDir) );
    double Llat = maxLat - minLat;
    double Lmaj = maxMaj - minMaj;
    
    int mveDir = 3 - (latDir + majDir);
    
    if(wayToApply=="internalCylinder")
    {
      // TODO this is just a quick solution, make it general
      scalar cL = constant::mathematical::pi * 198.0;
      forAll(pointFace, i)
      {
        double x,y;
        x = pointFace[i].x();
        y = pointFace[i].y();
        
        scalar phi = Foam::atan2(y,x);
        phi /= constant::mathematical::twoPi;
        if(phi<0) phi += 1.0;
        
        scalar curMaj = pointFace[i].z();
        
        if(curMaj<cL)
        {
          int curm = std::floor(curMaj / cL * (M-1));
          int curn = std::floor(phi * (N-1));
          int ind = curn + N * curm;
          wd[i] = sFn[ind] * ( 1 - curMaj / cL );
        }
      }
    }
    else
    {
      forAll(pointFace, i)
      {
        scalar curMaj = pointFace[i].component(majDir) - minMaj;
        scalar curLat = pointFace[i].component(latDir) - minLat;

        scalar sign = normlFace[i].component(mveDir) 
                / 
                mag(normlFace[i].component(mveDir));

        int curm = std::floor(curMaj / Lmaj * (M-1));
        int curn = std::floor(curLat / Llat * (N-1));
        int ind = curn + N * curm;

        if(wayToApply=="symmetric" || wayToApply=="oneSurface")
          wd[i] = sFn[ind];

        if(wayToApply=="oneSurfaceDecay")
        {
          //wd[i] = sFn[ind];
          double factor = ( 1 - curMaj / Llat );
          if( factor<0.0 ) factor = 0.0;
          wd[i] = sFn[ind] * factor;
        }
        
        if(wayToApply=="synchronous")
          wd[i] = sign * sFn[ind];

        if(wayToApply=="asymmetric")
        {
          if(sign<0)
            wd[i] = sFn[ind];
          else
            wd[i] = sFp[ind];
        }
      }
    }
    Info<<"Displacement calculated"<<nl;
  }

private:
  // converts indexes
  label index(label m, label n){ return n + N * m; }
  
  scalar pspec(int u)
  {
    scalar p = Foam::pow(u, -0.5 * (fractalParam+1) );
    p *= Foam::exp(-smoothing*u);
    return p;
  }

  void fftDisp(scalarField& disp)
  {
    unsigned int MN = M*N;

    if ( MN & (MN - 1) )
    {
        FatalErrorIn
        (
             "getSurfaceDisplacement  "
        )   << "number of elements is not a power of 2" << endl
            << "    Number of elements = " << MN
            << abort(FatalError);
    }

    std::vector<std::complex<double> > f, F;
    f.resize(MN);
    F.resize(MN);

    Random rnd( seed );
    scalar TwoPi = constant::mathematical::twoPi;

    Info<< "Displacement calc starts...."<<nl;
    /*
     *   --- ---
     *  | 1 | 2 |
     *   --- ---
     *  | 3 | 4 |
     *   --- ---
     */
    // calculating 1 and 4
    for(int m=0; m<M/2+1; ++m)
    {
      for(int n=0; n<N/2+1; ++n)
      {
        scalar p = TwoPi * rnd.scalar01();
        int u = N * m*m / static_cast<double>(M) + M * n*n / static_cast<double>(N);

        scalar rad = 0.0;
        if(u == 0)
          rad = 0.0;
        else
          rad = pspec(u) * rnd.GaussNormal();

        f[ index(m,n) ] = 
                rad * std::complex<double>(Foam::cos(p),  Foam::sin(p));
        f[ index(((M-m)%M),(N-n)%N) ] = 
                rad * std::complex<double>(Foam::cos(p), -Foam::sin(p));
      }
    }

    f[ index(M/2,0)   ].imag(0.0);
    f[ index(0,  N/2) ].imag(0.0);
    f[ index(M/2,N/2) ].imag(0.0);

    for(int m=1; m<M/2; ++m)
    {
      for(int n=1; n<N/2; ++n)
      {
        scalar p = TwoPi * rnd.scalar01();
        int u = N * m*m / static_cast<double>(M) + M * n*n / static_cast<double>(N);
        scalar rad = pspec(u) * rnd.GaussNormal();

        f[ index(  m, N-n) ] = 
                rad * std::complex<double>(Foam::cos(p),  Foam::sin(p));
        f[ index(M-m,   n) ] = 
                rad * std::complex<double>(Foam::cos(p), -Foam::sin(p));
      }
    }

    fftw_plan plan;
    plan = fftw_plan_dft_2d(M, N,
                             reinterpret_cast<fftw_complex*>(&f[0]),
                             reinterpret_cast<fftw_complex*>(&F[0]),
                             FFTW_BACKWARD,
                             FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);


    scalarField sF(MN);
    forAll(sF, ii)
    {
      sF[ii] = F[ii].real();
    }

    scalarField sF2 = sqr(sF);
    scalar avSF     = average(sF);
    scalar avSF2    = average(sF2);

    scalar factor = rgh / Foam::sqrt( mag(avSF2 - sqr(avSF)) );

    sF *= factor;

    disp = sF;
  }
};


int main(int argc, char *argv[])
{
  #include "setRootCase.H"
  #include "createTime.H"
  #include "createDynamicFvMesh.H"

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  // reading dictionary surfRoughGenDict
  IOdictionary surfRoughGenDict
  (
    IOobject
    (
      "surfRoughGenDict",
      runTime.system(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    )
  );
  
  word patchName;
  if( !surfRoughGenDict.readIfPresent<word>("patchName", patchName) )
  {
    SeriousErrorIn("main")
        << "There is no `patchName` parameter in surfRoughGenDict dictionary"
        << exit(FatalError);
  }
  
  int seed;
  if( !surfRoughGenDict.readIfPresent<int>("seed", seed) ){
    SeriousErrorIn("main")
        <<"There is no `seed` parameter in surfRoughGenDict dictionary"
        <<exit(FatalError);
  }

  int M;
  if( !surfRoughGenDict.readIfPresent<int>("sizeMaj", M) ){
    SeriousErrorIn("main")
        <<"There is no `sizeMaj` parameter in surfRoughGenDict dictionary"
        <<exit(FatalError);
  }
  int N;
  if( !surfRoughGenDict.readIfPresent<int>("sizeLat", N) ){
    SeriousErrorIn("main")
        <<"There is no `sizeLat` parameter in surfRoughGenDict dictionary"
        <<exit(FatalError);
  }
  
  int majDir;
  if( !surfRoughGenDict.readIfPresent<int>("majDir", majDir) ){
    SeriousErrorIn("main")
        <<"There is no `majDir` parameter in surfRoughGenDict dictionary"
        <<exit(FatalError);
  }
  int latDir;
  if( !surfRoughGenDict.readIfPresent<int>("latDir", latDir) ){
    SeriousErrorIn("main")
        <<"There is no `latDir` parameter in surfRoughGenDict dictionary"
        <<exit(FatalError);
  }
  
  scalar rgh;
  if( !surfRoughGenDict.readIfPresent<scalar>("roughness", rgh) ){
    SeriousErrorIn("main")
        <<"There is no `roughness` parameter in surfRoughGenDict dictionary"
        <<exit(FatalError);
  }
  // symmetric, synchronous, asymmetric
  word wayToApply;
  if( !surfRoughGenDict.readIfPresent<word>("apply", wayToApply) ){
    SeriousErrorIn("main")
        <<"There is no `synchronous` parameter in surfRoughGenDict dictionary"
        <<exit(FatalError);
  }
  
  double fractalParam;
  if( !surfRoughGenDict.readIfPresent<double>("fractalParam", fractalParam) ){
    SeriousErrorIn("main")
        <<"There is no `fractalParam` parameter in surfRoughGenDict dictionary"
        <<exit(FatalError);
  }
  
  double smoothing;
  if( !surfRoughGenDict.readIfPresent<double>("smoothing", smoothing) ){
    SeriousErrorIn("main")
        <<"There is no `smoothing` parameter in surfRoughGenDict dictionary"
        <<exit(FatalError);
  }
  
  Info<< "Patch:      " << patchName << endl;
  Info<< "Seed:       " << seed << endl;
  Info<< "M:          " << M << endl;
  Info<< "N:          " << N << endl;
  Info<< "roughness:  " << rgh << endl;
  Info<< "apply:      " << wayToApply << endl;
  Info<< "fractalParam:  " << fractalParam << endl;
  Info<< "smoothing:     " << smoothing << endl;

  Info<< "Setup RoughnessGenerator class" << endl;
  RoughnessGenerator rg(seed, M, N, rgh, wayToApply, fractalParam, smoothing);
  
  double cpuTime = runTime.elapsedCpuTime();
  //Info<< "Setup mesh relaxation class" << endl;
  //meshRelax mesh_rlx(mesh, args);
  
  // Get patch ID for boundaries we want to move ("walls" "inlet")
  label wallID  = mesh.boundaryMesh().findPatchID(patchName);
  if( wallID==-1 ){
    SeriousErrorIn("main")
        <<"patch "
        <<patchName
        <<" is missing"
        <<exit(FatalError);
  }
  
  coupledPatchInterpolation patchInterpolator( mesh.boundaryMesh()[wallID], mesh );
  
  const pointField& boundaryPoints = mesh.boundaryMesh()[wallID].localPoints();
  vectorField pointDispWall(boundaryPoints.size(), vector::zero);
  
  vectorField pointNface = mesh.boundaryMesh()[wallID].faceNormals();
  vectorField motionN = patchInterpolator.faceToPointInterpolate(pointNface);
  forAll(motionN, ii) motionN[ii]/=mag(motionN[ii]);
  
  scalarField faceDisp(pointNface.size(), 0.0);
  rg.getSurfaceDisplacement(mesh, faceDisp, wallID, majDir, latDir);
  scalarField pointDisp = patchInterpolator.faceToPointInterpolate(faceDisp);
  
  Info<<nl<< "Maximum and minimum face displacement in Y direction:"<<nl;
  Info<<nl<< "Max: " << max(faceDisp) << "  min: "<<min(faceDisp)<<nl<<nl;
  Info<<nl<< "Maximum and minimum point displacement in Y direction:"<<nl;
  Info<<     "Max: " << max(pointDisp) << "  min: "<<min(pointDisp)<<nl<<nl;
  
  forAll( pointDispWall, i )
    pointDispWall[i] = pointDisp[i] * motionN[i];
  
  pointDispWall /= runTime.deltaTValue();
  
  //bool auxSw = mesh_rlx.get_fixInletWallEdgeDispl();
  //mesh_rlx.set_fixInletWallEdgeDispl(false);
  //mesh_rlx.meshUpdate(pointDispWall, runTime);
  //mesh_rlx.set_fixInletWallEdgeDispl(auxSw);
  pointVectorField& pointVelocity = const_cast<pointVectorField&>
  (
    mesh.objectRegistry::lookupObject<pointVectorField>( "pointMotionU" )
  );
  pointVelocity.boundaryFieldRef()[wallID] == pointDispWall;
  mesh.update();
  
  cpuTime = runTime.elapsedCpuTime() - cpuTime;
  
  Info<<nl<<"Time statistics:"<<nl;
  
  int wlNP = mesh.boundaryMesh()[wallID].nPoints();
  Info<<"Total number of points:                   "<<mesh.nPoints()<<nl;
  Info<<"Number of points on the walls:            "<<wlNP<<nl;
  Info<<"Running time:                             "<<cpuTime<<nl<<endl;

  Info<<"Overwriting points in current time directory."<<nl;
  runTime.writeNow();
  
  Info<<"End"<<nl;
  return 0;
}

// **************************** End of the solver ******************************** //
