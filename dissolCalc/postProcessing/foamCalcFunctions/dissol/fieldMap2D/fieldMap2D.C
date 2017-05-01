/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.


Description
    contact info: vitaliy.starchenko@gmail.com
    
SourceFiles
    fieldMap2D.C

\*---------------------------------------------------------------------------*/

#define NUMBER_OF_COLUMNS 10

#include "fieldMap2D.H"

#include "addToRunTimeSelectionTable.H"
#include "memInfo.H"

#include "OFstreamMod.H"

#include "interpolation.H"

#include "triSurfaceSearch.H"
#include "meshSearch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
  namespace calcTypes
  {
    defineTypeNameAndDebug(fieldMap2D, 0);
    addToRunTimeSelectionTable(calcType, fieldMap2D, dictionary);
  }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::calcTypes::fieldMap2D::fieldMap2D()
:
  calcType()
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::calcTypes::fieldMap2D::~fieldMap2D()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::calcTypes::fieldMap2D::init()
{
  argList::validArgs.append("fieldMap2D");
  argList::validArgs.append("processingType");
  argList::validArgs.append("N");
  argList::validArgs.append("M");
}

void Foam::calcTypes::fieldMap2D::preCalc
(
  const argList& args,
  const Time& runTime,
  const fvMesh& mesh
)
{
  
  const word dictName("fieldMap2Ddict");
  IOdictionary dict
  (
    IOobject
    (
      dictName,
      runTime.system(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    )
  );

  if( !dict.readIfPresent<word>("patchName", patchName) )
  {
    SeriousErrorIn("preCalc")
          << "There is no patchName parameter in fieldMap2Ddict dictionary"
          << exit(FatalError);
  }
  
  if( !dict.readIfPresent<word>("geometry", geometry) )
  {
    SeriousErrorIn("preCalc")
          << "There is no geometry parameter in fieldMap2Ddict dictionary"
          << exit(FatalError);
  }
  
  point minPoint;
  if( !dict.readIfPresent<point>("minPoint", minPoint) )
  {
    SeriousErrorIn("preCalc")
          << "There is no minPoint parameter in fieldMap2Ddict dictionary"
          << exit(FatalError);
  }
  point maxPoint;
  if( !dict.readIfPresent<point>("maxPoint", maxPoint) )
  {
    SeriousErrorIn("preCalc")
          << "There is no maxPoint parameter in fieldMap2Ddict dictionary"
          << exit(FatalError);
  }

  int integrationDir;
  if( !dict.readIfPresent<int>("integrationDirection", integrationDir) )
  {
    SeriousErrorIn("preCalc")
          << "There is no integrationDirection parameter "
             "in fieldMap2Ddict dictionary"
          << exit(FatalError);
  }
  intDir = integrationDir;
  
  int flowDir;
  if( !dict.readIfPresent<int>("flowDirection", flowDir) )
  {
    SeriousErrorIn("preCalc")
          << "There is no flowDirection parameter "
             "in fieldMap2Ddict dictionary"
          << exit(FatalError);
  }
  majDir = flowDir;

  int expectedNumberOfIntersections;
  if
  (
    !dict.readIfPresent<int>
          (
            "expectedNumberOfIntersections", 
            expectedNumberOfIntersections
          ) 
  )
  {
    SeriousErrorIn("preCalc")
          << "There is no expectedNumberOfIntersections parameter "
             "in fieldMap2Ddict dictionary"
          << exit(FatalError);
  }
  expNI = expectedNumberOfIntersections;
  
  int numberOfIntegrationPoints;
  if
  (
    !dict.readIfPresent<int>
          (
            "numberOfIntegrationPoints",
            numberOfIntegrationPoints
          ) 
  )
  {
    SeriousErrorIn("preCalc")
          << "There is no numberOfIntegrationPoints parameter "
             "in fieldMap2Ddict dictionary"
          << exit(FatalError);
  }
  
  
  #ifdef FOAM_DEV
    const word& processingType = args.additionalArgs()[1];
    const word& Nword = args.additionalArgs()[2];
    const word& Mword = args.additionalArgs()[3];
    const word& Kword = args.additionalArgs()[4];

    N_ = std::atoi( Nword.c_str() );
    M_ = std::atoi( Mword.c_str() );
    Info << "FOAM_DEV true: "<< FOAM_DEV <<nl;
  #else
    const word processingType = args[2];
    N_ = std::atoi( args[3].c_str() );
    M_ = std::atoi( args[4].c_str() );
    Info << "FOAM_DEV false: " <<nl;
  #endif
  processingType_ = const_cast<word&>(processingType);
  latDir = 3 - (intDir + majDir);
  
  Info<<"* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"
          <<nl
          <<"Post-processing parameters:"<<nl
          <<"patchName:                 "<< patchName<<nl
          <<"geometry:                  "<< geometry<<nl
          <<"minPoint:                  "<< minPoint<<nl
          <<"maxPoint:                  "<< maxPoint<<nl
          <<"major direction:           "<< majDir<<nl
          <<"lateral direction:         "<< latDir<<nl
          <<"integration direction:     "<< intDir<<nl
          <<"number of intersections:   "<<expNI<<nl
          <<"number of integration points: "<< numberOfIntegrationPoints
          <<endl;
  
  if
  (
    processingType_!="ccAll" && geometry == "concentricCylinders"
  )
  {
    SeriousErrorIn("preCalc")
          << "Concentric cylinders geometry should have proc type ccAll"
          << exit(FatalError);
  }

  N1_ = N_+1;
  M1_ = M_+1;
  K_  = numberOfIntegrationPoints;
  K1_ = numberOfIntegrationPoints+1;
  
  N1M1 = N1_*M1_;

  // need temporary arguments in order to add a 0 time directory
  /*
  argList argsTmp = args;
  argsTmp.setOption("time", "0");
  Foam::Time timeTmp(Foam::Time::controlDictName, argsTmp);
  Foam::instantList timeDirs = Foam::timeSelector::select0(timeTmp, argsTmp);
  timeTmp.setTime(timeDirs[0], 0);
  
  Foam::fvMesh meshTmp
  (
      Foam::IOobject
      (
          Foam::fvMesh::defaultRegion,
          timeTmp.timeName(),
          timeTmp,
          Foam::IOobject::MUST_READ
      )
  );
  if( timeTmp.timeName()!="0" )
  {
    SeriousErrorIn("fieldOperations::getInletAreaT0")
            <<"There is no 0 time directory. Check your decomposition as well!"<<nl
            <<"TimeTmp: "<< timeTmp.timeName()<<nl
            <<exit(FatalError);
  }
  */
  
  //Info << "Calculating the max and min coordinates"<<nl;
  //const pointField &allPoints = meshTmp.points();
  //minPosX = min( allPoints.component(vector::X) );
  //minPosY = min( allPoints.component(vector::Y) );
  //minPosZ = min( allPoints.component(vector::Z) );
  
  if(geometry=="flat")
  {
    minPosMaj = minPoint.component(majDir);
    minPosLat = minPoint.component(latDir);
    minPosInt = minPoint.component(intDir);

    maxPosMaj = maxPoint.component(majDir);
    maxPosLat = maxPoint.component(latDir);
    maxPosInt = maxPoint.component(intDir);

    // z becomes x; x becomes y
    dx = (maxPosMaj - minPosMaj) / static_cast<scalar>(N_);
    dy = (maxPosLat - minPosLat) / static_cast<scalar>(M_);
  }
  else if(geometry=="concentricCylinders")
  {
    minPosMaj = minPoint.component(majDir);
    minPosLat = 0.0;
    minPosInt = 0.0;

    maxPosMaj = maxPoint.component(majDir);
    maxPosLat = constant::mathematical::twoPi;
    maxPosInt = maxPoint.component(intDir);

    // in case of concentric cylinders we go around the circle
    dx = (maxPosMaj - minPosMaj) / static_cast<scalar>(N_);
    dy = (maxPosLat - minPosLat) / static_cast<scalar>(M_);
  }
  else
  {
    FatalError<<"Unknown geometry"<<nl<<exit(FatalError);
  }
  
  
  // calculate the size of the current pattern being not bigger then 10^7
  thisTimeSize = min(N1M1, maxNumProcPoints);
  curNum = 0;
  curEnd = thisTimeSize;
  
  totNumLoop = (N1M1-1) / maxNumProcPoints + 1;
}


// TODO make it more general
void Foam::calcTypes::fieldMap2D::calc
(
  const argList& args,
  const Time& runTime,
  const fvMesh& mesh
)
{
  // coordinates of points on the surface of the fracture walls
  //memInfo mf;
  
  // print what is calculated
  if(processingType_ == "all")
    Info << "Processing all fields..." << endl;
  else if(processingType_ == "surf")
    Info << "Calculating csurf and h" << endl;
  else if(processingType_ == "int")
    Info << "Calculating C, qx and qy" << endl;
  else if(processingType_ == "h")
    Info << "Processing the aperture of the fracture..." << endl;
  else if(processingType_ == "U")
    Info << "Processing the flux qx and qy..." << endl;
  else if(processingType_ == "p")
    FatalError<<"p processing is not implemented yet"<<nl<<exit(FatalError);
  else if(processingType_ == "C")
    Info << "Processing the concentration field..." << endl;
  else if(processingType_ == "csurf")
    Info << "Calculating concentration on the surface" << endl;
  else if(processingType_ == "temp")
    Info << "Running temporary function..." << endl;
  else if(processingType_ == "ccAll")
    Info << "Processing all fields for concentric cylinder geometry" << endl;
  else
    FatalError<<"Unable to process "<<processingType_<<nl<<exit(FatalError);
  
  for(int cI=0; cI<totNumLoop; cI++)
  {
    curNum = cI;
    curBlock = thisTimeSize * curNum;
    
    sizeAA = thisTimeSize;
    if(cI==totNumLoop-1){
      sizeAA = N1M1 - (totNumLoop-1)*thisTimeSize;
    }

    Info << "Find the points on the surface"<<nl;
    pointsXYonsurface.clear();
    pointsXYonsurface.setSize( expNI * sizeAA );
    
    if(geometry=="flat")
    {
      build_surface_points( mesh );
    }
    else if(geometry=="concentricCylinders")
    {
      build_surface_pointsCC( mesh );
    }
    Info << "build_surface_points done"<<nl;
    
    fileName current_dissolCalc_dir;
    current_dissolCalc_dir = "postProcessing/dissolCalc" / runTime.timeName();
    if ( !isDir(current_dissolCalc_dir) ) mkDir(current_dissolCalc_dir);

    if(processingType_ == "all"){
      write_all(mesh, runTime);
    }
    else if(processingType_ == "surf"){
      write_surf(mesh, runTime);
    }
    else if(processingType_ == "int"){
      write_int(mesh, runTime);
    }
    else if(processingType_ == "h"){
      write_h(mesh, runTime);
    }
    else if(processingType_ == "U"){
      write_q(mesh, runTime);
    }
    else if(processingType_ == "C"){
      write_Ccup(mesh, runTime);
    }
    else if(processingType_ == "csurf"){
      write_csurf(mesh, runTime);
    }
    else if(processingType_ == "temp"){
      write_temp(mesh, runTime);
    }
    else if(processingType_ == "ccAll"){
      write_ccAll(mesh, runTime);
    }
    else{
      FatalError<<"Unable to process "<<processingType_<<nl<<nl<<exit(FatalError);
    }
  }
}

scalar Foam::calcTypes::fieldMap2D::primitive_simpson_integration
(
  scalarField& x,
  scalarField& y  
)
{
  scalar result = 0.0;
  
  // TODO check len y == len x
  int step = 2;
  
  scalarField dlt_x( x.size()-1 );
  
  forAll(dlt_x, i){
    dlt_x = x[i+1]-x[i];
  }
  
  int len_lab_list = static_cast<int>( std::floor( dlt_x.size() / 2.0 ) );
  
  labelList slice0( len_lab_list );
  labelList slice1( len_lab_list );
  labelList slice2( len_lab_list );
  
  int count = 0;
  forAll(slice0, i){
    slice0[i] = count;
    slice1[i] = count+1;
    slice2[i] = count+2;
    count += step;
  }
  
  scalarField dx0(len_lab_list);
  scalarField dx1(len_lab_list);
  forAll(dx0, j){
    dx0[j] = dlt_x[slice0[j]];
    dx1[j] = dlt_x[slice1[j]];
  }

  scalarField hsum(len_lab_list);     //= h0 + h1
  scalarField hprod(len_lab_list);    //= h0 * h1
  scalarField h0divh1(len_lab_list);  //= h0 / h1

  
  forAll(dx0, i){
    hsum[i] = dx0[i] + dx1[i];
    hprod[i] = dx0[i] * dx1[i];
    h0divh1[i] = dx0[i] / dx1[i];
  }
  
  forAll(hsum, i){
    result += hsum[i]/6.0 *
                          (
                            y[slice0[i]] * (2-1.0/h0divh1[i]) +
                            y[slice1[i]] * hsum[i]*hsum[i]/hprod[i] +
                            y[slice2[i]] * (2-h0divh1[i])
                          );
  }
  
  return result;
}



void Foam::calcTypes::fieldMap2D::build_surface_points
(
  const fvMesh& mesh
)
{
  // triangulation
  // @TODO in case of parallel calculations see surfaceMeshTriangulate.C
  // @TODO add the error handling for "walls"
  label patchID = mesh.boundaryMesh().findPatchID(patchName);
  
  if(patchID==-1)
  {
    SeriousErrorIn("build_surface_points")
            <<"There is no "
            << patchName
            << exit(FatalError);
  }
  
  labelHashSet includePatches(1);
  includePatches.insert(patchID);
  
  triSurface wallTriSurface
  (
    triSurfaceTools::triangulate( mesh.boundaryMesh(), includePatches )
  );
  
  // intersection
  const triSurfaceSearch querySurf(wallTriSurface);
  const indexedOctree<treeDataTriSurface>& tree = querySurf.tree();
  
  scalar curMaj = 0.0;
  scalar curLat = 0.0;
  
  const pointField &localPoints = mesh.boundaryMesh()[patchID].localPoints();
  scalar maxMaj = max( localPoints.component(majDir) );
  scalar minMaj = min( localPoints.component(majDir) );
  scalar maxLat = max( localPoints.component(latDir) );
  scalar minLat = min( localPoints.component(latDir) );
  
  // N+1 and M+1 because we need points at x=minPosX and x=maxPosX
  for(int ijk=0; ijk<sizeAA; ijk++)
  {
    int totIJK = ijk + curNum * thisTimeSize;
    label i = totIJK / M1_;
    label j = totIJK % M1_;
    curMaj = minPosMaj + i*dx;
    curLat = minPosLat + j*dy;
    label ind = ijk;
      
    scalar eps=1e-2;
    point searchStart;
    searchStart.component(majDir) = curMaj;
    searchStart.component(latDir) = curLat;
    searchStart.component(intDir) = minPosInt;
    
    while(searchStart.component(majDir) >= maxMaj ) 
      searchStart.component(majDir) -= eps;
    while(searchStart.component(majDir) <= minMaj ) 
      searchStart.component(majDir) += eps;
    while(searchStart.component(latDir) >= maxLat ) 
      searchStart.component(latDir) -= eps;
    while(searchStart.component(latDir) <= minLat ) 
      searchStart.component(latDir) += eps;
    
    point searchEnd = searchStart;
    searchEnd.component(intDir) = maxPosInt +
            mag(maxPosInt-minPosInt)/(maxPosInt-minPosInt);

    point hitPoint(0.0, 0.0, 0.0);
    
    pointIndexHit pHit = tree.findLine(searchStart, searchEnd);
    
    if ( pHit.hit() )
    {
      hitPoint = pHit.hitPoint();
    }
    
    label hitWallLabel = expNI*ind;
    
    pointsXYonsurface[hitWallLabel] = hitPoint;
    
    if(expNI>1)
    {
      // search for second intersection
      searchStart.component(intDir) = hitPoint.component(intDir) + eps;
      pHit = tree.findLine(searchStart, searchEnd);
      label secondHitWallLabel = expNI * ind + 1;
      if ( pHit.hit() )
      {
        pointsXYonsurface[secondHitWallLabel] = pHit.hitPoint();
      }
      else
      {
        pointsXYonsurface[hitWallLabel] = point::zero;
        pointsXYonsurface[secondHitWallLabel] = point::zero;
      }

    }

  }
}

void Foam::calcTypes::fieldMap2D::build_surface_pointsCC
(
  const fvMesh& mesh
)
{
  // triangulation
  // @TODO in case of parallel calculations see surfaceMeshTriangulate.C
  // @TODO add the error handling for "walls"
  label patchID = mesh.boundaryMesh().findPatchID(patchName);
  
  if(patchID==-1)
  {
    SeriousErrorIn("build_surface_points")
            <<"There is no "
            << patchName
            << exit(FatalError);
  }
  
  labelHashSet includePatches(1);
  includePatches.insert(patchID);
  
  triSurface wallTriSurface
  (
    triSurfaceTools::triangulate( mesh.boundaryMesh(), includePatches )
  );
  
  // intersection
  const triSurfaceSearch querySurf(wallTriSurface);
  const indexedOctree<treeDataTriSurface>& tree = querySurf.tree();
  
  scalar curMaj = 0.0;
  scalar curLat = 0.0;
  
  // N+1 and M+1 because we need points at x=minPosX and x=maxPosX
  for(int ijk=0; ijk<sizeAA; ijk++)
  {
    int totIJK = ijk + curNum * thisTimeSize;
    label i = totIJK / M1_;
    label j = totIJK % M1_;
    curMaj = minPosMaj + i*dx;
    curLat = minPosLat + j*dy; // here it is an angle [0;2Pi]
    label ind = ijk;
      
    point searchStart;
    searchStart.component(majDir) = curMaj;
    searchStart.component(latDir) = 0.0;
    searchStart.component(intDir) = 0.0;
    
    
    point searchEnd = searchStart;
    searchEnd.component(intDir) = maxPosInt +
            mag(maxPosInt-minPosInt)/(maxPosInt-minPosInt);
    
    searchEnd.component(latDir) = maxPosInt * Foam::cos(curLat);
    searchEnd.component(intDir) = maxPosInt * Foam::sin(curLat);

    point hitPoint(0.0, 0.0, 0.0);
    
    pointIndexHit pHit = tree.findLine(searchStart, searchEnd);
    
    if ( pHit.hit() )
    {
      hitPoint = pHit.hitPoint();
    }
    
    label hitWallLabel = expNI*ind;
    
    pointsXYonsurface[hitWallLabel] = hitPoint;
  }
}




// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// an empty function for something extremely not general
void Foam::calcTypes::fieldMap2D::write_temp
(
  const fvMesh& mesh,
  const Time& runTime
)
{
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void Foam::calcTypes::fieldMap2D::write_csurf
(
  const fvMesh& mesh,
  const Time& runTime
)
{
  fileName current_file_path;
  current_file_path = "postProcessing/dissolCalc"/runTime.timeName()/"csurf";
  
  ios_base::openmode mode = (curNum==0) ? 
    ios_base::out|ios_base::trunc 
          : 
    ios_base::out|ios_base::app;
  
  OFstreamMod aFile( current_file_path, mode);
  
  if( curNum==0 )
  {
    aFile << N_ << "   " << M_ << "   " << runTime.value() 
            << "   " << dx << "   " << dy << endl;
  }

  typedef GeometricField<scalar, fvPatchField, volMesh> fieldC;
  IOobject headerC
  (
      "C",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ
  );
  fieldC field_c(headerC, mesh);
  autoPtr<interpolation<scalar> >interpolatorC
  (
    interpolation<scalar>::New("cellPoint",field_c)
  );
  meshSearch searchEngine(mesh);
  
  int count = 0;
  int i = 0;
  while (i<pointsXYonsurface.size() )
  {
    point first = pointsXYonsurface[i];
    
    label cellI = searchEngine.findCell( first );
    if (cellI==-1) cellI = searchEngine.findNearestCell( first );
    scalar ac = interpolatorC->interpolate(first, cellI, -1);
    
    aFile << ac << "  ";
    
    i+=expNI;
    count++;
    if(count>=NUMBER_OF_COLUMNS){ 
      aFile <<"\n";
      count=0;
    }
  }
}


void Foam::calcTypes::fieldMap2D::write_h
(
  const fvMesh& mesh,
  const Time& runTime
)
{
  fileName current_file_path;
  current_file_path = "postProcessing/dissolCalc" / runTime.timeName() / "h";
  
  ios_base::openmode mode = (curNum==0) ? 
    ios_base::out|ios_base::trunc 
          : 
    ios_base::out|ios_base::app;
  
  OFstreamMod apertureMapFile( current_file_path, mode);
  
  if( curNum==0 )
  {
    apertureMapFile << N_ << "   " << M_ << "   " << runTime.value() 
            << "   " << dx << "   " << dy << endl;
  }

  int count = 0;
  int i = 0;
  while (i<pointsXYonsurface.size() )
  {
    point first = pointsXYonsurface[i];
    point second = first;
    
    if(expNI == 1)
    {
      second.component(intDir) = minPosInt;
    }
    else
    {
      second = pointsXYonsurface[i+1];
    }
    
    apertureMapFile << mag( first - second ) << "  ";
    
    i+=expNI;
    count++;
    if(count>=NUMBER_OF_COLUMNS){ 
      apertureMapFile <<"\n";
      count=0;
    }
  }
}

void Foam::calcTypes::fieldMap2D::write_q
(
  const fvMesh& mesh,
  const Time& runTime
)
{
  typedef GeometricField<vector, fvPatchField, volMesh> fieldU;

  IOobject headerU
  (
      "U",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ
  );
  
  //if (headerU.headerOk())
  {
    fieldU field_u(headerU, mesh);

    fileName current_file_path_qx, current_file_path_qy;
    current_file_path_qx =
            "postProcessing/dissolCalc" / runTime.timeName() / "qx";
    current_file_path_qy =
            "postProcessing/dissolCalc" / runTime.timeName() / "qy";
  
    autoPtr<interpolation<vector> >interpolatorU
    (
      interpolation<vector>::New("cellPoint",field_u)
    );

    ios_base::openmode mode =
            (curNum==0) ? 
              ios_base::out|ios_base::trunc : 
              ios_base::out|ios_base::app;  
    OFstreamMod mapXYqx( current_file_path_qx, mode );
    OFstreamMod mapXYqy( current_file_path_qy, mode );
    
    if(curNum==0)
    {
      mapXYqx << N_ << "   " << M_ << "   " << runTime.value()
              << "   " << dx << "   " << dy << endl;
      mapXYqy << N_ << "   " << M_ << "   " << runTime.value()
              << "   " << dx << "   " << dy << endl;
    }
    
    meshSearch searchEngine(mesh);
        
    int count = 0;
    int iSurf = 0;
    while (iSurf<pointsXYonsurface.size() )
    {
      point first = pointsXYonsurface[iSurf];
      point second = first;

      if(expNI == 1)
      {
        second.component(intDir) = minPosInt;
      }
      else
      {
        second = pointsXYonsurface[iSurf+1];
      }
      
      vector dirc =  first - second;
      
      scalar integratedUx = 0.0;
      scalar integratedUy = 0.0;
      if( mag(dirc) > SMALL )
      {
        scalarField variable(K1_);
        scalarField Ux(K1_);
        scalarField Uy(K1_);
        
        vector dr = 1.0 / static_cast<scalar>(K_) * dirc;

        for(int i=0; i<K1_; i++)
        {
          point samp_point = second + i*dr;
          label cellI = searchEngine.findCell( samp_point );
          if (cellI==-1)
          {
            cellI = searchEngine.findNearestCell( samp_point );
          }
          label faceI = -1;

          vector interp_fieldU =
                  interpolatorU->interpolate(samp_point, cellI, faceI);
          // velocity = 0 on the surface
          if(i==0 || i==K_)
          {
            interp_fieldU = vector::zero;
          }
          
          variable[i] = samp_point.component(intDir);
          Ux[i] = interp_fieldU.component(majDir);
          Uy[i] = interp_fieldU.component(latDir);
        }
        
        integratedUx = primitive_simpson_integration(variable, Ux);
        integratedUy = primitive_simpson_integration(variable, Uy);
      }
      
      int ind1 = iSurf + curBlock;
      int prcnt  = 100 * (ind1    ) / N1M1;
      int prcnt1 = 100 * (ind1 - 1) / N1M1;
      if( prcnt%1==0 && prcnt!=prcnt1 )
      {
        Info<<"\r"<< prcnt << "%  "<<flush;
      }
      
      mapXYqx << integratedUx << "  ";
      mapXYqy << integratedUy << "  ";

      iSurf += expNI;
      count++;
      if(count >= NUMBER_OF_COLUMNS)
      { 
        mapXYqx <<"\n";
        mapXYqy <<"\n";
        count=0;
      }
    }  
  }
  /*
  else
  {
    FatalError<<"There is no U field"<<nl<<nl<<exit(FatalError);
  }
  */
}

void Foam::calcTypes::fieldMap2D::write_Ccup
(
  const fvMesh& mesh,
  const Time& runTime
)
{
  typedef GeometricField<vector, fvPatchField, volMesh> fieldU;
  typedef GeometricField<scalar, fvPatchField, volMesh> fieldC;

  IOobject headerU
  (
      "U",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ
  );
  IOobject headerC
  (
      "C",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ
  );
  
  //if (headerU.headerOk() && headerC.headerOk())
  {
    fieldU field_u(headerU, mesh);
    fieldC field_c(headerC, mesh);

    fileName current_file_path_cup = 
            "postProcessing/dissolCalc" / runTime.timeName() / "c";
  
    autoPtr<interpolation<vector> >interpolatorU
    (
      interpolation<vector>::New("cellPoint",field_u)
    );
    autoPtr<interpolation<scalar> >interpolatorC
    (
      interpolation<scalar>::New("cellPoint",field_c)
    );

    ios_base::openmode mode = 
            (curNum==0) ? 
              ios_base::out|ios_base::trunc : 
              ios_base::out|ios_base::app;  
    OFstreamMod mapXYccup( current_file_path_cup, mode );
    
    if(curNum==0)
    {
      mapXYccup << N_ << "   " << M_ << "   " << runTime.value() 
              << "   " << dx << "   " << dy << endl;
    }
    
    meshSearch searchEngine(mesh);
        
    int count = 0;
    int iSurf = 0;
    while (iSurf<pointsXYonsurface.size() )
    {
      point first = pointsXYonsurface[iSurf];
      point second = first;

      if(expNI == 1)
      {
        second.component(intDir) = minPosInt;
      }
      else
      {
        second = pointsXYonsurface[iSurf+1];
      }
      
      vector dirc =  first - second;
      
      scalar integratedUx = 0.0;
      scalar integratedUy = 0.0;
      scalar integratedUCx = 0.0;
      scalar integratedUCy = 0.0;
      scalar Ccup = 0.0;
      if( mag(dirc) > SMALL )
      {
        scalarField variable(K1_);
        scalarField Ux(K1_);
        scalarField Uy(K1_);
        scalarField UCx(K1_);
        scalarField UCy(K1_);
        
        vector dr = 1.0 / static_cast<scalar>(K_) * dirc;

        for(int i=0; i<K1_; i++)
        {
          point samp_point = second + i*dr;
          label cellI = searchEngine.findCell( samp_point );
          if (cellI==-1)
          {
            cellI = searchEngine.findNearestCell( samp_point );
          }
          label faceI = -1;

          vector interp_fieldU =
                  interpolatorU->interpolate(samp_point, cellI, faceI);
          // velocity = 0 on the surface
          if(i==0 || i==K_)
          {
            interp_fieldU = vector::zero;
          }
          
          variable[i] = samp_point.component(intDir);
          Ux[i] = interp_fieldU.component(majDir);
          Uy[i] = interp_fieldU.component(latDir);
          
          scalar interp_fieldC = 
                  interpolatorC->interpolate(samp_point, cellI, faceI);
          UCx[i] = Ux[i] * interp_fieldC;
          UCy[i] = Uy[i] * interp_fieldC;
        }
        
        integratedUx = primitive_simpson_integration(variable, Ux);
        integratedUy = primitive_simpson_integration(variable, Uy);
        integratedUCx = primitive_simpson_integration(variable, UCx);
        integratedUCy = primitive_simpson_integration(variable, UCy);
        scalar qSqr = integratedUx*integratedUx+integratedUy*integratedUy;
        if( qSqr > SMALL )
        {
          Ccup = std::sqrt
                 (
                  (integratedUCx*integratedUCx+integratedUCy*integratedUCy)
                  /
                  qSqr
                 );
        }
      }
      
      int ind1 = iSurf + curBlock;
      int prcnt  = 100 * (ind1    ) / N1M1;
      int prcnt1 = 100 * (ind1 - 1) / N1M1;
      if( prcnt%1==0 && prcnt!=prcnt1 )
      {
        Info<<"\r"<< prcnt << "%  "<<flush;
      }
      
      mapXYccup << Ccup << "  ";

      iSurf += expNI;
      count++;
      if(count >= NUMBER_OF_COLUMNS)
      { 
        mapXYccup <<"\n";
        count=0;
      }
    }  
  }
  
  /*
  else if( !headerU.headerOk() && headerC.headerOk() )
  {
    FatalError<<"There is no U field"<<nl<<nl<<exit(FatalError);
  }
  else if( headerU.headerOk() && !headerC.headerOk() )
  {
    FatalError<<"There is no C field"<<nl<<nl<<exit(FatalError);
  }
  else
  {
    FatalError<<"There is no U and C field"<<nl<<nl<<exit(FatalError);
  }
  */
}



void Foam::calcTypes::fieldMap2D::write_surf
(
  const fvMesh& mesh,
  const Time& runTime
)
{
  typedef GeometricField<scalar, fvPatchField, volMesh> fieldC;

  IOobject headerC
  (
      "C",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ
  );
  
  //if (headerC.headerOk())
  {
    fieldC field_c(headerC, mesh);

    fileName current_file_path_h =
            "postProcessing/dissolCalc" / runTime.timeName() / "h";
    fileName current_file_path_csurf =
            "postProcessing/dissolCalc" / runTime.timeName() / "csurf";
  
    autoPtr<interpolation<scalar> >interpolatorC
    (
      interpolation<scalar>::New("cellPoint",field_c)
    );

    ios_base::openmode mode = 
            (curNum==0) ?
              ios_base::out|ios_base::trunc : 
              ios_base::out|ios_base::app;  
    OFstreamMod mapXYcsurf( current_file_path_csurf, mode);
    OFstreamMod apertureMapFile( current_file_path_h, mode);
    
    if(curNum==0)
    {
      mapXYcsurf << N_ << "   " << M_ << "   " << runTime.value() 
              << "   " << dx << "   " << dy << endl;
      apertureMapFile << N_ << "   " << M_ << "   " << runTime.value() 
              << "   " << dx << "   " << dy << endl;
    }
    
    meshSearch searchEngine(mesh);
        
    int count = 0;
    int iSurf = 0;
    while (iSurf<pointsXYonsurface.size() )
    {
      point first = pointsXYonsurface[iSurf];
      point second = first;

      if(expNI == 1)
      {
        second.component(intDir) = minPosInt;
      }
      else
      {
        second = pointsXYonsurface[iSurf+1];
      }
      
      vector dirc =  first - second;
      scalar dist = mag( dirc );
      
      scalar csurf = 0.0;
      if( mag(dirc) > SMALL ){
        label cellI = searchEngine.findCell( first );
        if (cellI==-1)
        {
          cellI = searchEngine.findNearestCell( first );
        }
                  
        csurf = interpolatorC->interpolate(first, cellI, -1);
      }
      
      int ind1 = iSurf + curBlock;
      int prcnt  = 100 * (ind1    ) / N1M1;
      int prcnt1 = 100 * (ind1 - 1) / N1M1;
      if( prcnt%1==0 && prcnt!=prcnt1 )
      {
        Info<<"\r"<< prcnt << "%  "<<flush;
      }
      
      mapXYcsurf << csurf << "  ";
      apertureMapFile << dist << "  ";

      iSurf += expNI;
      count++;
      if(count >= NUMBER_OF_COLUMNS)
      { 
        mapXYcsurf <<"\n";
        apertureMapFile <<"\n";
        count=0;
      }
    }  
  }
  /*
  else
  {
    FatalError<<"There is no C field"<<nl<<nl<<exit(FatalError);
  }
  */
}


void Foam::calcTypes::fieldMap2D::write_int
(
  const fvMesh& mesh,
  const Time& runTime
)
{
  typedef GeometricField<vector, fvPatchField, volMesh> fieldU;
  typedef GeometricField<scalar, fvPatchField, volMesh> fieldC;

  IOobject headerU
  (
      "U",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ
  );
  IOobject headerC
  (
      "C",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ
  );
  
  //if (headerU.headerOk() && headerC.headerOk())
  {
    fieldU field_u(headerU, mesh);
    fieldC field_c(headerC, mesh);

    fileName current_file_path_cup =
            "postProcessing/dissolCalc" / runTime.timeName() / "c";
    fileName current_file_path_qx, current_file_path_qy;
    current_file_path_qx =
            "postProcessing/dissolCalc" / runTime.timeName() / "qx";
    current_file_path_qy =
            "postProcessing/dissolCalc" / runTime.timeName() / "qy";
  
    autoPtr<interpolation<vector> >interpolatorU
    (
      interpolation<vector>::New("cellPoint",field_u)
    );
    autoPtr<interpolation<scalar> >interpolatorC
    (
      interpolation<scalar>::New("cellPoint",field_c)
    );

    ios_base::openmode mode =
            (curNum==0) ? 
              ios_base::out|ios_base::trunc : 
              ios_base::out|ios_base::app;  
    OFstreamMod mapXYccup( current_file_path_cup, mode );
    OFstreamMod mapXYqx( current_file_path_qx, mode );
    OFstreamMod mapXYqy( current_file_path_qy, mode );
    
    if(curNum==0)
    {
      mapXYccup << N_ << "   " << M_ << "   " << runTime.value() 
              << "   " << dx << "   " << dy << endl;
      mapXYqx << N_ << "   " << M_ << "   " << runTime.value()
              << "   " << dx << "   " << dy << endl;
      mapXYqy << N_ << "   " << M_ << "   " << runTime.value()
              << "   " << dx << "   " << dy << endl;
    }
    
    meshSearch searchEngine(mesh);
        
    int count = 0;
    int iSurf = 0;
    while (iSurf<pointsXYonsurface.size() )
    {
      point first = pointsXYonsurface[iSurf];
      point second = first;

      if(expNI == 1)
      {
        second.component(intDir) = minPosInt;
      }
      else
      {
        second = pointsXYonsurface[iSurf+1];
      }
      
      vector dirc =  first - second;
      
      scalar integratedUx = 0.0;
      scalar integratedUy = 0.0;
      scalar integratedUCx = 0.0;
      scalar integratedUCy = 0.0;
      scalar Ccup = 0.0;
      if( mag(dirc) > SMALL )
      {
        scalarField variable(K1_);
        scalarField Ux(K1_);
        scalarField Uy(K1_);
        scalarField UCx(K1_);
        scalarField UCy(K1_);
        
        vector dr = 1.0 / static_cast<scalar>(K_) * dirc;

        for(int i=0; i<K1_; i++)
        {
          point samp_point = second + i*dr;
          label cellI = searchEngine.findCell( samp_point );
          if (cellI==-1)
          {
            cellI = searchEngine.findNearestCell( samp_point );
          }
          label faceI = -1;

          vector interp_fieldU =
                  interpolatorU->interpolate(samp_point, cellI, faceI);
          // velocity = 0 on the surface
          if(i==0 || i==K_)
          {
            interp_fieldU = vector::zero;
          }
          
          variable[i] = samp_point.component(intDir);
          Ux[i] = interp_fieldU.component(majDir);
          Uy[i] = interp_fieldU.component(latDir);
          
          scalar interp_fieldC = 
                  interpolatorC->interpolate(samp_point, cellI, faceI);
          UCx[i] = Ux[i] * interp_fieldC;
          UCy[i] = Uy[i] * interp_fieldC;
        }
        
        integratedUx = primitive_simpson_integration(variable, Ux);
        integratedUy = primitive_simpson_integration(variable, Uy);
        integratedUCx = primitive_simpson_integration(variable, UCx);
        integratedUCy = primitive_simpson_integration(variable, UCy);
        scalar qSqr = integratedUx*integratedUx+integratedUy*integratedUy;
        if( qSqr > SMALL )
        {
          Ccup = std::sqrt
                 (
                  (integratedUCx*integratedUCx+integratedUCy*integratedUCy)
                  /
                  qSqr
                 );
        }
      }
      
      int ind1 = iSurf + curBlock;
      int prcnt  = 100 * (ind1    ) / N1M1;
      int prcnt1 = 100 * (ind1 - 1) / N1M1;
      if( prcnt%1==0 && prcnt!=prcnt1 )
      {
        Info<<"\r"<< prcnt << "%  "<<flush;
      }
      
      mapXYccup << Ccup << "  ";
      mapXYqx << integratedUx << "  ";
      mapXYqy << integratedUy << "  ";

      iSurf += expNI;
      count++;
      if(count >= NUMBER_OF_COLUMNS)
      { 
        mapXYccup <<"\n";
        mapXYqx <<"\n";
        mapXYqy <<"\n";
        count=0;
      }
    }  
  }
  /*
  else if( !headerU.headerOk() && headerC.headerOk() )
  {
    FatalError<<"There is no U field"<<nl<<nl<<exit(FatalError);
  }
  else if( headerU.headerOk() && !headerC.headerOk() )
  {
    FatalError<<"There is no C field"<<nl<<nl<<exit(FatalError);
  }
  else
  {
    FatalError<<"There is no U and C field"<<nl<<nl<<exit(FatalError);
  }
  */
}

void Foam::calcTypes::fieldMap2D::write_all
(
  const fvMesh& mesh,
  const Time& runTime
)
{
  typedef GeometricField<vector, fvPatchField, volMesh> fieldU;
  typedef GeometricField<scalar, fvPatchField, volMesh> fieldC;

  IOobject headerU
  (
      "U",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ
  );
  IOobject headerC
  (
      "C",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ
  );
  
  //if (headerU.headerOk() && headerC.headerOk())
  {
    fieldU field_u(headerU, mesh);
    fieldC field_c(headerC, mesh);

    fileName current_file_path_cup =
            "postProcessing/dissolCalc" / runTime.timeName() / "c";
    fileName current_file_path_qx, current_file_path_qy;
    current_file_path_qx =
            "postProcessing/dissolCalc" / runTime.timeName() / "qx";
    current_file_path_qy =
            "postProcessing/dissolCalc" / runTime.timeName() / "qy";
    fileName current_file_path_h =
            "postProcessing/dissolCalc" / runTime.timeName() / "h";
    fileName current_file_path_csurf =
            "postProcessing/dissolCalc" / runTime.timeName() / "csurf";
  
    autoPtr<interpolation<vector> >interpolatorU
    (
      interpolation<vector>::New("cellPoint",field_u)
    );
    autoPtr<interpolation<scalar> >interpolatorC
    (
      interpolation<scalar>::New("cellPoint",field_c)
    );

    ios_base::openmode mode =
            (curNum==0) ? 
              ios_base::out|ios_base::trunc : 
              ios_base::out|ios_base::app;  
    OFstreamMod mapXYccup( current_file_path_cup, mode );
    OFstreamMod mapXYqx( current_file_path_qx, mode );
    OFstreamMod mapXYqy( current_file_path_qy, mode );
    OFstreamMod mapXYcsurf( current_file_path_csurf, mode);
    OFstreamMod apertureMapFile( current_file_path_h, mode);
    
    if(curNum==0)
    {
      mapXYccup << N_ << "   " << M_ << "   " << runTime.value() 
              << "   " << dx << "   " << dy << endl;
      mapXYcsurf << N_ << "   " << M_ << "   " << runTime.value() 
              << "   " << dx << "   " << dy << endl;
      mapXYqx << N_ << "   " << M_ << "   " << runTime.value()
              << "   " << dx << "   " << dy << endl;
      mapXYqy << N_ << "   " << M_ << "   " << runTime.value()
              << "   " << dx << "   " << dy << endl;
      apertureMapFile << N_ << "   " << M_ << "   " << runTime.value() 
              << "   " << dx << "   " << dy << endl;
    }
    
    meshSearch searchEngine(mesh);
        
    int count = 0;
    int iSurf = 0;
    while (iSurf<pointsXYonsurface.size() )
    {
      point first = pointsXYonsurface[iSurf];
      point second = first;

      if(expNI == 1)
      {
        second.component(intDir) = minPosInt;
      }
      else
      {
        second = pointsXYonsurface[iSurf+1];
      }
      
      vector dirc =  first - second;
      scalar dist = mag( dirc );
      
      scalar csurf = 0.0;
      scalar integratedUx = 0.0;
      scalar integratedUy = 0.0;
      scalar integratedUCx = 0.0;
      scalar integratedUCy = 0.0;
      scalar Ccup = 0.0;
      if( mag(dirc) > SMALL )
      {
        scalarField variable(K1_);
        scalarField Ux(K1_);
        scalarField Uy(K1_);
        scalarField UCx(K1_);
        scalarField UCy(K1_);
        
        vector dr = 1.0 / static_cast<scalar>(K_) * dirc;

        for(int i=0; i<K1_; i++)
        {
          point samp_point = second + i*dr;
          label cellI = searchEngine.findCell( samp_point );
          if (cellI==-1)
          {
            cellI = searchEngine.findNearestCell( samp_point );
          }
          label faceI = -1;

          vector interp_fieldU =
                  interpolatorU->interpolate(samp_point, cellI, faceI);
          // velocity = 0 on the surface
          if(i==0 || i==K_)
          {
            interp_fieldU = vector::zero;
          }
          
          variable[i] = samp_point.component(intDir);
          Ux[i] = interp_fieldU.component(majDir);
          Uy[i] = interp_fieldU.component(latDir);
          
          scalar interp_fieldC = 
                  interpolatorC->interpolate(samp_point, cellI, faceI);
          UCx[i] = Ux[i] * interp_fieldC;
          UCy[i] = Uy[i] * interp_fieldC;
          if(i==K_)
          {
            csurf = interp_fieldC;
          }
        }
        
        integratedUx = primitive_simpson_integration(variable, Ux);
        integratedUy = primitive_simpson_integration(variable, Uy);
        integratedUCx = primitive_simpson_integration(variable, UCx);
        integratedUCy = primitive_simpson_integration(variable, UCy);
        scalar qSqr = integratedUx*integratedUx+integratedUy*integratedUy;
        if( qSqr > SMALL )
        {
          Ccup = std::sqrt
                 (
                  (integratedUCx*integratedUCx+integratedUCy*integratedUCy)
                  /
                  qSqr
                 );
        }
      }
      
      int ind1 = iSurf + curBlock;
      int prcnt  = 100 * (ind1    ) / N1M1;
      int prcnt1 = 100 * (ind1 - 1) / N1M1;
      if( prcnt%1==0 && prcnt!=prcnt1 )
      {
        Info<<"\r"<< prcnt << "%  "<<flush;
      }
      
      mapXYccup << Ccup << "  ";
      mapXYcsurf << csurf << "  ";
      mapXYqx << integratedUx << "  ";
      mapXYqy << integratedUy << "  ";
      apertureMapFile << dist << "  ";

      iSurf += expNI;
      count++;
      if(count >= NUMBER_OF_COLUMNS)
      { 
        mapXYccup <<"\n";
        mapXYcsurf <<"\n";
        mapXYqx <<"\n";
        mapXYqy <<"\n";
        apertureMapFile <<"\n";
        count=0;
      }
    }  
  }
  /*
  else if( !headerU.headerOk() && headerC.headerOk() )
  {
    FatalError<<"There is no U field"<<nl<<nl<<exit(FatalError);
  }
  else if( headerU.headerOk() && !headerC.headerOk() )
  {
    FatalError<<"There is no C field"<<nl<<nl<<exit(FatalError);
  }
  else
  {
    FatalError<<"There is no U and C field"<<nl<<nl<<exit(FatalError);
  }
  */
}

void Foam::calcTypes::fieldMap2D::write_ccAll
(
  const fvMesh& mesh,
  const Time& runTime
)
{
  typedef GeometricField<vector, fvPatchField, volMesh> fieldU;
  typedef GeometricField<scalar, fvPatchField, volMesh> fieldC;

  IOobject headerU
  (
      "U",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ
  );
  IOobject headerC
  (
      "C",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ
  );
  
  //if (headerU.headerOk() && headerC.headerOk())
  {
    fieldU field_u(headerU, mesh);
    fieldC field_c(headerC, mesh);

    fileName current_file_path_cup =
            "postProcessing/dissolCalc" / runTime.timeName() / "c";
    fileName current_file_path_qx, current_file_path_qy;
    current_file_path_qx =
            "postProcessing/dissolCalc" / runTime.timeName() / "qx";
    current_file_path_qy =
            "postProcessing/dissolCalc" / runTime.timeName() / "qy";
    fileName current_file_path_h =
            "postProcessing/dissolCalc" / runTime.timeName() / "h";
    fileName current_file_path_csurf =
            "postProcessing/dissolCalc" / runTime.timeName() / "csurf";
  
    autoPtr<interpolation<vector> >interpolatorU
    (
      interpolation<vector>::New("cellPoint",field_u)
    );
    autoPtr<interpolation<scalar> >interpolatorC
    (
      interpolation<scalar>::New("cellPoint",field_c)
    );

    ios_base::openmode mode =
            (curNum==0) ? 
              ios_base::out|ios_base::trunc : 
              ios_base::out|ios_base::app;  
    
    OFstreamMod mapXYccup( current_file_path_cup, mode );
    OFstreamMod mapXYqx( current_file_path_qx, mode );
    OFstreamMod mapXYqy( current_file_path_qy, mode );
    OFstreamMod mapXYcsurf( current_file_path_csurf, mode);
    OFstreamMod apertureMapFile( current_file_path_h, mode);
    
    if(curNum==0)
    {
      double dyL = dy * maxPosInt;
      mapXYccup << N_ << "   " << M_ << "   " << runTime.value() 
              << "   " << dx << "   " << dyL << endl;
      mapXYcsurf << N_ << "   " << M_ << "   " << runTime.value() 
              << "   " << dx << "   " << dyL << endl;
      mapXYqx << N_ << "   " << M_ << "   " << runTime.value()
              << "   " << dx << "   " << dyL << endl;
      mapXYqy << N_ << "   " << M_ << "   " << runTime.value()
              << "   " << dx << "   " << dyL << endl;
      apertureMapFile << N_ << "   " << M_ << "   " << runTime.value() 
              << "   " << dx << "   " << dyL << endl;
    }
    
    meshSearch searchEngine(mesh);
        
    int count = 0;
    int iSurf = 0;
    while (iSurf<pointsXYonsurface.size() )
    {
      point first = pointsXYonsurface[iSurf];
      point second = first;

      int totIJK = iSurf + curNum * thisTimeSize;
      label j = totIJK % M1_;
      scalar angle = minPosLat + j*dy;
      
      second.component(latDir) = maxPosInt * Foam::cos(angle);
      second.component(intDir) = maxPosInt * Foam::sin(angle);
      
      vector dirc =  second - first;
      
      scalar dist = mag( dirc );
      
      Info<<"  f: "<<first<<"  s: "<<second<<endl;
      
      scalar csurf = 0.0;
      scalar integratedUx = 0.0;
      scalar integratedUy = 0.0;
      scalar integratedUCx = 0.0;
      scalar integratedUCy = 0.0;
      scalar Ccup = 0.0;
      if( mag(dirc) > SMALL )
      {
        scalarField variable(K1_);
        scalarField Ux(K1_);
        scalarField Uy(K1_);
        scalarField UCx(K1_);
        scalarField UCy(K1_);
        
        vector dr = 1.0 / static_cast<scalar>(K_) * dirc;
        scalar magDr = mag(dr);
        
        vector flDir = vector::zero;
        flDir.component(majDir) = 1.0;
        vector tan = flDir ^ dr;
        tan = tan / mag(tan);

        for(int i=0; i<K1_; i++)
        {
          point samp_point = first + i*dr;
          label cellI = searchEngine.findCell( samp_point );
          if (cellI==-1)
          {
            cellI = searchEngine.findNearestCell( samp_point );
          }
          label faceI = -1;

          vector interp_fieldU =
                  interpolatorU->interpolate(samp_point, cellI, faceI);
          // velocity = 0 on the surface
          if(i==0 || i==K_)
          {
            interp_fieldU = vector::zero;
          }
          
          variable[i] = i * magDr;
          Ux[i] = interp_fieldU & dr / magDr;
          Uy[i] = interp_fieldU & tan;  // tan is already normalized
          
          scalar interp_fieldC = 
                  interpolatorC->interpolate(samp_point, cellI, faceI);
          UCx[i] = Ux[i] * interp_fieldC;
          UCy[i] = Uy[i] * interp_fieldC;
          if(i==K_)
          {
            csurf = interp_fieldC;
          }
        }
        
        integratedUx = primitive_simpson_integration(variable, Ux);
        integratedUy = primitive_simpson_integration(variable, Uy);
        integratedUCx = primitive_simpson_integration(variable, UCx);
        integratedUCy = primitive_simpson_integration(variable, UCy);
        scalar qSqr = integratedUx*integratedUx+integratedUy*integratedUy;
        if( qSqr > SMALL )
        {
          Ccup = std::sqrt
                 (
                  (integratedUCx*integratedUCx+integratedUCy*integratedUCy)
                  /
                  qSqr
                 );
        }
      }
      
      int ind1 = iSurf + curBlock;
      int prcnt  = 100 * (ind1    ) / N1M1;
      int prcnt1 = 100 * (ind1 - 1) / N1M1;
      if( prcnt%1==0 && prcnt!=prcnt1 )
      {
        Info<<"\r"<< prcnt << "%  "<<flush;
      }
      
      mapXYccup << Ccup << "  ";
      mapXYcsurf << csurf << "  ";
      mapXYqx << integratedUx << "  ";
      mapXYqy << integratedUy << "  ";
      apertureMapFile << dist << "  ";

      iSurf += expNI;
      count++;
      if(count >= NUMBER_OF_COLUMNS)
      { 
        mapXYccup <<"\n";
        mapXYcsurf <<"\n";
        mapXYqx <<"\n";
        mapXYqy <<"\n";
        apertureMapFile <<"\n";
        count=0;
      }
    }  
  }
  /*
  else if( !headerU.headerOk() && headerC.headerOk() )
  {
    FatalError<<"There is no U field"<<nl<<nl<<exit(FatalError);
  }
  else if( headerU.headerOk() && !headerC.headerOk() )
  {
    FatalError<<"There is no C field"<<nl<<nl<<exit(FatalError);
  }
  else
  {
    FatalError<<"There is no U and C field"<<nl<<nl<<exit(FatalError);
  }
  */
}



// ************************************************************************* //

