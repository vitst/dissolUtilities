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
    
    Write a description

SourceFiles
    fieldMap2d.C

\*---------------------------------------------------------------------------*/

#ifndef fieldMap2d_H
#define fieldMap2d_H

#include "calcType.H"
#include <map>
#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace calcTypes
{

/*---------------------------------------------------------------------------*\
                           Class fieldMap2d Declaration
\*---------------------------------------------------------------------------*/

class fieldMap2d
:
    public calcType
{
    // Private Member Functions
  
        const static int maxNumProcPoints = 1000000; // max size of the pattern
        
        int thisTimeSize;
        int totNumLoop;
        
        int sizeAA;
        
        int curNum;
        
        int curBlock;
        
        int curEnd;
        
        //- Disallow default bitwise copy construct
        fieldMap2d(const fieldMap2d&);

        //- Disallow default bitwise assignment
        void operator=(const fieldMap2d&);
        
        word processingType_;
        label N_, M_, K_;   // N - major, M - lateral
        label N1_, M1_, K1_;
        
        int N1M1;
        
        word patchName;
        word geometry;
        
        int expNI;  //   expected number of intersections
        
        //scalar maxPosX, minPosX; // maximum and minimum in X direction
        //scalar maxPosY, minPosY; // maximum and minimum in Y direction
        //scalar maxPosZ, minPosZ; // maximum and minimum in Z direction
        
        // directions
        label majDir, latDir, intDir;
        
        scalar maxPosMaj, minPosMaj;
        scalar maxPosLat, minPosLat;
        // in case of one surface minPosIntegr is a surface confining 
        // a reservoir instead of second wall
        scalar maxPosInt, minPosInt;
        
        scalar dx, dy; // dx - major, dy - lateral
        
        pointField pointsXYonsurface;
        
protected:

    // Member Functions

        // Calculation routines

            //- Initialise - typically setting static variables,
            //  e.g. command line arguments
            virtual void init();

            //- Pre-time loop calculations
            virtual void preCalc
            (
                const argList& args,
                const Time& runTime,
                const fvMesh& mesh
            );

            //- Time loop calculations
            virtual void calc
            (
                const argList& args,
                const Time& runTime,
                const fvMesh& mesh
            );

            void build_surface_points
            (
              const fvMesh& mesh
            );
            
            void build_surface_pointsCC
            (
              const fvMesh& mesh
            );
            
            scalar primitive_simpson_integration
            (
              scalarField& x,
              scalarField& y  
            );

        // I-O
            //- Write component fields
            void write_temp
            (
              const fvMesh& mesh,
              const Time& runTime
            );
            void write_csurf
            (
              const fvMesh& mesh,
              const Time& runTime
            );
            void write_h
            (
              const fvMesh& mesh,
              const Time& runTime
            );
            void write_q
            (
              const fvMesh& mesh,
              const Time& runTime
            );
            void write_Ccup
            (
              const fvMesh& mesh,
              const Time& runTime
            );

            void write_all
            (
              const fvMesh& mesh,
              const Time& runTime
            );
            void write_surf
            (
              const fvMesh& mesh,
              const Time& runTime
            );
            void write_int
            (
              const fvMesh& mesh,
              const Time& runTime
            );
            
            void write_ccAll
            (
              const fvMesh& mesh,
              const Time& runTime
            );
public:

    //- Runtime type information
    TypeName("fieldMap2d");

    // Constructor
        fieldMap2d();

    //- Destructor
        virtual ~fieldMap2d();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace calcTypes
} // End namespace Foam

#endif

// ************************************************************************* //
