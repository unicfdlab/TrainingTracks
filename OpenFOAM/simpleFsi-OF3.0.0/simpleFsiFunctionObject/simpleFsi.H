/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Class
    Foam::simpleFsi

Group
    grpsimpleFsiFunctionObjects

Description

SeeAlso


SourceFiles


\*---------------------------------------------------------------------------*/

#ifndef simpleFsi_H
#define simpleFsi_H

#include "functionObjectFile.H"
#include "coordinateSystem.H"
#include "coordinateSystems.H"
#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "HashSet.H"
#include "Tuple2.H"
#include "OFstream.H"
#include "Switch.H"
#include "writer.H"
#include "forces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class simpleFsi Declaration
\*---------------------------------------------------------------------------*/

class simpleFsi
:
    public forces
{

protected:

        //- cylinder mass
        scalar M_;
        
        //- damping coefficient
        scalar C_;
        
        //- rigidity coefficient
        scalar K_;
        
        //- ratio of real cylinder length to depth of domain
        scalar R_;
        
        //- current state of system (y,Vy)
        Pair<scalar> Y_; 
        
        //- old state of system (y, Vy)
        Pair<scalar> Yold_;
        
        //-
        autoPtr<OFstream> yD_;

protected:

        //- Disallow default bitwise copy construct
        simpleFsi(const simpleFsi&);

        //- Disallow default bitwise assignment
        void operator=(const simpleFsi&);


public:

    //- Runtime type information
    TypeName("simpleFsi");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        simpleFsi
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false,
            const bool readFields = true
        );


        //- Construct from components
        simpleFsi
        (
            const word& name,
            const objectRegistry&,
            const labelHashSet& patchSet,
            const word& pName,
            const word& UName,
            const word& rhoName,
            const scalar rhoInf,
            const scalar pRef,
            const coordinateSystem& coordSys
        );


    //- Destructor
    virtual ~simpleFsi();


    // Member Functions

        //- Read the simpleFsi data
        virtual void read(const dictionary&);

        //- Write the simpleFsi
        virtual void write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
