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
    Foam::helloWorld

Group
    grphelloWorldFunctionObjects

Description
    

SeeAlso
    

SourceFiles


\*---------------------------------------------------------------------------*/

#ifndef helloWorld_H
#define helloWorld_H

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
                           Class helloWorld Declaration
\*---------------------------------------------------------------------------*/

class helloWorld
:
    public forces
{

protected:
    
    //- Disallow default bitwise copy construct
    helloWorld(const helloWorld&);

    //- Disallow default bitwise assignment
    void operator=(const helloWorld&);

public:

    //- Runtime type information
    TypeName("helloWorld");

    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        helloWorld
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false,
            const bool readFields = true
        );


    //- Destructor
        virtual ~helloWorld();


    // Member Functions

        //- Read the helloWorld data
        virtual void read(const dictionary&);

        //- Write the helloWorld
        virtual void write();
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
