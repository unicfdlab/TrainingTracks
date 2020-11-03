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
    Foam::weaklyCoupled3DoFFsi

Group
    grpweaklyCoupled3DoFFsiFunctionObjects

Description
    

SeeAlso
    

SourceFiles


\*---------------------------------------------------------------------------*/

#ifndef weaklyCoupled3DoFFsi_H
#define weaklyCoupled3DoFFsi_H


#include "forces.H"
#include "volFieldsFwd.H"
#include "Tuple2.H"
#include "OFstream.H"
#include "coordinateSystem.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class weaklyCoupled3DoFFsi Declaration
\*---------------------------------------------------------------------------*/

class weaklyCoupled3DoFFsi
:
    public forces
{

protected:

        //- cylinder mass
        scalar M_;
        
        //- damping coefficient
        vector C_;
        
        //- rigidity coefficient
        vector K_;
        
        //- ratio of real cylinder length to depth of domain
        scalar R_;

        //- maximum amplitude of displacement
        vector Ymax_;

        //- true if after restart data should be appended to previous log 
        //  false if log should be overwritten
        bool append_;

        //- current state of system (y,Vy)
        Pair<vector> Y_; 
        
        //- old state of system (y, Vy)
        Pair<vector> Yold_;
        
        //-
        autoPtr<coordinateSystem> coordSys_;

protected:

        //- Disallow default bitwise copy construct
        weaklyCoupled3DoFFsi(const weaklyCoupled3DoFFsi&);

        //- Disallow default bitwise assignment
        void operator=(const weaklyCoupled3DoFFsi&);
        
        //-
        void createFsiOutFile(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("weaklyCoupled3DoFFsi");


    // Constructors

        //- Construct from given Time and dictionary.
        weaklyCoupled3DoFFsi
        (
            const word& name,
            const Time&,
            const dictionary&
        );
        
        //- Construct from given objectRegistry and dictionary.
        weaklyCoupled3DoFFsi
        (
            const word& name,
            const objectRegistry&,
            const dictionary&
        );

    //- Destructor
    virtual ~weaklyCoupled3DoFFsi();


    // Member Functions

        //- Read the weaklyCoupled3DoFFsi data
        virtual bool read(const dictionary&);
        
        //-
        virtual bool execute();
        
        //- Write the weaklyCoupled3DoFFsi
        virtual bool write();
        
        //- Distributes displacements between slave processes
        // and sets cellDisplacement field Y component on patch
        void setDisplacements(volVectorField& yDispl);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
