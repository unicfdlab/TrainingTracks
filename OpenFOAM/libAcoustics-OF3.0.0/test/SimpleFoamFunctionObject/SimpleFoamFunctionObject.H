/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::SimpleFoamFunctionObject

Description
    Simple example of the Foam::functionObject API usage

SourceFiles
    SimpleFoamFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef SimpleFoamFunctionObject_H
#define SimpleFoamFunctionObject_H

#include "functionObject.H"
#include "dictionary.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class SimpleFoamFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class SimpleFoamFunctionObject
:
    public functionObject
{
    // Private data

        //- Reference to the time database
        const Time& time_;

        //- Input dictionary
        dictionary dict_;
        
        //- greetings address
        word helloAddress_;

    // Private Member Functions

        //- Read relevant dictionary entries
        void readDict();

        //- Disallow default bitwise copy construct
        SimpleFoamFunctionObject(const SimpleFoamFunctionObject&);

        //- Disallow default bitwise assignment
        void operator=(const SimpleFoamFunctionObject&);


public:

    //- Runtime type information
    TypeName("SimpleFoamFunctionObject");


    // Constructors

        //- Construct from components
        SimpleFoamFunctionObject
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    // Member Functions

        // Access

            //- Return time database
            virtual const Time& time() const
            {
                return time_;
            }

            //- Return the input dictionary
            virtual const dictionary& dict() const
            {
                return dict_;
            }


        // Function object control


            //- Called at the start of the time-loop
            virtual bool start();

            //- Called at each ++ or += of the time-loop
            virtual bool execute(bool);

            //- Called when Time::run() determines that the time-loop exits
            virtual bool end();

            //- Read and set the function object if its data have changed
            virtual bool read(const dictionary&);
            
	// In OpenFOAM-2.2
	 //- called whesh topology was changed
	 virtual void updateMesh(const Foam::mapPolyMesh&)
	 {
	 };
	 
	 //- called when mesh points were moved
	 virtual void movePoints(const Foam::polyMesh&)
	 {
	 }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
