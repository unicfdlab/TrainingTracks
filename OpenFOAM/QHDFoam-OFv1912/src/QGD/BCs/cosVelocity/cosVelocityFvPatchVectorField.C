/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) 2016-2019 ISP RAS (www.ispras.ru) UniCFD Group (www.unicfd.ru)
-------------------------------------------------------------------------------
License
    This file is part of QGDsolver library, based on OpenFOAM+.

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
\*---------------------------------------------------------------------------*/

#include "cosVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cosVelocityFvPatchVectorField::cosVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(p, iF),
    A_(0.0, 0.0, 0.0),
    H_(1.0),
    Hdirection_(0.0, 0.0, 0.0),
    omega0_(0.0),
    phi0_(0.0),
    minZ_(0.0)
{
    this->refValue() = pTraits<vector>::zero;
    this->refGrad() = pTraits<vector>::zero;
    this->valueFraction() = 0.0;
    this->setPatchVelocities(this->refValue());
    
    Info << "Executing cosVelocityFvPatchVectorField(fvPatch, DimensionedField) " << endl;
}


Foam::cosVelocityFvPatchVectorField::cosVelocityFvPatchVectorField
(
    const cosVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<vector>(ptf, p, iF, mapper),
    A_(ptf.A_),
    H_(ptf.H_),
    Hdirection_(ptf.Hdirection_),
    omega0_(ptf.omega0_),
    phi0_(ptf.phi0_),
    minZ_(ptf.minZ_)
{
    this->setPatchVelocities(this->refValue());
    mixedFvPatchVectorField::updateCoeffs();

    Info << "Executing cosVelocityFvPatchVectorField(cosVelocity, fvPatch, DimensionedField, fvPatchFieldMapper) " << endl;
}


Foam::cosVelocityFvPatchVectorField::cosVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<vector>(p, iF)
{
    
    dict.lookup("A") >> A_;
    dict.lookup("H") >> H_;
    dict.lookup("Hdirection") >> Hdirection_;
    dict.lookup("omega0") >> omega0_;
    dict.lookup("phi0") >> phi0_;
    dict.lookup("minZ") >> minZ_;

    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            Field<vector>("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<vector>::operator=(this->patchInternalField());
    }

    this->refValue() = pTraits<vector>::zero;
    this->refGrad() = pTraits<vector>::zero;
    this->valueFraction() = 1.0;
    this->setPatchVelocities(this->refValue());
    mixedFvPatchVectorField::updateCoeffs();

    Info << "Executing cosVelocityFvPatchVectorField(fvPatch, DimensionedField, dictionary) " << endl;
}


Foam::cosVelocityFvPatchVectorField::cosVelocityFvPatchVectorField
(
    const cosVelocityFvPatchVectorField& ptpsf
)
:
    mixedFvPatchVectorField(ptpsf),
    A_(ptpsf.A_),
    H_(ptpsf.H_),
    Hdirection_(ptpsf.Hdirection_),
    omega0_(ptpsf.omega0_),
    phi0_(ptpsf.phi0_),
    minZ_(ptpsf.minZ_)
{
    this->setPatchVelocities(this->refValue());
    Info << "Executing cosVelocityFvPatchVectorField(const cosVelocity&) " << endl;
}


Foam::cosVelocityFvPatchVectorField::cosVelocityFvPatchVectorField
(
    const cosVelocityFvPatchVectorField& ptpsf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(ptpsf, iF),
    A_(ptpsf.A_),
    H_(ptpsf.H_),
    Hdirection_(ptpsf.Hdirection_),
    omega0_(ptpsf.omega0_),
    phi0_(ptpsf.phi0_),
    minZ_(ptpsf.minZ_)
{
    this->setPatchVelocities(this->refValue());

    Info << "Executing cosVelocityFvPatchVectorField(const cosVelocity, const DimensionedField) " << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cosVelocityFvPatchVectorField::autoMap(const fvPatchFieldMapper& m)
{
    mixedFvPatchField<vector>::autoMap(m);

    Info << "Executing autoMap " << endl;
}

void Foam::cosVelocityFvPatchVectorField::rmap(const fvPatchVectorField& pf, const labelList& ll)
{
    mixedFvPatchField<vector>::rmap(pf, ll);
    Info << "Executing rmap " << endl;
}

void Foam::cosVelocityFvPatchVectorField::setPatchVelocities (vectorField& pV)
{
    scalarField z = (this->patch().Cf() & Hdirection_) - minZ_;
    
    forAll(pV, iFace)
    {
	scalar t = this->patch().boundaryMesh().mesh().time().value();
	this->refGrad()[iFace] = vector (0.0, 0.0, 0.0);
	this->valueFraction()[iFace] = 1.0;
	pV[iFace] = 
	    A_ * cos(Foam::constant::mathematical::pi * z[iFace] / H_) * (-omega0_) * sin(omega0_ * t + phi0_);
	if (z[iFace] < 0.0 || z[iFace] > H_)
	{
	    pV[iFace] = A_*0.0;
	}
    }
}

void Foam::cosVelocityFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    
    this->setPatchVelocities(this->refValue());
    
    Info << "max/min velocity: " << gMax(this->refValue()) << "/" << gMin(this->refValue()) << endl;

    mixedFvPatchField<vector>::updateCoeffs();
}


void Foam::cosVelocityFvPatchVectorField::write(Ostream& os) const
{

    mixedFvPatchVectorField::write(os);

    os.writeKeyword("A") << A_ << token::END_STATEMENT << nl;
        
    os.writeKeyword("H") << H_ << token::END_STATEMENT << nl;
    
    os.writeKeyword("minZ") << minZ_ << token::END_STATEMENT << nl;
    
    os.writeKeyword("Hdirection") << Hdirection_ << token::END_STATEMENT << nl;
    
    os.writeKeyword("omega0") << omega0_ << token::END_STATEMENT << nl;
    
    os.writeKeyword("phi0") << phi0_ << token::END_STATEMENT << nl;
    
}

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        cosVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
