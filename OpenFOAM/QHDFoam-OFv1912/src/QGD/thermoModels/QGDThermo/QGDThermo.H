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
    This file is part of QGDsolver, based on OpenFOAM library.

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
Group
     QGDThermo
Class
    Foam::QGDThermo

Description
    Abstract base class for classes implementing thermophysical properties
    of gases and fluids governed by regularized equations (QGD & QHD)

SourceFiles
    QGDThermo.C

\*---------------------------------------------------------------------------*/

#ifndef QGDTHERMO_H
#define QGDTHERMO_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "Switch.H"

namespace Foam
{
class fvMesh;
class dictionary;
namespace qgd
{
  class QGDCoeffs;
}

class QGDThermo
{

  //-
  const fvMesh& mesh_;

  //-
  const dictionary& dict_;

  //-
  autoPtr<Foam::qgd::QGDCoeffs> qgdCoeffsPtr_;

  //-
  Switch implicitDiffusion_;

private:

  //-
  QGDThermo();

  //-
  QGDThermo(const QGDThermo& );

protected:

  //-
  qgd::QGDCoeffs& qgdCoeffs();

  //-
  void correctQGD(volScalarField& mu, volScalarField& alphau);

public:

  TypeName ("QGDThermo");

  //- Construct from mesh and dictionary
  QGDThermo(const fvMesh& mesh, const dictionary& dict);

  virtual ~QGDThermo();

  //-
  virtual const volScalarField& c() const = 0;

  //-
  virtual const volScalarField& p() const = 0;

  //-
  virtual tmp<volScalarField> rho() const = 0;

  //-
  virtual tmp<volScalarField> mu() const = 0;

  //-
  const volScalarField& hQGD() const;

  //-
  const volScalarField& tauQGD() const;

  //-
  const surfaceScalarField& hQGDf() const;

  //-
  const surfaceScalarField& tauQGDf() const;

  //-
  const volScalarField& muQGD() const;

  //-
  const volScalarField& alphauQGD() const;

  //-
  const volScalarField& aQGD() const;

  //-
  const surfaceScalarField& aQGDf() const;

  //-
  Switch implicitDiffusion() const;

  //-
  virtual bool read();

};
}

#endif

//
//END-OF-FILE
//
