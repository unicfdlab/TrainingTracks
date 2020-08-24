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

Group
    grpPsiQGDReactionThermo

\*---------------------------------------------------------------------------*/

#include "makeReactionThermo.H"

#include "psiQGDReactionThermo.H"
#include "hePsiQGDThermo.H"

#include "specie.H"
#include "perfectGas.H"
#include "hConstThermo.H"
#include "janafThermo.H"
#include "sensibleEnthalpy.H"
#include "thermo.H"
#include "constTransport.H"
#include "sutherlandTransport.H"

#include "homogeneousMixture.H"
#include "inhomogeneousMixture.H"
#include "veryInhomogeneousMixture.H"
#include "multiComponentMixture.H"
#include "reactingMixture.H"
#include "singleStepReactingMixture.H"
#include "singleComponentMixture.H"

#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// constTransport, hConstThermo

makeReactionThermos
(
    psiQGDThermo,
    psiQGDReactionThermo,
    hePsiQGDThermo,
    homogeneousMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeReactionThermos
(
    psiQGDThermo,
    psiQGDReactionThermo,
    hePsiQGDThermo,
    inhomogeneousMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeReactionThermos
(
    psiQGDThermo,
    psiQGDReactionThermo,
    hePsiQGDThermo,
    veryInhomogeneousMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);


// sutherlandTransport, hConstThermo

makeReactionThermos
(
    psiQGDThermo,
    psiQGDReactionThermo,
    hePsiQGDThermo,
    homogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeReactionThermos
(
    psiQGDThermo,
    psiQGDReactionThermo,
    hePsiQGDThermo,
    inhomogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeReactionThermos
(
    psiQGDThermo,
    psiQGDReactionThermo,
    hePsiQGDThermo,
    veryInhomogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);


// sutherlandTransport, janafThermo

makeReactionThermos
(
    psiQGDThermo,
    psiQGDReactionThermo,
    hePsiQGDThermo,
    homogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    perfectGas,
    specie
);

makeReactionThermos
(
    psiQGDThermo,
    psiQGDReactionThermo,
    hePsiQGDThermo,
    inhomogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    perfectGas,
    specie
);

makeReactionThermos
(
    psiQGDThermo,
    psiQGDReactionThermo,
    hePsiQGDThermo,
    veryInhomogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    perfectGas,
    specie
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Multi-component thermo for sensible enthalpy

makeThermoPhysicsReactionThermos
(
    psiQGDThermo,
    psiQGDReactionThermo,
    hePsiQGDThermo,
    multiComponentMixture,
    constGasHThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    psiQGDThermo,
    psiQGDReactionThermo,
    hePsiQGDThermo,
    multiComponentMixture,
    gasHThermoPhysics
);


// Multi-component thermo for internal energy

makeThermoPhysicsReactionThermos
(
    psiQGDThermo,
    psiQGDReactionThermo,
    hePsiQGDThermo,
    multiComponentMixture,
    constGasEThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    psiQGDThermo,
    psiQGDReactionThermo,
    hePsiQGDThermo,
    multiComponentMixture,
    gasEThermoPhysics
);


// Reaction thermo for sensible enthalpy

makeThermoPhysicsReactionThermos
(
    psiQGDThermo,
    psiQGDReactionThermo,
    hePsiQGDThermo,
    reactingMixture,
    constGasHThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    psiQGDThermo,
    psiQGDReactionThermo,
    hePsiQGDThermo,
    reactingMixture,
    gasHThermoPhysics
);


// Single-step reaction thermo for sensible enthalpy

makeThermoPhysicsReactionThermos
(
    psiQGDThermo,
    psiQGDReactionThermo,
    hePsiQGDThermo,
    singleStepReactingMixture,
    gasHThermoPhysics
);


// Reaction thermo for internal energy

makeThermoPhysicsReactionThermos
(
    psiQGDThermo,
    psiQGDReactionThermo,
    hePsiQGDThermo,
    reactingMixture,
    constGasEThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    psiQGDThermo,
    psiQGDReactionThermo,
    hePsiQGDThermo,
    reactingMixture,
    gasEThermoPhysics
);


// Single-step reaction thermo for internal energy

makeThermoPhysicsReactionThermos
(
    psiQGDThermo,
    psiQGDReactionThermo,
    hePsiQGDThermo,
    singleStepReactingMixture,
    gasEThermoPhysics
);


// Single-component thermo for sensible enthalpy

makeThermoPhysicsReactionThermo
(
    psiQGDReactionThermo,
    hePsiQGDThermo,
    singleComponentMixture,
    constGasHThermoPhysics
);

makeThermoPhysicsReactionThermo
(
    psiQGDReactionThermo,
    hePsiQGDThermo,
    singleComponentMixture,
    gasHThermoPhysics
);


// Single-component thermo for internal energy

makeThermoPhysicsReactionThermo
(
    psiQGDReactionThermo,
    hePsiQGDThermo,
    singleComponentMixture,
    constGasEThermoPhysics
);

makeThermoPhysicsReactionThermo
(
    psiQGDReactionThermo,
    hePsiQGDThermo,
    singleComponentMixture,
    gasEThermoPhysics
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
