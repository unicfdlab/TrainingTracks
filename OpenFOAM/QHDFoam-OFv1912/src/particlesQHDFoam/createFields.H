Info<< "Reading thermophysical properties\n" << endl;

autoPtr<rhoQGDThermo> pThermo
(
    rhoQGDThermo::New(mesh)
);
rhoQGDThermo& thermo = pThermo();
thermo.correct();

volScalarField& e = thermo.he();

volScalarField& p = thermo.p();
const volScalarField& T = thermo.T();
const surfaceScalarField& hQGDf = thermo.hQGDf();
const surfaceScalarField& tauQGDf = thermo.tauQGDf();

Info << "Thermo corrected" << endl;

autoPtr<compressible::turbulenceModel> turbulence;

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

volVectorField W
(
    IOobject
    (
        "W",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    U
);




Info<< "\nReading gravitationalProperties" << endl;
  
IOdictionary gravitationalProperties
(
    IOobject
    (
        "gravitationalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const dimensionedVector g(gravitationalProperties.lookup("g"));

dimensionedScalar beta
(
    "beta",
    dimless/dimTemperature,
    thermo.subDict("mixture").subDict("transport")
);


surfaceScalarField phiu
(
    "phiu",
    mesh.Sf() & linearInterpolate(U)
);

surfaceScalarField  phiwo
(
    "phiwo",
    mesh.Sf() & linearInterpolate(W)
);

surfaceScalarField phi
(
    "phi",
    mesh.Sf() & (linearInterpolate(U) - linearInterpolate(W))
);

volVectorField BdFrc
(
    "BdFrc",
    T*g*beta
);

Switch implicitDiffusion (thermo.implicitDiffusion());

surfaceScalarField phiRhof
(
    "phiRhof",
    linearInterpolate(rho)*phi
);

Info<< "Creating turbulence model\n" << endl;
turbulence.reset
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phiRhof,
        thermo
    ).ptr()
);

#include "createIcoZeroSources.H"

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, thermo.subDict("QGD"), pRefCell, pRefValue);

SLGThermo slgThermo(mesh, thermo);

#include "createClouds.H"

//
//END-OF-FILE
//

