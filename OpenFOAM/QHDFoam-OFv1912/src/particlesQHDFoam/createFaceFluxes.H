
/*
  Field for debugging parallel execution
  and tau-terms evaluation
*/
/*
volScalarField cellNo
(
    IOobject
    (
        "cellNo",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rho*0.0
);

forAll(cellNo, celli)
{
    cellNo[celli] = mesh.C()[celli].z();
}
forAll(cellNo.boundaryField(), patchi)
{
    forAll(cellNo.boundaryField()[patchi], facei)
    {
        cellNo.boundaryFieldRef()[patchi][facei] = mesh.Cf().boundaryField()[patchi][facei].z();
    }
}
surfaceVectorField cellNoGradf
(
    "cellNoGradf", fvsc::grad(cellNo)
);

cellNo.write();
cellNoGradf.write();
return 0;
*/
//Gradients and divergence
//---------Start---------


surfaceVectorField gradPf 
(
    "gradPf", fvsc::grad(p)
);

surfaceVectorField gradTf
(
    "gradTf", fvsc::grad(T)
);

surfaceTensorField gradUf
(
    "gradUf",
    fvsc::grad(U)
);

surfaceTensorField gradWf
(
    "gradWf",
    fvsc::grad(W)
);

surfaceVectorField gradRhof
(
    "gradRhof",
    fvsc::grad(rho)
);

surfaceScalarField divUf
(
    "divUf",
    tr(gradUf)
);

//---------End---------
//P equation fluxes
//---------Start---------
//---------End---------

// Fluxes for momentum balance equation
//---------Start---------
surfaceVectorField phiUf
(
    "phiUf",
    phi * Uf + (mesh.Sf() & (Wf * Wf))
);

surfaceScalarField Pbyrho
(
    "Pbyrho",
    pf/rhof
);

surfaceScalarField phiTf
(
    "phiTf",
    phi*Tf
);

surfaceVectorField phiPi
(
    "phiPi",
    muf/rhof * mesh.Sf() & (gradUf+gradUf.T())
);

//---------End---------

// Fluxes for T balance equation
//---------Start---------
surfaceScalarField phiUmWT
(
    "phiUmWT",
    ((Uf-Wf)*Tf) & mesh.Sf()
);

surfaceScalarField Hif
(
    "Hif",
    alphaf/rhof
);
Info << "Creating Face Fluxes is done"<< endl;

//---------End---------
