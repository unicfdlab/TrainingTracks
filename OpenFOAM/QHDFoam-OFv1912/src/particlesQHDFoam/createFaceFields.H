    /*
     *
     * Linear interpolation of fields from volumes to face centers
     *
     */
    
    // Density
    surfaceScalarField rhof
    (
        "rhof",
        linearInterpolate(rho)
    );
    
    // Velocity
    surfaceVectorField Uf
    (
        "Uf",
        linearInterpolate(U)
    );

    surfaceScalarField pf
    (
        "pf",
	linearInterpolate(p)
    );
    
    //
    surfaceVectorField Wf
    (
        "Wf",
        linearInterpolate(W)
    );
    
    surfaceVectorField BdFrcf
    (
        "BdFrcf",
        linearInterpolate(BdFrc)
    );
    
    surfaceScalarField Tf
    (
        "Tf",
	linearInterpolate(T)
    );
    
    // Heat capacity at constant pressure
    surfaceScalarField Cpf
    (
        "Cpf",
        linearInterpolate(thermo.Cp())
    );
    
    surfaceScalarField muf
    (
        "muf",
        linearInterpolate(thermo.mu()*0.0)
    );
    
    /*
     *
     * QGD coefficients
     *
     */
    // OpenFOAM delta coeffs can be negative for coupled patches
    surfaceScalarField taubyrhof
    (
        "taubyrhof",
        tauQGDf/rhof
    );
    
    surfaceScalarField alphaf
    (
        "alphaf",
        linearInterpolate(thermo.alpha())
    );


