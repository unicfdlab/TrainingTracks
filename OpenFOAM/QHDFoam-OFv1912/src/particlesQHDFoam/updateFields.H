/*
 *
 * Gradients at faces
 *
 */
gradUf = fvsc::grad(U);

gradWf = fvsc::grad(W);

gradTf = fvsc::grad(T);

/*
 *
 * QGD coefficients now updated in thermo library
 *
 */
muf    = qgdInterpolate(turbulence->muEff());
alphaf = qgdInterpolate(thermo.alpha());

/*
 *
 * Linear interpolation of fields from volumes to face centers
 *
 */

// Density
rhof = qgdInterpolate(rho);

// Velocity
Uf   = qgdInterpolate(U);
    
// T
Tf   = qgdInterpolate(T);

//
BdFrc  = -beta*T*g;
BdFrcf = qgdInterpolate(BdFrc);

//Hi
Hif = alphaf/rhof;

// Heat capacity at constant pressure
Cpf = qgdInterpolate(thermo.Cp());




