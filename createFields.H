#include "createRDeltaT.H"
Info<< "Reading thermophysical properties\n" << endl;

autoPtr<psiThermo> pThermo
(
    psiThermo::New(mesh)
);
psiThermo& thermo = pThermo();

volScalarField& e = thermo.he();

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// alpha to be defined in transport properties
scalar alphaMType = transportProperties.lookupOrDefault<scalar>("alphaMType",1);
scalar alphaM = transportProperties.lookupOrDefault<scalar>("alphaM",0);
Info<< "\n alphaM: "<< alphaM << endl;
Info<< "\n alphaMType: "<< alphaMType << endl;

// Initial mass-velocity field to be read and auto-write
Info<< "Reading field U\n" << endl;
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


//Density field
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
    thermo.rho()   //calculates density from thermo physical model
);


//rho*U
volVectorField rhoU
(
    IOobject
    (
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U
);


volScalarField rhoE
(
    IOobject
    (
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*(e + 0.5*magSqr(U))
);


surfaceScalarField pos
(
    IOobject
    (
        "pos",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("pos", dimless, 1.0)
);

surfaceScalarField neg
(
    IOobject
    (
        "neg",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("neg", dimless, -1.0)
);

// flux of rho*U
surfaceScalarField phi("phi", fvc::flux(rhoU));


Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

 volScalarField lnr
 (
    log(rho*dimensionedScalar("one", dimless/dimDensity, 1))
 );

// kappa_m calculation 

volScalarField Km
(
    IOobject
    (
        "Km",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
   alphaM*turbulence->muEff()/rho
);

if(alphaMType==2)
{ 
	Km = alphaM*turbulence->alphaEff()/rho;
}

volVectorField Jv
(
    IOobject
    (
        "Jv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    Km*fvc::grad(lnr)
);

volVectorField Uns
(
    IOobject
    (
        "Uns",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    U-Jv
);

        
Info<< "\n Max(Km)\n"<< max(Km) << endl;
Info<< "\n Min(Km)\n"<< min(Km) << endl;