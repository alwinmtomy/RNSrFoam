/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Application
    RNSFoam

Description
    Density-based compressible flow solver based on central-upwind schemes of
    Kurganov and Tadmor.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "RNSfixedRhoFvPatchScalarField.H"
#include "directionInterpolate.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createTimeControls.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * **  * * * * * * * //

    #include "readFluxScheme.H"

    dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;
  
    Info<< "\nStarting time loop\n" << endl;
          
     while (runTime.run())
    {
        // --- Directed interpolation of primitive fields onto faces

        surfaceScalarField rho_pos(interpolate(rho, pos));
        surfaceScalarField rho_neg(interpolate(rho, neg));

        surfaceVectorField rhoU_pos(interpolate(rhoU, pos, U.name()));
        surfaceVectorField rhoU_neg(interpolate(rhoU, neg, U.name()));

        volScalarField rPsi("rPsi", 1.0/psi); //psi from psiThermo.H [ p * psi = rho ; psi = comprresibility]
        surfaceScalarField rPsi_pos(interpolate(rPsi, pos, T.name())); 
        surfaceScalarField rPsi_neg(interpolate(rPsi, neg, T.name()));

        surfaceScalarField e_pos(interpolate(e, pos, T.name()));
        surfaceScalarField e_neg(interpolate(e, neg, T.name()));

        surfaceVectorField U_pos("U_pos", rhoU_pos/rho_pos);
        surfaceVectorField U_neg("U_neg", rhoU_neg/rho_neg);

        surfaceScalarField p_pos("p_pos", rho_pos*rPsi_pos);
        surfaceScalarField p_neg("p_neg", rho_neg*rPsi_neg);

        surfaceScalarField phiv_pos("phiv_pos", U_pos & mesh.Sf());
        surfaceScalarField phiv_neg("phiv_neg", U_neg & mesh.Sf());

        volScalarField c("c", sqrt(thermo.Cp()/thermo.Cv()*rPsi));
        surfaceScalarField cSf_pos
        (
            "cSf_pos",
            interpolate(c, pos, T.name())*mesh.magSf()  
        );
        surfaceScalarField cSf_neg
        (
            "cSf_neg",
            interpolate(c, neg, T.name())*mesh.magSf()
        );

        surfaceScalarField ap
        (
            "ap",
            max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero)
        );
        surfaceScalarField am
        (
            "am",
            min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero)
        );

        surfaceScalarField a_pos("a_pos", ap/(ap - am));

        surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));

        surfaceScalarField aSf("aSf", am*a_pos);

        if (fluxScheme == "Tadmor")
        {
            aSf = -0.5*amaxSf;
            a_pos = 0.5;
        }

        surfaceScalarField a_neg("a_neg", 1.0 - a_pos);

        phiv_pos *= a_pos;
        phiv_neg *= a_neg;

        surfaceScalarField aphiv_pos("aphiv_pos", phiv_pos - aSf);
        surfaceScalarField aphiv_neg("aphiv_neg", phiv_neg + aSf);

        // Reuse amaxSf for the maximum positive and negative fluxes
        // estimated by the central scheme
        amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));

        #include "centralCourantNo.H"
        #include "readTimeControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "setDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        phi = aphiv_pos*rho_pos + aphiv_neg*rho_neg;

        surfaceVectorField phiUp    //u.rhoU+pI
        (
            (aphiv_pos*rhoU_pos + aphiv_neg*rhoU_neg)
          + (a_pos*p_pos + a_neg*p_neg)*mesh.Sf()
        );

        surfaceScalarField phiEp    //u(E+p)
        (
            "phiEp",
            aphiv_pos*(rho_pos*(e_pos + 0.5*magSqr(U_pos)) + p_pos)
          + aphiv_neg*(rho_neg*(e_neg + 0.5*magSqr(U_neg)) + p_neg)
          + aSf*p_pos - aSf*p_neg
        );

        volScalarField muEff("muEff", turbulence->muEff()); 

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

        volScalarField lnr (
        	IOobject
		    (
		        "lnr",
		        runTime.timeName(),
		        mesh,
		        IOobject::NO_READ,
		        IOobject::AUTO_WRITE
		    ),
		    log(rho*dimensionedScalar("one", dimless/dimDensity, 1)) );

        //volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U))));  //part of stress tensor NS
    	// Foam::T() <-- transposes matrix 
        //stress tensor  modified NS
        volTensorField tau_v1
        (   "tau_v1",
            muEff*dev2(Foam::T(fvc::grad(U))) 
            -2*muEff*Km*fvc::grad(fvc::grad(lnr))
            +(2/3)*Km*fvc::laplacian(muEff,lnr)*I
		); 
        //tauMC = tau_v1;
       //Modified NS-stress tensor tau --(transformation)--> tau_v
        volTensorField tau_v("tau_v", muEff*fvc::grad(U)+ tau_v1 );
        // stress Tensor part of tau_RNS without visosity terms
       // tau_RNS = tau_V + tau_v2 ; 
        volTensorField tau_v2
        (   "tau_v2",
             tau_v1 
            -Km*Km*(fvc::grad(rho)*fvc::grad(rho))/rho 
            +Km*(U*fvc::grad(rho)) 
            +Km*(fvc::grad(rho)*U) 
        );  

        volTensorField tau_RNS
        (   
            IOobject
            (
                "tau_RNS",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            muEff*fvc::grad(U) + tau_v2
        );
     
        // --- Solve density RNS
        solve(
		fvm::ddt(rho) + fvc::div(phi) - fvc::laplacian(Km,rho)
	     );

        // --- Solve momentum RNS
        //solve(fvm::ddt(rhoU) + fvc::div(phiUp));
        solve(
                fvm::ddt(rhoU) + fvc::div(phiUp)    //phiUp = UU+pI
            );

        U.ref() =  rhoU() /rho();
        U.correctBoundaryConditions();
        rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();


 
        if (!inviscid)
        {
            solve
            (
                fvm::ddt(rho,U) - fvc::ddt(rho, U)
              - fvm::laplacian(muEff, U)
              - fvc::div(tau_v2)
              - Km*Km*fvc::grad(fvc::laplacian(rho))
              + Km*fvc::grad(fvc::div(rhoU))
            ); // 
            rhoU = rho*U;
        }

        // --- Solve energy
        surfaceScalarField sigmaDotU
        (
            "sigmaDotU",
            (
                fvc::interpolate(muEff)*mesh.magSf()*fvc::snGrad(U)
              + fvc::dotInterpolate(mesh.Sf(), -tau_v2) //tauMCS
            )
          & (a_pos*U_pos + a_neg*U_neg)
        );
	   // Below -> N1, N2, N3, N4, N5 and N6 based on RNS paper energy equation : Appendix
        volVectorField N1("N1", -(U & fvc::grad(rho))*U - 0.5*(U&U)*fvc::grad(rho));
        volVectorField N2("N2",(U&fvc::grad(rho))*(fvc::grad(lnr))+0.5*magSqr(fvc::grad(rho))*U/rho);
        volVectorField N3("N3",-0.5*magSqr(fvc::grad(rho))*((fvc::grad(lnr))/rho)) ;
        volScalarField N41("N41",fvc::div(rho*(U*U) + rho*I/psi - tau_RNS )&fvc::grad(lnr));
        volScalarField N42("N42",-U&(fvc::grad(lnr)*fvc::div(rhoU) - fvc::grad(fvc::div(rhoU))));
        volScalarField N4("N4",N41+N42);
        volScalarField N5("N5",
                            fvc::laplacian(rho)*(U&fvc::grad(lnr))
                            - (U&fvc::grad(fvc::laplacian(rho))) 
                            +  0.5 *magSqr(fvc::grad(rho))*((fvc::div(rhoU))/rho)/rho       
                        );        
        volScalarField N6
        (
            "N6",
            -0.5*magSqr(fvc::grad(rho))*((fvc::laplacian(rho))/rho)/rho 
		);

 
        solve
        (
            fvm::ddt(rhoE)                                                    //e0
          + fvc::div(phiEp)        // div[( rho E + 0.5 rho u^2 )*U + p*U ]   //e1
          - fvc::div(sigmaDotU)
          + fvc::div(Km*tau_v & fvc::grad(lnr))     // div[ tau_v & Uv ]                      // e2
          // - fvc::div((Km/rho)*tau_v & fvc::grad(rho))
          // + fvc::div(Km*N1 + Km*Km*N2 + Km*Km*Km*N3)
          // - fvc::div(Km*(rho*e*fvc::grad(rho)/rho - ((rho*rPsi*I)&(fvc::grad(rho)/rho))))
          // + Km*N4
          // + Km*Km*N5
          // + Km*Km*Km*N6
        );
 
        e = rhoE/rho - 0.5*magSqr(U);
        e.correctBoundaryConditions();
        thermo.correct();
        rhoE.boundaryFieldRef() ==
            rho.boundaryField()*
            (
                e.boundaryField() + 0.5*magSqr(U.boundaryField())
                   
            );

        if (!inviscid)
        {
            solve
            (
                  fvm::ddt(rho, e) - fvc::ddt(rho, e)
                - fvm::laplacian(turbulence->alphaEff(),e)               //qRNS1  eq 28
                - fvc::div(Km*rho*e*fvc::grad(lnr))        //qRNS2
                - fvc::div((Km*rho*I/psi)&fvc::grad(lnr))  //qRNS3
                + fvc::div(Km*N1 + Km*Km*N2 + Km*Km*Km*N3)  // e4
                + Km*N4
                + Km*Km*N5
                + Km*Km*Km*N6 //e6
            );
            thermo.correct(); 
            rhoE = rho*(e + 0.5*magSqr(U));
        }
       
        p.ref() =
            rho()
           /psi();
        p.correctBoundaryConditions();
        rho.boundaryFieldRef() == psi.boundaryField()*p.boundaryField();

        turbulence->correct();

	    runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //