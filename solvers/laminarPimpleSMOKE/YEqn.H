// Mass Fluxes
#include "fluxes.H"

// Discretization schemes
tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi)")
    )
);


{
    double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
    
    volScalarField Yt = 0.0*Y[0];

    for (label i=0; i<Y.size(); i++)
    {
        if (i != inertIndex)
        {
           	volScalarField& Yi = Y[i];
	    	volScalarField& Dmixi = Dmix[i];

		// Basic equations
		fvScalarMatrix YiEqn
		(
		     fvm::ddt(rho, Yi)
		   + mvConvection->fvmDiv(phi, Yi)
		   - fvm::laplacian(rho*Dmixi, Yi)
		     ==
                   - fvm::div(Jc,Yi, "div(Jc,Yi)") + 						
		     fvOptions(rho, Yi)
		);

		// Add Soret effect
		if (soretEffect == true)
		{ 
			if (soretEffectList[i] == true)
				YiEqn -= fvc::laplacian(rho*Dsoret[indexSoret[i]]/T, T, "laplacian(teta,Yi)");
		}

		// Add thermophoretic effect
		if (thermophoreticEffect == true)
		{
			if (thermophoreticEffectList[i] == true)
				YiEqn -= fvc::laplacian(0.55*mu/T*Yi, T, "laplacian(teta,Yi)");
		}

		// Solve
		YiEqn.relax();
		fvOptions.constrain(YiEqn);
		YiEqn.solve(mesh.solver("Yi"));
		fvOptions.correct(Yi);	
			
		// Sum of mass fractions
            	Yi.max(0.0);
           	Yt += Yi;
        }
    }

    Info << "Inert species is " << Y[inertIndex].name() << " with local index equal to " << inertIndex << endl;
    Y[inertIndex] = scalar(1.0) - Yt;
    Y[inertIndex].max(0.0);

    double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
	
    Info << "Transport equations of species solved in " << tEnd - tStart << " s " << endl;
}

