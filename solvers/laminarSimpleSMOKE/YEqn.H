
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

if (thermophoreticEffect == true)
{
	Info << "Calculating thermophoretic effect... " << endl;
	rhoVT = -0.55*mu*fvc::grad(T)/T;
	VT = rhoVT/rho;
}
word thermophoreticDivScheme("div(rhoVT,Yi)");
word soretDivScheme("soretEffect");

{
    double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
    
    volScalarField Yt = 0.0*Y[0];

    for (label i=0; i<Y.size(); i++)
    {
        if (i != inertIndex)
        {
           	volScalarField& Yi = Y[i];
	    	volScalarField& Dmixi = Dmix[i];

	    	if (diffusivityCorrection != DIFFUSIVITY_CORRECTION_NONE)
			{
				dimensionedScalar Mi( "Mi", dimensionSet(1,0,0,0,-1,0,0), thermodynamicsMapXML->MW()[i+1] ); 

				sumDiffusionCorrections *= 0.;
				for (label j=0; j<Y.size(); j++)
				{
					volScalarField& Yj = Y[j];	
					dimensionedScalar Mj( "Mj", dimensionSet(1,0,0,0,-1,0,0), thermodynamicsMapXML->MW()[j+1] );
					sumDiffusionCorrections += fvc::laplacian(rho*MWmix/Mj*Dmixi*Yi, Yj);
				
					if ( (i==j) && (diffusivityCorrection == DIFFUSIVITY_CORRECTION_SEMI_IMPLICIT) )
						sumDiffusionCorrections -= fvc::laplacian(rho*MWmix/Mj*Dmixi*Yi, Yj);
				}

				if (diffusivityCorrection == DIFFUSIVITY_CORRECTION_SEMI_IMPLICIT)
				{
					if (thermophoreticEffectList[i] == true)
					{
						fvScalarMatrix YiEqn
						(
						    mvConvection->fvmDiv(phi, Yi)
						  - fvm::laplacian(rho*Dmixi, Yi)
						  + fvm::laplacian(rho*Dmixi*MWmix*Yi/Mi, Yi)
						  + sumDiffusionCorrections
						  + fvc::div(rhoVT*Yi,thermophoreticDivScheme)
						 ==
							sourceExplicit[i] + fvm::Sp(sourceImplicit[i],Yi) +	
							fvOptions(rho, Yi)
						);

						YiEqn.relax();
						fvOptions.constrain(YiEqn);
						YiEqn.solve(mesh.solver("Yi"));
						fvOptions.correct(Yi);
					}
					else 
					{
						if (soretEffectList[i] == true)
						{
							fvScalarMatrix YiEqn
							(
							    mvConvection->fvmDiv(phi, Yi)
							  - fvm::laplacian(rho*Dmixi, Yi)
							  + fvm::laplacian(rho*Dmixi*MWmix*Yi/Mi, Yi)
							  + sumDiffusionCorrections
							  + fvc::div(rho*Dsoret[indexSoret[i]]/T*fvc::grad(T),soretDivScheme)
							 ==
								sourceExplicit[i] + fvm::Sp(sourceImplicit[i],Yi) +	
								fvOptions(rho, Yi)
							);

							YiEqn.relax();
							fvOptions.constrain(YiEqn);
							YiEqn.solve(mesh.solver("Yi"));
							fvOptions.correct(Yi);
						}
						else
						{
							fvScalarMatrix YiEqn
							(
							    mvConvection->fvmDiv(phi, Yi)
							  - fvm::laplacian(rho*Dmixi, Yi)
							  + fvm::laplacian(rho*Dmixi*MWmix*Yi/Mi, Yi)
							  + sumDiffusionCorrections
							 ==
								sourceExplicit[i] + fvm::Sp(sourceImplicit[i],Yi) +	
								fvOptions(rho, Yi)
							);

							YiEqn.relax();
							fvOptions.constrain(YiEqn);
							YiEqn.solve(mesh.solver("Yi"));
							fvOptions.correct(Yi);
						}				
					}						
				}
				else if (diffusivityCorrection == DIFFUSIVITY_CORRECTION_EXPLICIT)
				{		
					if (thermophoreticEffectList[i] == true)
					{			
						fvScalarMatrix YiEqn
						(
						    mvConvection->fvmDiv(phi, Yi)
						  - fvm::laplacian(rho*Dmixi, Yi)
						  + sumDiffusionCorrections
						  + fvc::div(rhoVT*Yi,thermophoreticDivScheme)
						 ==
							sourceExplicit[i] + fvm::Sp(sourceImplicit[i],Yi) +	
							fvOptions(rho, Yi)
						);

						YiEqn.relax();
						fvOptions.constrain(YiEqn);
						YiEqn.solve(mesh.solver("Yi"));
						fvOptions.correct(Yi);			
					}
					else
					{			
						fvScalarMatrix YiEqn
						(
					            mvConvection->fvmDiv(phi, Yi)
						  - fvm::laplacian(rho*Dmixi, Yi)
						  + sumDiffusionCorrections
						 ==
							sourceExplicit[i] + fvm::Sp(sourceImplicit[i],Yi) +	
							fvOptions(rho, Yi)
						);

						YiEqn.relax();
						fvOptions.constrain(YiEqn);
						YiEqn.solve(mesh.solver("Yi"));
						fvOptions.correct(Yi);			
					}
				}
			}
			else
			{
				if (thermophoreticEffectList[i] == true)
				{
					fvScalarMatrix YiEqn
					(
					    mvConvection->fvmDiv(phi, Yi)
					  - fvm::laplacian(rho*Dmixi, Yi)
					  + fvc::div(rhoVT*Yi,thermophoreticDivScheme)
					 ==
						sourceExplicit[i] + fvm::Sp(sourceImplicit[i],Yi) +	
						fvOptions(rho, Yi)
					);

					YiEqn.relax();
					fvOptions.constrain(YiEqn);
					YiEqn.solve(mesh.solver("Yi"));
					fvOptions.correct(Yi);	
				}
				else
				{
					fvScalarMatrix YiEqn
					(
					    mvConvection->fvmDiv(phi, Yi)
					  - fvm::laplacian(rho*Dmixi, Yi)
					 ==
					     sourceExplicit[i] + fvm::Sp(sourceImplicit[i],Yi) +						
					     fvOptions(rho, Yi)
					);

					YiEqn.relax();
					fvOptions.constrain(YiEqn);
					YiEqn.solve(mesh.solver("Yi"));
					fvOptions.correct(Yi);	
				}
			}

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
