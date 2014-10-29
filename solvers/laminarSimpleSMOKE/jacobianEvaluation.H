{
	Info<< "Jacobian evaluation..." << endl;

	double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();	
	{
		const unsigned int NC = thermodynamicsMapXML->NumberOfSpecies();

		const scalarField& TCells = T.internalField();
		const scalarField& pCells = p.internalField(); 
	
		Eigen::MatrixXd J(NC+1,NC+1);
		OpenSMOKE::OpenSMOKEVectorDouble Source(thermodynamicsMapXML->NumberOfSpecies()+1);	
		OpenSMOKE::OpenSMOKEVectorDouble y(thermodynamicsMapXML->NumberOfSpecies()+1);

		forAll(TCells, celli)
		{
			for(int i=0;i<NC;i++)
				y[i+1] = Y[i].internalField()[celli];
			y[NC+1] = TCells[celli];

			reactionSourceTerms(	*thermodynamicsMapXML, *kineticsMapXML, y, pCells[celli], Source);
			reactionJacobian( *thermodynamicsMapXML, *kineticsMapXML, y, pCells[celli], J );

			for(int i=0;i<NC+1;i++)
				sourceImplicit[i].internalField()[celli] = J(i,i);
			
			for(int i=0;i<NC+1;i++)
				sourceExplicit[i].internalField()[celli] = Source[i+1] - J(i,i)*y[i+1];
		}
	}
	double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

	Info << "   Jacobian evaluations in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per cell)" << endl;
}
