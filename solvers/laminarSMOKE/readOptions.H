Info<< "Reading Solver options\n" << endl;
IOdictionary solverOptions
(
	IOobject
	(
		"solverOptions",
		U.time().constant(),
		U.db(),
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	)
);

label inertIndex;
OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>*			thermodynamicsMapXML; 
OpenSMOKE::KineticsMap_CHEMKIN<double>* 				kineticsMapXML;
OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>* 		transportMapXML;

const dictionary& kineticsDictionary = solverOptions.subDict("Kinetics");
{
	//- Kinetics folder
	Foam::string kinetics_folder = kineticsDictionary.lookup("folder");

	// Reading homogeneous data
	Info<< " * loading homogeneous kinetics...\n" << endl;
	{	
		word kinetics("opensmoke");
		boost::filesystem::path path_kinetics = kinetics_folder;

		rapidxml::xml_document<> doc;
		std::vector<char> xml_string;
		OpenSMOKE::OpenInputFileXML(doc, xml_string,path_kinetics / "kinetics.xml");

		double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
		
		thermodynamicsMapXML = new OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>(doc); 
		transportMapXML = new OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>(doc); 
		kineticsMapXML = new OpenSMOKE::KineticsMap_CHEMKIN<double>(*thermodynamicsMapXML, doc); 
							
		double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
		std::cout << " * Time to read XML file: " << tEnd-tStart << std::endl;
	}

	//- Inert species
	word inertSpecies(kineticsDictionary.lookup("inertSpecies"));
	inertIndex = thermodynamicsMapXML->IndexOfSpecies(inertSpecies)-1;
}

//- Mass fractions tolerance
scalar massFractionsTol = 1e-10;

//- Internal models
enum { STRANG_MOMENTUM_TRANSPORT_REACTION, STRANG_MOMENTUM_REACTION_TRANSPORT} strangAlgorithm;
enum { DIFFUSIVITY_MODEL_MULTICOMPONENT, DIFFUSIVITY_MODEL_LEWIS_NUMBERS} diffusivityModel;

// Physical model
Switch energyEquation;
Switch homogeneousReactions;
Switch iMassDiffusionInEnergyEquation;
Switch constPressureBatchReactor;
Switch includeDpDt;
Switch soretEffect;
Eigen::VectorXd LewisNumbers;
std::vector<bool> soretEffectList;
std::vector<int> indexSoret;
scalar thermodynamicPressure = 101325.;

const dictionary& physicalModelDictionary = solverOptions.subDict("PhysicalModel");
{
	energyEquation = Switch(physicalModelDictionary.lookup(word("energyEquation")));
	homogeneousReactions = Switch(physicalModelDictionary.lookup(word("homogeneousReactions")));
	includeDpDt = Switch(physicalModelDictionary.lookup(word("includeDpDt")));
	iMassDiffusionInEnergyEquation = Switch(physicalModelDictionary.lookup(word("massDiffusionInEnergyEquation")));
	
	// Soret effect
	soretEffect = Switch(physicalModelDictionary.lookup(word("soretEffect")));
	soretEffectList.resize(thermodynamicsMapXML->NumberOfSpecies());
	indexSoret.resize(thermodynamicsMapXML->NumberOfSpecies());

	for(unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
		soretEffectList[i] = false;	
	for(unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
		indexSoret[i] = -1;
	
	// Algorithm
	{
		word strang(physicalModelDictionary.lookup("strangAlgorithm"));
		if (strang == "MomentumTransportReaction")	strangAlgorithm = STRANG_MOMENTUM_TRANSPORT_REACTION;
		else if (strang == "MomentumReactionTransport")	strangAlgorithm = STRANG_MOMENTUM_REACTION_TRANSPORT;
		else
		{
			Info << "Wrong strangAlgorithm option: MomentumTransportReaction || MomentumReactionTransport" << endl;
			abort();
		}
	}
	
	// Diffusivity of species
	{
		word diffusivity(physicalModelDictionary.lookup("diffusivityModel"));
		if (diffusivity == "multi-component")		diffusivityModel = DIFFUSIVITY_MODEL_MULTICOMPONENT;
		else if (diffusivity == "lewis-numbers")	diffusivityModel = DIFFUSIVITY_MODEL_LEWIS_NUMBERS;
		else
		{
			Info << "Wrong diffusivityModel option: multi-component || lewis-numbers" << endl;
			abort();
		}
		
		if (diffusivityModel == DIFFUSIVITY_MODEL_LEWIS_NUMBERS)
		{
			LewisNumbers.resize(thermodynamicsMapXML->NumberOfSpecies());
			
			const dictionary& lewisNumberDictionary = solverOptions.subDict("LewisNumbers");
			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				LewisNumbers(i) = readScalar(lewisNumberDictionary.lookup(thermodynamicsMapXML->NamesOfSpecies()[i]));
			
			Info << endl;
			Info << "Diffusivities of species will be calculated using the following Lewis numbers: " << endl;
			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				Info << " " << thermodynamicsMapXML->NamesOfSpecies()[i] << "\t" << LewisNumbers(i) << endl;
			Info << endl;
		}
	}
	
	//- Constant pressure reactors
	constPressureBatchReactor = Switch(physicalModelDictionary.lookup(word("constPressureBatchReactor")));

	// Themodynamic pressure [Pa]
	thermodynamicPressure = readScalar(physicalModelDictionary.lookup("thermodynamicPressure"));
}

//- Detect spark
const dictionary& sparkDictionary = solverOptions.subDict("Spark");
sparkModel spark;
spark.Read(sparkDictionary);

Eigen::VectorXd outputFormationRatesIndices;
Eigen::VectorXd outputDiffusivitiesIndices;
const dictionary& outputDictionary = solverOptions.subDict("Output"); 
{
	{
		Switch outputFormationRates = Switch(outputDictionary.lookup(word("formationRates")));
		if (outputFormationRates == true)
		{
			List<word>  listFormationRates(outputDictionary.lookup("listFormationRates"));
			outputFormationRatesIndices.resize(listFormationRates.size());
			for (int i=0;i<listFormationRates.size();i++)
				outputFormationRatesIndices(i) = thermodynamicsMapXML->IndexOfSpecies(listFormationRates[i])-1;
		}
	}

	{
		Switch outputDiffusivities = Switch(outputDictionary.lookup(word("diffusivities")));
		if (outputDiffusivities == true)
		{
			List<word>  listDiffusivities(outputDictionary.lookup("listDiffusivities"));
			outputDiffusivitiesIndices.resize(listDiffusivities.size());
			for (int i=0;i<listDiffusivities.size();i++)
				outputDiffusivitiesIndices(i) = thermodynamicsMapXML->IndexOfSpecies(listDiffusivities[i])-1;
		}
	}
}

// Radiative heat transfer
radiativeHeatTransferClass radiativeHeatTransfer;

const dictionary& radiativeHeatTransferDictionary = solverOptions.subDict("Radiation");
{
	word radiation(radiativeHeatTransferDictionary.lookup("radiationModel"));
	     if (radiation == "none") 		radiativeHeatTransfer.SetRadiationModel(RADIATION_NONE);
	else if (radiation == "optically-thin") radiativeHeatTransfer.SetRadiationModel(RADIATION_OPTICALLY_THIN);
	else
	{
		Info << "Wrong radiationModel option: none || optically-thin" << endl;
		abort();
	}

	if (radiativeHeatTransfer.radiationModel() == RADIATION_OPTICALLY_THIN)
	{
		radiativeHeatTransfer.SetThermodynamics(*thermodynamicsMapXML);
		radiativeHeatTransfer.SetEnvironmentTemperature( readScalar(radiativeHeatTransferDictionary.lookup("environmentTemperature")) );
	}
}

// Thermophoretic effect
Switch thermophoreticEffect;
std::vector<bool> thermophoreticEffectList;
const dictionary& thermophoreticEffectDictionary = solverOptions.subDict("ThermophoreticEffect");
{ 
	thermophoreticEffect = Switch(thermophoreticEffectDictionary.lookup(word("thermophoreticEffect")));
	thermophoreticEffectList.resize(thermodynamicsMapXML->NumberOfSpecies());
	for(unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
		thermophoreticEffectList[i] = false;

	// Populate the list of thermophoretic species
	if (thermophoreticEffect == true)
	{
		word tag = thermophoreticEffectDictionary.lookup(word("thermophoreticEffectTag"));

		unsigned int count = 0;
		for(unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
		{
			std::string name = thermodynamicsMapXML->NamesOfSpecies()[i];
			std::string prefix = name.substr(0,tag.size());
			if (prefix == tag)
			{
				thermophoreticEffectList[i] = true;
				count++;
			}
		}
			
		if (count == 0)
			std::cout << "Warning: No " << tag << " species were found in the kinetic scheme! Are you sure you want to account for the thermophoretic effect?" << std::endl;
		else
			std::cout << "Thermophoretic effect: Number of " << tag << " species found equal to " << count << std::endl;
	}
} 

/*********************************************************************************************************************/ 

#if STEADYSTATE != 1

// Batch reactor homogeneous: ode parameters
const dictionary& odeHomogeneousDictionary = solverOptions.subDict("OdeHomogeneous");
OpenSMOKE::ODE_Parameters odeParameterBatchReactorHomogeneous;
{
	//- Mass fractions tolerance
	scalar relTolerance = readScalar(odeHomogeneousDictionary.lookup("relTolerance"));
	scalar absTolerance = readScalar(odeHomogeneousDictionary.lookup("absTolerance"));
	odeParameterBatchReactorHomogeneous.SetRelativeTolerance(relTolerance);
	odeParameterBatchReactorHomogeneous.SetAbsoluteTolerance(absTolerance);
	
	//- Full pivoting (only for OpenSMOKE solver)
	Switch fullPivoting(odeHomogeneousDictionary.lookup("fullPivoting"));
	odeParameterBatchReactorHomogeneous.SetFullPivoting(fullPivoting);
	
	//- Maximum order of integration (only for OpenSMOKE solver)
	label maximumOrder = readLabel(odeHomogeneousDictionary.lookup("maximumOrder"));
	odeParameterBatchReactorHomogeneous.SetMaximumOrder(maximumOrder);
	
	//- 
	
	// Type
	word homogeneousODESolverString(odeHomogeneousDictionary.lookup("odeSolver"));
	if (	homogeneousODESolverString != "OpenSMOKE" 	&& homogeneousODESolverString != "DVODE"  && 
			homogeneousODESolverString != "DLSODE" 		&& homogeneousODESolverString != "DLSODA" && 
			homogeneousODESolverString != "CVODE" 		&& homogeneousODESolverString != "DASPK"  &&
			homogeneousODESolverString != "MEBDF" 		&& homogeneousODESolverString != "RADAU5"  
	   )
	{
		Info << "Wrong homogeneous ODE Solver: OpenSMOKE || DVODE || DLSODE || DLSODA || CVODE || DASPK || MEBDF || RADAU5" << endl;
		abort();
	}

	if (homogeneousODESolverString == "OpenSMOKE")	odeParameterBatchReactorHomogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE);
	if (homogeneousODESolverString == "DVODE") 		odeParameterBatchReactorHomogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DVODE);
	if (homogeneousODESolverString == "DLSODE") 	odeParameterBatchReactorHomogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE);
	if (homogeneousODESolverString == "DLSODA") 	odeParameterBatchReactorHomogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODA);
	if (homogeneousODESolverString == "CVODE") 		odeParameterBatchReactorHomogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE);
	if (homogeneousODESolverString == "DASPK") 		odeParameterBatchReactorHomogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DASPK);
	if (homogeneousODESolverString == "MEBDF") 		odeParameterBatchReactorHomogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_MEBDF);
	if (homogeneousODESolverString == "RADAU5") 	odeParameterBatchReactorHomogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_RADAU5);	
	
}


// Check type
{
	if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DVODE)
	{
		#if OPENSMOKE_USE_DVODE != 1
		{
			Info << "The solver was compiled without the DVODE support. Please choose a different ODE solver." << endl;
			abort();
		}
		#endif
	}
	else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE)
	{
		#if OPENSMOKE_USE_SUNDIALS != 1
		{
			Info << "The solver was compiled without the CVODE support. Please choose a different ODE solver." << endl;
			abort();
		}
		#endif
	}
	else if ( odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE ||
	          odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODA )
	{
		#if OPENSMOKE_USE_ODEPACK != 1
		{
			Info << "The solver was compiled without the ODEPACK support. Please choose a different ODE solver." << endl;
			abort();
		}
		#endif
	}
	else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_RADAU5)
	{
		#if OPENSMOKE_USE_RADAU != 1
		{
			Info << "The solver was compiled without the RADAU support. Please choose a different ODE solver." << endl;
			abort();
		}
		#endif
	}	
	else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DASPK)
	{
		#if OPENSMOKE_USE_DASPK != 1
		{
			Info << "The solver was compiled without the DASPK support. Please choose a different ODE solver." << endl;
			abort();
		}
		#endif
	}
	else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_MEBDF)
	{
		#if OPENSMOKE_USE_MEBDF != 1
		{
			Info << "The solver was compiled without the MEBDF support. Please choose a different ODE solver." << endl;
			abort();
		}
		#endif
	}			
}
#endif

#if STEADYSTATE != 1
#if OPENSMOKE_USE_ISAT == 1
    #include "readOptions_ISAT.H"
#endif
#endif

Info<< "Solver options correctly read\n" << endl;
