//- Reading ISAT parameters
Switch isatCheck(solverOptions.subDict("ISAT").lookup("ISAT_IntegrationMethod"));

// mandatory if ISAT ON
scalar epsTol_HOM = solverOptions.subDict("ISAT").lookupOrDefault<double>("epsTolISAT_HOM", epsTol_HOM);
scalar maxSizeMRU = solverOptions.subDict("ISAT").lookupOrDefault<int>("maxSizeMRU", 100);
scalar maxSizeMFU = solverOptions.subDict("ISAT").lookupOrDefault<int>("maxSizeMFU", 100);
scalar maxSizeBT = solverOptions.subDict("ISAT").lookupOrDefault<int>("maxSizeBT", 50000);
scalar maxSearchMRU = solverOptions.subDict("ISAT").lookupOrDefault<int>("maxSearchMRU", 10);
scalar maxSearchMFU = solverOptions.subDict("ISAT").lookupOrDefault<int>("maxSearchMFU", 30);
Switch searchMRU(solverOptions.subDict("ISAT").lookupOrDefault<Switch>("searchMRU", true)); 
Switch searchMFU(solverOptions.subDict("ISAT").lookupOrDefault<Switch>("searchMFU", true));
Switch searchBF(solverOptions.subDict("ISAT").lookupOrDefault<Switch>("searchBruteForce", false));
Switch clearIfFull(solverOptions.subDict("ISAT").lookupOrDefault<Switch>("clearingIfFull", false));
scalar maxTimeOldCoeff = solverOptions.subDict("ISAT").lookupOrDefault<double>("maxTimeOldCoeff", 0.7);
scalar maxGrowCoeff = solverOptions.subDict("ISAT").lookupOrDefault<double>("maxGrowCoeff", 0.5);
scalar maxHeightCoeff = solverOptions.subDict("ISAT").lookupOrDefault<double>("maxHeightCoeff", 20.);
scalar maxUsedCoeff = solverOptions.subDict("ISAT").lookupOrDefault<double>("maxUsedCoeff", 20.);
Switch printStats(solverOptions.subDict("ISAT").lookupOrDefault("stats", false));
word typeMapGradSolver = solverOptions.subDict("ISAT").lookupOrDefault<word>("typeMapGradSolver","PartialPiv");

int typeSolver = 0;
if (typeMapGradSolver != "PartialPiv" && typeMapGradSolver != "FullPiv")
{
	Info << "Wrong typeMapGradSolver options: PartialPiv || FullPiv" << endl;
	abort();
} 
else 
{
	if(typeMapGradSolver == "FullPiv") 
		typeSolver = 0;
	else
		typeSolver = 1;		
}


//- ISAT scale factor 
Eigen::VectorXd scaleFactor_HOM;
scaleFactor_HOM.setConstant(thermodynamicsMapXML->NumberOfSpecies()+1,1.);
if(isatCheck == true) 
{
	const unsigned int NC  = thermodynamicsMapXML->NumberOfSpecies();

	massFractionsTol = epsTol_HOM;
	
	dictionary scaleFact(solverOptions.subDict("ISAT").subDict("ScaleFactor"));
	for(unsigned int i=0;i<NC;i++) 
	{	
		if(scaleFact.found(thermodynamicsMapXML->NamesOfSpecies()[i])) 
			scaleFactor_HOM(i) = 1./readScalar(scaleFact.lookup(thermodynamicsMapXML->NamesOfSpecies()[i]));
		else 
			scaleFactor_HOM(i) = 1./readScalar(scaleFact.lookup("otherSpecies"));
	}
	scaleFactor_HOM(NC) = 1./readScalar(scaleFact.lookup("T"));
}
//- ISAT scaling error factors
Eigen::VectorXd scaleError_HOM;
scaleError_HOM.setConstant(thermodynamicsMapXML->NumberOfSpecies()+1,1.);
if(isatCheck == true) 
{
	const unsigned int NC  = thermodynamicsMapXML->NumberOfSpecies();
	
	dictionary scaleErr(solverOptions.subDict("ISAT").subDict("ScalingError"));
	
	for(unsigned int i=0;i<NC;i++) 
	{
		if(scaleErr.found(thermodynamicsMapXML->NamesOfSpecies()[i])) 
			scaleError_HOM(i) = readScalar(scaleErr.lookup(thermodynamicsMapXML->NamesOfSpecies()[i]));
		else 
			scaleError_HOM(i) = readScalar(scaleErr.lookup("otherSpecies"));
	}
	scaleError_HOM(NC) = readScalar(scaleErr.lookup("T"));
}

//- ISAT HOM
ISAT *isat_HOM = new ISAT(scaleError_HOM, epsTol_HOM,thermodynamicsMapXML->NumberOfSpecies()+1);	

// - Setting ISAT_HOM param
isat_HOM->setMaxSizeBT(maxSizeBT);
isat_HOM->setMaxSizeMRU(maxSizeMRU);
isat_HOM->setMaxSizeMFU(maxSizeMFU);
isat_HOM->setMaxSearchMRU(maxSearchMRU);
isat_HOM->setMaxSearchMFU(maxSearchMFU);
isat_HOM->setFlagSearchMRU(searchMRU);
isat_HOM->setFlagSearchMFU(searchMFU);
isat_HOM->setFlagClearingIfFull(clearIfFull);
isat_HOM->setMaxGrowCoeff(maxGrowCoeff);
isat_HOM->setMaxHeightCoeff(maxHeightCoeff);
isat_HOM->setMaxTimeOldCoeff(maxTimeOldCoeff);
isat_HOM->setMaxUsedCoeff(maxUsedCoeff);

if(isatCheck == true) 
{
	const unsigned int NC  = thermodynamicsMapXML->NumberOfSpecies();
	
	Info << endl << "ISAT parameters " << endl;
	Info << "   homogeneous         : " << epsTol_HOM << endl;
	Info << "   typeMagGradSolver   : " << typeMapGradSolver << endl; 	
	
	Info << "   scaleFactors        : " << endl;
	for(unsigned int i=0;i<NC;i++)  
	{
		if(scaleFactor_HOM(i) != 1.) 
			Info << "      " << thermodynamicsMapXML->NamesOfSpecies()[i] << "\t       : " << scaleFactor_HOM(i) << endl;
	}
	Info << "       otherSpecies    : " << readScalar(solverOptions.subDict("ISAT").subDict("ScaleFactor").lookup("otherSpecies")) << endl;
	Info << "       T               : " << scaleFactor_HOM(NC) << endl;


	Info << "   clear if BT full    : " << clearIfFull << endl;

	Info << "   search in MRU       : " << searchMRU << endl;
	Info << "   search in MFU       : " << searchMFU << endl;
	Info << "   search BF           : " << searchBF << endl;

	Info << "   dimension parameters: " << endl;
	Info << "      max size BT      : " << maxSizeBT << endl;	
	if(searchMRU == true) 
	{
		Info << "      max size MRU     : " << maxSizeMRU << endl;
		Info << "      max search MRU   : " << maxSearchMRU << endl;
	}
	if(searchMFU == true) 
	{
		Info << "      max size MFU     : " << maxSizeMRU << endl;
		Info << "      max search MFU   : " << maxSearchMFU << endl;  
	}
	
	Info << "   balance parameters  : " << endl;
	Info << "      maxHeightCoeff   : " << maxHeightCoeff << endl;	
	Info << "      maxGrowCoeff     : " << maxGrowCoeff << endl;
	Info << "      maxUsedCoeff     : " << maxUsedCoeff << endl;
	Info << "      maxTimeOldCoeff  : " << maxTimeOldCoeff << endl; 

	Info << "   scaling error       : " << endl;
	for(unsigned int i=0;i<NC;i++)  
	{
		if(scaleError_HOM(i) != 1.) 
			Info << "      " << thermodynamicsMapXML->NamesOfSpecies()[i] << "\t       : " << scaleError_HOM(i) << endl;
	}
	Info << "       otherSpecies    : " << readScalar(solverOptions.subDict("ISAT").subDict("ScalingError").lookup("otherSpecies")) << endl;
	Info << "       T               : " << scaleError_HOM(NC) << endl;
	
	Info << endl;
	Info << "   stats               : " << printStats << endl << endl;
}

