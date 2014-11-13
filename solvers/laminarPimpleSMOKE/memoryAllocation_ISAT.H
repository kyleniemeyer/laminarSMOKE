
//- HOM
Eigen::VectorXd phiISAT_HOM;
Eigen::VectorXd RphiISAT_HOM;
Eigen::MatrixXd mapGrad_HOM;

if(isatCheck == true)
{
	phiISAT_HOM.resize(thermodynamicsMapXML->NumberOfSpecies()+1);
	RphiISAT_HOM.resize(thermodynamicsMapXML->NumberOfSpecies()+1);
	mapGrad_HOM.resize(thermodynamicsMapXML->NumberOfSpecies()+1,thermodynamicsMapXML->NumberOfSpecies()+1);
}
