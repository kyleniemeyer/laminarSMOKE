//- Reading velocity [m/s]
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

//- Reading pressure [Pa]
Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//- Allocating memory for p/rho [s2/m2]
Info<< "Creating field psi\n" << endl;
volScalarField psi
(
    IOobject
    (
        "psi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionSet(0, -2, 2, 0, 0)
);

#if STEADYSTATE == 1
label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);

dimensionedScalar rhoMax(simple.dict().lookup("rhoMax"));
dimensionedScalar rhoMin(simple.dict().lookup("rhoMin"));

//- Reading pressure [Pa]
volScalarField DpDt
(
    IOobject
    (
        "DpDt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("DpDt", dimensionSet(1, -1, -3, 0, 0), 0.0)
);

#endif


//- Reading temperature [K]
Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//- Allocating memory for mean molecular weight [kg/mol]
Info<< "Creating field MWmix\n" << endl;
volScalarField MWmix
(
    IOobject
    (
        "gas::MW",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(1, 0, 0, 0, -1)
);

//- Allocating memory for total concentration [kmol/m3]
Info<< "Creating field cTot\n" << endl;
volScalarField cTot
(
    IOobject
    (
        "gas::cTot",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(0, -3, 0, 0, 1)
);

//- Allocating memory for dynamic viscosity [kg/m/s]
Info<< "Creating field mu\n" << endl;
volScalarField mu
(
    IOobject
    (
        "gas::mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(1, -1, -1, 0, 0)
);

//- Allocating memory for thermal conductivity [W/m/K]
Info<< "Creating field lambda\n" << endl;
volScalarField lambda
(
    IOobject
    (
        "gas::lambda",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(1, 1, -3, -1, 0)
);

//- Allocating memory for radiative heat transfer term [W/m3]
volScalarField Qrad
(
    IOobject
    (
        "rad::Qrad",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("rad::Qrad", dimensionSet(1, -1, -3, 0, 0), 0.0)
);

//- Allocating memory for Planck mean absorption coefficient [1/m]
volScalarField aPlanck
(
    IOobject
    (
        "rad::aPlanck",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("rad::aPlanck", dimensionSet(0, 0, -1, 0, 0), 0.0)
);

//- Allocating memory for reaction heat [W/m3]
volScalarField Q
(
    IOobject
    (
        "gas::Q",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("gas::Q", dimensionSet(1, -1, -3, 0, 0), 0.0)
); 

//- Allocating memory for mean specific heat [J/kg/K]
Info<< "Creating field Cp\n" << endl;
volScalarField cp
(
    IOobject
    (
        "gas::Cp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(0, 2, -2, -1, 0)
);

//- Allocating memory for mean specific heat [J/kg/K]
Info<< "Creating field Cv\n" << endl;
volScalarField cv
(
    IOobject
    (
        "gas::Cv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionSet(0, 2, -2, -1, 0)
);

//- Allocating memory for total enthalpy [m2/s2]
Info<< "Creating field enthalpy\n" << endl;
volScalarField h
(
    IOobject
    (
        "gas::h",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(0, 2, -2, 0, 0)
);

Info<< "Creating massDiffusionInEnergyEquation\n" << endl;
volScalarField massDiffusionInEnergyEquation
(
    IOobject
    (
		"gas::massDiffusionInEnergyEquation",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
    ),
    mesh,
    dimensionSet(1, -1, -3, 0, 0)
);
massDiffusionInEnergyEquation = dimensionedScalar("zero", dimensionSet(1, -1, -3, 0, 0), 0.);
