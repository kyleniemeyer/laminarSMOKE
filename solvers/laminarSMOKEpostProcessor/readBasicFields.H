Info<< "Reading field T..." << endl;
volScalarField T
(
    IOobject
    (
	"T",
	runTime.timeName(),
	mesh,
	IOobject::MUST_READ,
	IOobject::NO_WRITE
    ),
    mesh
);

Info<< "Reading field p..." << endl;
volScalarField p
(
    IOobject
    (
	"p",
	runTime.timeName(),
	mesh,
	IOobject::MUST_READ,
	IOobject::NO_WRITE
    ),
    mesh
);
