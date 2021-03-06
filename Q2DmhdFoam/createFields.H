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

#   include "createPhi.H"
    
    Info<< "Reading field B\n" << endl;
    volScalarField B
    (
        IOobject
        (
            "B",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, piso.dict(), pRefCell, pRefValue);
    mesh.schemesDict().setFluxRequired(p.name());

    Info<< "Reading/calculating field rho\n" << endl;
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        rho0*(scalar(1) - beta*(T - T0))
    );
    

    Info<< "Calculating field sourceT\n" << endl;
    // Infer the sourceT BCs from the pressure (set all to zeroGradient)
    wordList sourceTBCTypes
    (
        p.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );
    volScalarField y = mesh.C().component(vector::Y);
    if (coord == 1)
    {
        y = mesh.C().component(vector::X);
    }
    if (coord == 2)
    {
        y = mesh.C().component(vector::Y);
    }
    if (coord == 3)
    {
        y = mesh.C().component(vector::Z);
    }
    volScalarField sourceT
    (
        IOobject
        (
            "sourceT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        q0*exp(- (y/a + 1) * m),
        sourceTBCTypes
    );


    Info<< "Creating field vorticity\n" << endl;
    volVectorField vorticity
    (
        IOobject
        (
            "vorticity",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("vorticity", dimVelocity/dimLength, vector(0.0, 0.0, 0.0))
    );

    
    Info<< "Creating gradP field for fully developed case\n" << endl;
    dimensionedScalar gradP
    (
        "gradP",
        dimensionSet(0, 1, -2, 0, 0), 
        0.0
    );

    // Temperature averaging variables
    dimensionedScalar tempVolAverage("tempVolAverage", dimensionSet(0,0,0,1,0,0,0), 0);
    dimensionedScalar magUbarStar("magUbarStar", dimVelocity, 0);
    dimensionedScalar gradPplus("gradPplus", dimensionSet(0, 1, -2, 0, 0), 0);
    
    // Other
    volScalarField rUA
    (
        IOobject
        (
            "rUA",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("rUA", dimTime, 0.0)
    );
