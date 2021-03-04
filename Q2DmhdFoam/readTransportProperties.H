    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar rho0
    (
        transportProperties.lookup("rho0")
    );

    dimensionedScalar T0
    (
        transportProperties.lookup("T0")
    );

    dimensionedScalar nu
    (
        transportProperties.lookup("nu")
    );

    dimensionedScalar Cp
    (
        transportProperties.lookup("Cp")
    );

    dimensionedScalar DT
    (
        dimensionedScalar(transportProperties.lookup("k"))/
        (rho0*Cp)
    );

    dimensionedScalar beta
    (
        transportProperties.lookup("beta")
    );

    dimensionedScalar a
    (
        transportProperties.lookup("a")
    );

    dimensionedScalar b
    (
        transportProperties.lookup("b")
    );
    
    dimensionedScalar sigma
    (
        transportProperties.lookup("sigma")
    );

    dimensionedScalar tauHa
    (
        b*sqrt(rho0/(sigma*nu))
    );

    dimensionedScalar q0
    (
        transportProperties.lookup("q0")
    );

    dimensionedScalar m
    (
        transportProperties.lookup("m")
    );

    scalar coord
    (
        readScalar(transportProperties.lookup("coord"))
    );

    //  Read centerline velocity for channel simulations
    dimensionedVector Ubar(transportProperties.lookup("Ubar"));

    //  Magnitude of Ubar
    dimensionedScalar magUbar = mag(Ubar);
