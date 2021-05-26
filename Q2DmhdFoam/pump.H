// Relaxation factor ready if needed
//scalar relaxGradP = 1;
// Get X component (non-default)
// ONLY VALID FOR X DIRECTION OF THE FLOW
//volScalarField Ux = U.component(vector::X);

// Extract the velocity in the flow direction
magUbarStar = (flowDirection & U)().weightedAverage(mesh.V());

// Calculate the pressure gradient increment needed to
// adjust the average flow-rate to the correct value
gradPplus =
  (magUbar - magUbarStar)/rUA.weightedAverage(mesh.V());

// Get new values
U    += gradPplus*rUA*flowDirection;
gradP += gradPplus;

Info<< "Uncorrected Ubar = " << magUbarStar.value() << tab
    << "pressure gradient = " << gradP.value() << endl;
