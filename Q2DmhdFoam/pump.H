// Relaxation factor ready if needed
//scalar relaxGradP = 1;
// Get X component (non-default)
// ONLY VALID FOR X DIRECTION OF THE FLOW
volScalarField Ux = U.component(vector::X);

// Extract the velocity in the flow direction
dimensionedScalar magUbarStar = Ux.weightedAverage(mesh.V());

// Calculate the pressure gradient increment needed to
// adjust the average flow-rate to the correct value
dimensionedScalar gradPplus =
  (magUbar - magUbarStar)/rUA.weightedAverage(mesh.V());

// Get new values
Ux    += rUA*gradPplus;
gradP += gradPplus;
gradP.value()  = max(0,gradP.value());

// Need to add Ux back into  U
U = volVectorField
    (
        vector(1, 0, 0)*Ux +
        vector(0, 1, 0)*Ux*0 +
        vector(0, 0, 1)*Ux*0
    );
