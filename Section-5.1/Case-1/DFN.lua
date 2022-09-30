Simulation
 {
  name  = 'DFN',
  options = {
  femAutoResultAttributePrefix = true,
  disabledWarnings = {'FixedBcForInvalidDof'},
  },     
 }
dofile('$SIMULATIONDIR/$SIMULATIONNAME_model.lua')
dofile('$SIMULATIONDIR/$SIMULATIONNAME_solution.lua')
