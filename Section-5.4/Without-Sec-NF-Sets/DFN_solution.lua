NumericalSolver {
  id          = 'solver',
  typeName    = 'ArmadilloSolver',
  description = 'Simple matrix solver',
  sparse = true,
}
------------------------------------- 
-- Physical Methods 
------------------------------------- 
PhysicalMethod { 
  id       = 'Geostatic',
  typeName = 'HydroFemPhysics.Flow', 
  type     = 'fem',

  mesh = 'mesh',
  elementGroups = {'gcep_rock',}, 
  materials = {'saturated'},
  isoParametric = false, 
  ruleSet = 1,

  properties = { material = 'materialH' },
} 

PhysicalMethod { 
  id       = 'Rock',
  typeName = 'HydroFemPhysics.Flow',
  type     = 'fem',

  mesh = 'mesh',
  elementGroups = {'gcep_rock',}, 
  materials = {'saturated'}, 
  isoParametric = false, 
  boundaryConditions={'bpm',},
  ruleSet = 1,

  volumeFluxMode     = 'node', 
  properties = { material = 'materialH' },
} 

PhysicalMethod { 
  id       = 'Interface',
  typeName = 'HydroFemPhysics.Interface', 
  type     = 'fem',

  mesh = 'mesh',
  elementGroups = {'Tri_coh_gcep_rock','Quad_cohp_iniopen','Quad_coh_gcep_rock','Quad_cohp_NF',}, 
  materials = {'interfaceFlow'},
  permeabilityUpdate = true,
  --boundaryConditions={'bfm',},
  volumeFluxMode     = 'node', 
  ruleSet = 1,

  properties = { material = 'materialH' },
} 

 dofile('$SIMULATIONDIR/AuxData.lua') 

local solverOptionsGeo = {                                        
  type               = 'transient nonlinear',                  
  timeMax            = 1,
  timeInitIncrement  = 1.000E+00,                              
  timeMinIncrement   = 1.000E-05,                              
  timeMaxIncrement   = 1,
  iterationsMax      = 15,                                     
  eulerTheta         = 1.000E+00,                              
  attemptMax         = 10,                                     
  incTimeMax         = 2,                                      
  tolerance          ={mechanic =1.000E-05,hydraulic =1e-5,},  
  Increment_Time_Factor = 1.5,                                 
  Frequency = 1,                                               
  --geostaticType      = 'fixedNode',                          
}                                                            

local solverOptions = {                                        
  type               = 'transient nonlinear',                  
  timeMax            = 86400,
  timeInitIncrement  = 60,                              
  timeMinIncrement   = 1.000E-05,                              
  timeMaxIncrement   = 2000,
  iterationsMax      = 15,                                     
  eulerTheta         = 1.000E+00,                              
  attemptMax         = 10,                                     
  incTimeMax         = 2,                                      
  tolerance          ={mechanic =1.000E-05,hydraulic =1e-5,},  
  Increment_Time_Factor = 2.0,                                 
  Frequency = 1,                                               
  --geostaticType      = 'fixedNode',                          
}                                                              

-------------------------------------------------------------  
--  Adaptative Time function                                   
-- Autor: Julio Rueda                                          
-------------------------------------------------------------  
function adaptativeTime(dt,err,conv,niter,nbroke,solverOptions)
  local TolPhy = solverOptions.tolerance["hydraulic"]
  local dtol = (TolPhy-err)/TolPhy*100                         
  local ftime = solverOptions.Increment_Time_Factor            
  if(conv)then                                                 
     if(dtol < 15 or niter > 5) then                           
        dt = dt                                                
     elseif(ftime*dt >= solverOptions.timeMaxIncrement ) then  
        dt = solverOptions.timeMaxIncrement                    
     else                                                     
        dt = ftime*dt                                         
     end                                                      
  else                                                        
     if(nbroke - 1 < solverOptions.attemptMax) then           
        dt = dt/4                                              
     elseif(nbroke > solverOptions.attemptMax)then            
        dt = solverOptions.timeMinIncrement                    
     end                                                      
  end                                                         
  io.print(tr('DTol = %1 %') :num(dtol))                    
                                                              
  return dt                                                   
 end                                                           

-------------------------------------------------------------  
--  Adaptative Processnction                                   
-- Autor: Julio Rueda                                          
-------------------------------------------------------------  
 function ProcessScript()
  dofile('$SCRIPTS/fileLib.lua')  
  fileLib.dirCreate('$SIMULATIONDIR/out') 
  local file = io.prepareMeshFile('mesh', '$SIMULATIONDIR/out/$SIMULATIONNAME', 'nf', {'P'},{'An'},{allstates = true, split = true,  })
  initPore()           
  local solver = fem.init({'Interface','Geostatic'}, 'solver', solverOptionsGeo) 
  local dt = solverOptionsGeo.timeInitIncrement 
  dt=fem.step(solver, dt)  
  setCurrentTime(0.0)  
  io.addResultToMeshFile(file, 0.0)   
  local wfileSta = io.open(translatePath('$SIMULATIONDIR/out/$SIMULATIONNAME.sta'),"w+")                                                
  wfileSta:write( ' SUMMARY OF JOB INFORMATION:', '\n')               
  wfileSta:write( '   INC ATT EQUIL  TOTAL     INC OF', '\n')          
  wfileSta:write( '           ITERS  TIME     TIME/LPF', '\n')          
                                                                              
  --==============================                                             
  -- start running                                                            
  --==============================                                              
  local solver = fem.init({'Rock','Interface'}, 'solver', solverOptions)                                                                      
  local wfile1 = io.open(translatePath('$SIMULATIONDIR/out/Results.dat'),"w+")
  io.output(wfile1) 
  local vFTm = {0.0, 0.0, 0.0, 0.0,0.0, 0.0} 
  local dt = solverOptions.timeInitIncrement               
  local FinalTime =  solverOptions.timeMax                 
  local Time      = 0                                      
  local LastStep  = false                                  
  local countP = 0                                         
  local conv  = true                                       
  local countAttemp = 0                                    
  local frequency = solverOptions.Frequency                
  while (Time < FinalTime) do                              
	    if (Time + dt > FinalTime) then 
		   dt   = FinalTime - Time
		end 
		UpdateOpening()
	    local newdt,err,conv,niter = fem.step(solver, dt, true) 
	    wfileSta:write(string.format('%s%5d%s%2d%s%5d%s%11.4e%s%11.4e\n',' ', countP,' ',countAttemp,' ',niter,' ',Time,' ', dt  ))   	        
	    if(conv)then                                           
        -- Adjust time to guarantee that the last iteration will be on the     
        -- requested final simulation time                                     
		   Time = Time +  dt                                                  
		   io.print(tr('Time = %1 s') :num(Time)) 	                        
		   countAttemp = 0                                                        
	       countP = countP + 1                                                    
	       vFTm = volumeReactFlux (countP, Time, wfile1, dt, vFTm) 
	       if( countP % frequency == 0 or Time==FinalTime)then                        
	          io.addResultToMeshFile(file, Time)                                  
	       end                                                                    
		   dt = adaptativeTime(dt,err,conv,niter,countAttemp,solverOptions)   
	    else                                                                           
         setCurrentTime(Time)                                  countAttemp = countAttemp + 1                                           
         dt = adaptativeTime(dt,err,conv,niter,countAttemp,solverOptions)        
     end	                                                                        
	                                                                                
  end                                                                            
  io.closeMeshFile(file)                      
  io.close(wfileSta)                              
 end                                               
