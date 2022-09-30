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

  properties = { material = 'materialM' },
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
  properties = { material = 'materialM' },
} 

PhysicalMethod { 
  id       = 'Interface',
  typeName = 'HydroFemPhysics.Interface', 
  type     = 'fem',

  mesh = 'mesh',
  elementGroups = {'cohp_iniopen','cohp_nf_7','cohp_nf_8','cohp_nf_9','cohp_nf_10','cohp_nf_11','cohp_nf_12','cohp_nf_13','cohp_nf_14','cohp_nf_15','cohp_nf_16','cohp_nf_17','cohp_nf_18','cohp_nf_19','cohp_nf_20','cohp_nf_21','cohp_nf_22','cohp_nf_23','cohp_nf_24','cohp_nf_25','cohp_nf_26','cohp_nf_27','cohp_nf_28','cohp_nf_29','cohp_nf_30','cohp_nf_31','cohp_nf_32','cohp_nf_33','cohp_nf_34','cohp_nf_35','cohp_nf_36','cohp_nf_37','cohp_nf_38','cohp_nf_39','cohp_nf_40','cohp_nf_41','cohp_nf_42','cohp_nf_43','cohp_nf_44','cohp_nf_45','cohp_nf_46','cohp_nf_47','cohp_nf_48','cohp_nf_49','cohp_nf_50','cohp_nf_51','cohp_nf_52','cohp_nf_53','cohp_nf_54','cohp_nf_55','cohp_nf_56','cohp_nf_57','cohp_nf_58','cohp_nf_59','cohp_nf_60','cohp_nf_61','cohp_nf_62','cohp_nf_63','cohp_nf_64','cohp_nf_65','cohp_nf_66','cohp_nf_67','cohp_nf_68','cohp_nf_69','cohp_nf_70','cohp_nf_71','cohp_nf_72','cohp_nf_73','cohp_nf_74','coh_gcep_rock',}, 
  materials = {'interfaceFlow'},
  permeabilityUpdate = true,
  --boundaryConditions={'bfm',},
  volumeFluxMode     = 'node', 
  ruleSet = 1,

  properties = { material = 'materialM' },
} 

dofile('$SIMULATIONDIR/AuxData.lua') 
                                                          
s_day = 1/86400
local solverOptions = {                                        
  type               = 'transient nonlinear',                  
  timeMax            = 86400*365*s_day,
  timeInitIncrement  = 1*86400*s_day,                              
  timeMinIncrement   = 1*86400*s_day,                              
  timeMaxIncrement   = 1*86400*s_day,
  iterationsMax      = 5,                                                                  
  attemptMax         = 10,                                     
  incTimeMax         = 2,                                      
  tolerance          ={hydraulic =1e-4},  
  Increment_Time_Factor = 1.5,                                 
  Frequency = 10,                                               
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
      elseif(dt/4 < solverOptions.timeMinIncrement)then
        dt = solverOptions.timeMinIncrement 	  
     end                                                      
  end                                                         
  io.print(tr('Erro = %1 %') :num(err))                    
                                                              
  return dt                                                   
 end                                                           

-------------------------------------------------------------  
--  Adaptative Processnction                                   
-- Autor: Julio Rueda                                          
-------------------------------------------------------------  
 function ProcessScript() 
  dofile('$SCRIPTS/fileLib.lua')
  fileLib.dirCreate('$SIMULATIONDIR/out')
  LengthAreaVolume()
  --[[   ]]
  local wfileMsh = io.open(translatePath('$SIMULATIONDIR/out/$SIMULATIONNAME.post.msh'), "w+")                                       
  local wfileRes = io.open(translatePath('$SIMULATIONDIR/out/$SIMULATIONNAME.post.res'),"w+")                                        
  initPore()           
  
  local solver = fem.init({'Interface','Geostatic'}, 'solver', solverOptions) 
  dt=fem.step(solver, 1)  
  setCurrentTime(0.0)  
  PostMesh (wfileMsh)              
  io.close(wfileMsh)  
  UpdateOpening()  
  PostRest (0, 0.0, wfileRes)       
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
  local FinalTime = solverOptions.timeMax                  
  local Time      = 0                                      
  local LastStep  = false                                  
  local countP = 0                                         
  local conv  = true                                       
  local countAttemp = 0                                    
  local frequency = solverOptions.Frequency                
  while (Time < FinalTime) do      
        -- Adjust time to guarantee that the last iteration will be on the     
        -- requested final simulation time   
	    if (Time + dt > FinalTime) then dt   = FinalTime - Time end 
	    local newdt, conv = fem.step(solver, dt)
	    --wfileSta:write(string.format('%s%5d%s%2d%s%5d%s%11.4e%s%11.4e\n',' ', countP,' ',countAttemp,' ',niter,' ',Time,' ', dt  ))   	        
	    if(conv)then                                                                               
		   Time = Time +  dt		   
		   io.print(tr('Time = %1 day') :num(Time)) 	                        
		   countAttemp = 0                                                        
	       countP = countP + 1 
           UpdateOpening()        
	       vFTm = volumeReactFlux (countP, Time, wfile1, dt, vFTm) 
	       if( countP % frequency == 0 or Time==FinalTime)then                        
	          PostRest(countP, Time, wfileRes)                                
	       end  
        else
		    countAttemp = countAttemp + 1
		    --if(countAttemp > solverOptions.attemptMax)then break end                                         
        end	 
        dt = newdt  		                                                                            
  end  	                                                                        
                                                                           
  io.close(wfileRes)                          
  io.close(wfileSta)
   --]]   
 end                                               
