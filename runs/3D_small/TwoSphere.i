[Mesh]
   type = GeneratedMesh
   dim = 3
   nx = 64
   ny = 64
   nz = 64
  
   xmax = 64
   ymax = 64
   zmax = 64

   elem_type = HEX8
[]


[Variables]

   [./cNi]
      order = FIRST
      family = LAGRANGE
      [./InitialCondition]
         type = Micstr
         file_name = mic.nc
      [../]
   [../]

   [./xNi]
      order = FIRST
      family = LAGRANGE
      [./InitialCondition]
         type = Micstr
         file_name = mic.nc
      [../]
   [../]

[]

[Materials]
   [./bulk_FE]
      type = DerivativeParsedMaterial
      args = cNi
      f_name = bulk_FE
      function = '0.5*cNi^2*(1.0-cNi)^2'
   [../]

   [./Interp_FE]
      type = DerivativeParsedMaterial
      args = cNi
      f_name = h_FE
      function = 'if(cNi>=1.0, 1.0, if(cNi<=0.0, 0.0, cNi^3*(10.0+cNi*(6.0*cNi-15.0))))'
   [../]

   [./Interp_c]
      type = DerivativeParsedMaterial
      args = cNi
      f_name = h_c
      function = 'if(cNi>=1.0, 1.0, if(cNi<=0.0, 0.0, cNi))'
   [../]

   [./Ni_diff]
      type = DerivativeParsedMaterial
      args = cNi
      f_name = Ni_diff
      constant_names = 'DNi DPore DSurf'
      constant_expressions = '0.00 1.0 0.0'
      function = 'if(cNi>=1.0, DNi, if(cNi<=0.0, DPore, DNi*cNi+DPore*(1.0-cNi)+DSurf*cNi*(1.0-cNi)))'
   [../]

   [./Ni_mobility]
      type = DerivativeParsedMaterial
      args = cNi
      f_name = c_mob
      constant_names = 'Mob'
      constant_expressions = '1.0'
      function = 'Mob'
   [../]

   [./kappa_c]
      type = GenericConstantMaterial
      prop_names = kappa_c
      prop_values = 0.5
   [../]
   [./rho_c]
      type = GenericConstantMaterial
      prop_names = rho_c
      prop_values = 0.5
   [../]

   [./x_seq]
      type = GenericConstantMaterial
      prop_names = x_seq
      prop_values = 0.999
   [../]
   [./x_peq]
      type = GenericConstantMaterial
      prop_names = x_peq
      prop_values = 1e-6
   [../]
   [./RTV]
      type = GenericConstantMaterial
      prop_names = RT_Vm
      prop_values = 0.001
   [../]

[]


[Kernels]

   [./xTime]
     type = TimeDerivative
     variable = xNi
   [../]
   [./xdiff]
      type = xDiffusion
      variable = xNi
      c_variable = cNi
      D_Ni = Ni_diff
   [../]
   [./xcCouple]
      type = xcCoupling
      variable = xNi
      c_variable = cNi
      x_seq = x_seq
      x_peq = x_peq
      D_Ni = Ni_diff
      h_c  = h_c
   [../]

   [./cTime]
     type = TimeDerivative
     variable = cNi
   [../]
   [./cDiff]
      type = cDiffusion
      variable = cNi
      M_c = c_mob
      kappa_c = kappa_c
   [../]
   [./bulkFE]
      type = bulkFE
      variable = cNi
      M_c = c_mob
      rho_c = rho_c
      f_FE = bulk_FE
   [../]
   [./ChemFE]
      type = ChemFE
      variable = cNi
      x_variable = xNi
      M_c = c_mob
      x_seq = x_seq
      x_peq = x_peq
      RTVm = RT_Vm
      h_FE = h_FE
      h_c = h_c
   [../]
[]

[BCs]
  [./Periodic]
    [./auto]
       variable = 'xNi cNi'
       auto_direction = 'x y'
    [../]
  [../]
[]

[Preconditioning]
   #active = fdp
   #[./fdp]
   # type = FDP
   # full = true
   # petsc_options_iname = '-mat_fd_coloring_err -mat_fd_type'
   # petsc_options_value = '1e-6  ds'
   #[../]
   #[./pbp]
   # type = PBP
   # solve_order = 'cNi xNi'
   # preconditioner = 'ILU AMG'
   # off_diag_row = 'xNi'
   # off_diag_column = 'cNi'
   #[../]
[]

[Executioner]
  type = Transient
  #type = Steady
  solve_type = 'PJFNK'
  #solve_type = 'JFNK'
  num_steps = 1000
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_grmes_restart'
  petsc_options_value = 'hypre    boomeramg  101 '

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.0
    optimal_iterations = 10
    linear_iteration_ratio = 5
  [../]

  [./Adaptivity]
     initial_adaptivity = 0
     refine_fraction = 0.5
     coarsen_fraction = 0.
     max_h_level = 1
  [../]

[]

[Postprocessors]
  [./dofs]
     type = NumDOFs
  [../]
[]

[Outputs]

  [./pgraph]
     type = PerfGraphOutput
     execute_on = 'timestep_end'
     level = 2
     heaviest_branch = true
     heaviest_sections = 7
  [../]

  [./console]
    type = Console
    interval = 1 
  [../]
  [./exodus]
    type = Exodus
    interval = 1
  [../]
[]
