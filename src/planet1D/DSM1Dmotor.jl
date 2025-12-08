@doc raw"""
     DSM1Dmotor.jl prepares the vertical grids, construt T and H (and gravity+wind in the near future) 
     and preforms LU factorisation of (ω² T - H) for the DSM computation.
  
    This should deal with GPU/MPI optimisation and parallelisation, 
    i.e. we need to write the copy for Metals.jl and CUDA.jl separately.

"""




myVerticalGrid = DSM1D.VerticalGridStructure(DSM1D.my1DDSMmodel,DSM1D.dsm1Dconfig.re,maximum(DSM1D.input.ωᵣ))

