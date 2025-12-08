@doc raw"""
     planet1Dmotor.jl prepares the vertical grids, construt T and H (and gravity+wind in the near future) 
     and preforms LU factorisation of (ω² T - H) for the DSM computation.
  
    This should deal with GPU/MPI optimisation and parallelisation, 
    i.e. we need to write the copy for Metals.jl and CUDA.jl separately.

"""




myVerticalGrid = planet1D.VerticalGridStructure(planet1D.my1DDSMmodel,planet1D.dsm1Dconfig.re,maximum(planet1D.input.ωᵣ))

