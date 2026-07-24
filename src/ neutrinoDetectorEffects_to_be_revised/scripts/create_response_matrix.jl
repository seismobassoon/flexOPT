import Pkg
using Distributions
using Plots
include("../utils/stat_utils.jl")
using .stat_utils

function create_response_matrix(EnuTs, EnuRs, E₁, E₂, θ₁, θ₂)
    # compute energy and angular resolutions
    energy_resolutions  = energy_resolution(EnuTs, E₁, E₂)
    angular_resolutions = angular_resolution(EnuTs, θ₁, θ₂) 
    # create normal distributions
    energy_normal_distr = Normal.(EnuTs, energy_resolutions)
    angle_normal_distr  = Normal.(EnuTs, angular_resolutions)
    # compute elements of the response matrix
    energy_TtoR = pdf.(energy_normal_distr, EnuRs')
    angle_TtoR  = pdf.(angle_normal_distr, EnuRs')
    return energy_TtoR, angle_TtoR
end

EnuTs = logrange(1, 100, 100)
EnuTs_centers = sqrt.(EnuTs[1:end-1] .* EnuTs[2:end])
energy_TtoR, angle_TtoR = create_response_matrix(EnuTs_centers, EnuTs_centers, 0.5, 0.5, 1.0, 1.0)
size(EnuTs_centers)
size(energy_TtoR)
size(angle_TtoR)
plot(EnuTs_centers, angle_TtoR)












# Evaluate the PDF at x = 0.0 (the peak)
#y = pdf(d, 0.0) 
#println(y) # O

