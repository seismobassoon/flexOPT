using Test,BenchmarkTools,CairoMakie

ParamFile = "../config/testparam.csv"

include("../src/planet1D.jl")
using .DSM1D

#println(DSM1D.my1DDSMmodel.solid_or_fluid)
#print(DSM1Dconfig)
#print("dsm1Dconfig: ", DSM1D.dsm1Dconfig)
arrayRadius, arrayParams=DSM1D.compute1DseismicParamtersFromPolynomialCoefficients(DSM1D.my1DDSMmodel,10)
f=Figure()
#lines(f[1,1],arrayRadius*DSM1D.my1DDSMmodel.averagedPlanetRadiusInKilometer, arrayParams.ρ, markersize=1,color=:red)
lines(f[1,1],arrayRadius*DSM1D.my1DDSMmodel.averagedPlanetRadiusInKilometer, arrayParams.ρ,color=:red)
scatter!(f[1,1],arrayRadius*DSM1D.my1DDSMmodel.averagedPlanetRadiusInKilometer, arrayParams.ρ, markersize=3,color=:blue)
display(f)

display(DSM1D.input.ωᵣ)
@testset "First series of tests" verbose=true begin
    #@test DSM1D.greet() === nothing
    @test DSM1D.dsm1Dconfig.re === 0.01
    #@test DSM1D.dsm1Dconfig.omegai === 0.01
    #println("DSM1D.input.omegai: ", DSM1D.input.omegai)
end