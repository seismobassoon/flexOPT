using Test,BenchmarkTools,CairoMakie

ParamFile = "../config/testparam.csv"

include("../src/planet1D.jl")
using .planet1D

#println(planet1D.my1DDSMmodel.solid_or_fluid)
#print(planet1Dconfig)
#print("dsm1Dconfig: ", planet1D.dsm1Dconfig)
arrayRadius, arrayParams=planet1D.compute1DseismicParamtersFromPolynomialCoefficients(planet1D.my1DDSMmodel,10)
f=Figure()
#lines(f[1,1],arrayRadius*planet1D.my1DDSMmodel.averagedPlanetRadiusInKilometer, arrayParams.ρ, markersize=1,color=:red)
lines(f[1,1],arrayRadius*planet1D.my1DDSMmodel.averagedPlanetRadiusInKilometer, arrayParams.ρ,color=:red)
scatter!(f[1,1],arrayRadius*planet1D.my1DDSMmodel.averagedPlanetRadiusInKilometer, arrayParams.ρ, markersize=3,color=:blue)
display(f)

display(planet1D.input.ωᵣ)
@testset "First series of tests" verbose=true begin
    #@test planet1D.greet() === nothing
    @test planet1D.dsm1Dconfig.re === 0.01
    #@test planet1D.dsm1Dconfig.omegai === 0.01
    #println("planet1D.input.omegai: ", planet1D.input.omegai)
end