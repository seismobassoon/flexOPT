module Neutrino_Cross_Sections


    using JSON


    include("neutrinoCrossSections/generate_neutrino_cross_sections.jl")

    export read_neutrino_cross_sections_info, evaluate_cubicspline

    
end