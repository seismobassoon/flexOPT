module Neutrino_Cross_Sections

    using JSON
    using CairoMakie

    include("neutrinoCrossSections/generate_neutrino_cross_sections.jl")

    export CubicSpline, CrossSectionSplines, NeutrinoCrossSectionTable
    export read_neutrino_cross_sections_info, evaluate_cubicspline
    export evaluate_neutrino_cross_sections
    export plot_neutrino_cross_sections

end
