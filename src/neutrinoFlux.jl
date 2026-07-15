module neutrinoFlux
    #  (c) Yael Deniz and Isabel Goos 2026


    using CSV
    using DataFrames
    using Interpolations
    include("neutrinoFlux/config.jl")
    include("neutrinoFlux/convert_txt_to_csv.jl")
    include("neutrinoFlux/generate_neutrino_flux.jl")
    export read_neutrino_flux_table

end
