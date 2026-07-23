module neutrinoFlux
    #  (c) Yael Deniz and Isabel Goos 2026


    using CSV,DataFrames,Interpolations
    using PythonCall,Downloads,CodecZlib,DelimitedFiles
    include("neutrinoFlux/config.jl")
    #include("neutrinoFlux/convert_txt_to_csv.jl")
    include("neutrinoFlux/generate_neutrino_flux.jl")
    export read_neutrino_flux_table, produce_neutrino_flux, set_default_neutrino_flux!

end
