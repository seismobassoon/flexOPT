
"""
    (c) Yael Deniz and Isabel Goos 2026

    read_neutrino_flux_table(filename::String, nEbins::Int, nθbins::Int, has_header::Bool)

Read a neutrino flux CSV file from the `../data` directory 
and create a neutrino flux array for each neutrino flavor and type.
The amount of tau neutrinos and antineutrinos is negligible.

# Arguments
- `filename`: The name of the data file (without the `.csv` extension).
- `nEbins`: The number of energy bins.
- `nθbins`: The number of zenith angle bins.
- `has_header`: Is true if the dataset has a header. 

# Returns
A tuple of four `nEbins` × `nθbins` Matrices representing the flux for:
1. `νe` - Electron neutrino
2. `νμ` - Muon neutrino
3. `antiνe` - Electron antineutrino
4. `antiνμ` - Muon antineutrino
"""


function read_neutrino_flux_table(params::Dict)
    # have a look on function makeOPTsemiSymbolic(params::Dict)
    @unpack filename, nEbins, nθbins = params
    # see how to deal with the case where you don't need smth3
    #bin_centers, nuflux_νe_interp, nuflux_νμ_interp, nuflux_antiνe_interp, nuflux_antiνμ_interp, energies, nuflux_νe, nuflux_νμ, nuflux_antiνe, nuflux_antiνμ = read_neutrino_flux_table(filename,nEbis, nθbins, has_header)
    output = @strdict .....
    return output
end

function read_neutrino_flux_table(filename::String, nEbins::Int, nθbins::Int, has_header::Bool; data_dir= joinpath(@__DIR__, "..", "data"))

    # construct path
    #data_dir = joinpath(@__DIR__, "..", "data")
    csv_path = joinpath(data_dir, "$(filename).csv")

    # read the data (has_header=false ensures the first row isn't taken to be the description)

    if !isfile(csv_path)
        txt_path = joinpath(data_dir, "$(filename).txt")
        if isfile(txt_path)
            convert_txt_to_csv(filename;data_dir=data_dir)
        else
            @error "no .txt no .csv found"
        end
    end
    raw_df = CSV.read(csv_path, DataFrame, header=has_header)


   # ensure the number of rows matches the expected grid size
    expected_rows = nEbins * nθbins
    if nrow(raw_df) != expected_rows
        error("Dimension mismatch: CSV has $(nrow(raw_df)) rows, but nEbins * nθbins = $expected_rows.")
    end

    # convert DataFrame to a 2D Matrix, then reshape to a 3D grid
    flux_3d = reshape(Matrix(raw_df), nEbins, nθbins, 5)
 
    # Slice the 3D grid to extract the 2D planes for each specific flavor
    energies      = flux_3d[:, 1, 1]
    nuflux_νμ     = flux_3d[:, :, 2]
    nuflux_antiνμ = flux_3d[:, :, 3]
    nuflux_νe     = flux_3d[:, :, 4]
    nuflux_antiνe = flux_3d[:, :, 5]

    # Compute values at bin centers
    bin_centers, nuflux_νμ_interp = interpolate_flux_at_bin_centers(energies, nuflux_νμ, true)
    -, nuflux_antiνμ_interp = interpolate_flux_at_bin_centers(energies, nuflux_antiνμ, true)
    -, nuflux_νe_interp     = interpolate_flux_at_bin_centers(energies, nuflux_νe, true)
    -, nuflux_antiνe_interp = interpolate_flux_at_bin_centers(energies, nuflux_antiνe, true)

    return bin_centers, nuflux_νe_interp, nuflux_νμ_interp, nuflux_antiνe_interp, nuflux_antiνμ_interp, energies, nuflux_νe, nuflux_νμ, nuflux_antiνe, nuflux_antiνμ

end

"""
    interpolate_flux_at_bin_centers(energies::Vector{Float64}, flux::Matrix{Float64}, energies_in_log::Bool)

Calculates energy bin centers, interpolates the log(flux) at log(energies).

# Arguments
- `energies`: The bin edges for the energy axis.
- `flux`: The flux data matrix.
- `energies_in_log`: If `true`, it assumes the input `energies` are evenly spaced in log-scale 
  andit calculates geometric centers. Otherwise, it calculates arithmetic centers.

# Returns
A tuple containing:
1. `bin_centers`: The calculated centers of the energy bins.
2. `interpolated_flux`: The flux evaluated at the bin centers.
"""


function interpolate_flux_at_bin_centers(energies::Vector{Float64}, flux::Matrix{Float64}, energies_in_log::Bool)

    # create 1D coordinate array for the y axis
    angles_indices = 1.0:size(flux, 2)

    # calculate the bin centers    
    if energies_in_log == true
        bin_centers = sqrt.(energies[1:end-1] .* energies[2:end]) 
    else
        bin_centers = 0.5 .* (energies[1:end-1] .+ energies[2:end])
    end

    # we interpolate log(flux) vs log(energies)
    interpolators = [interpolate(log.(energies), log.(col), FritschCarlsonMonotonicInterpolation()) for col in eachcol(flux)]
    # evaluate the interpolators at the log of the bin centers, then exponentiate back
    interpolated_flux = exp.(hcat([interp.(log.(bin_centers)) for interp in interpolators]...))   
    return bin_centers, interpolated_flux

end

function produce_neutrino_flux_table(model::String)

    #if model = "Daemonflux" etc

    return Nothing

end

