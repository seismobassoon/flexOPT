
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

"""
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
function read_neutrino_flux_table(filename::String, nEbins::Int, nθbins::Int, has_header::Bool)

    # construct path
    data_dir = joinpath(@__DIR__, "..", "data")
    csv_path = joinpath(data_dir, "$(filename).csv")

    # read the data (has_header=false ensures the first row isn't taken to be the description)
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

"""
    set_dflux_params(flux_obj)

Initialize a dictionary of Daemonflux parameters set to their default values.

# Arguments
- `flux_obj`: A Python object representing a Daemonflux instance.

# Returns
- `Dict{String, Float64}`: A dictionary mapping parameter names to `0.0`.
"""
function set_dflux_params(flux_obj)

    # Fetch all the Daemonflux parameter names 
    params = pyconvert(Vector{String}, flux_obj.params.known_parameters)

    # Set all parameters to 0 (meaning that their difference with respect to the
    # default parameters is zero)
    params_dict = Dict(p => 0.0 for p in params) 
    
    return params_dict

end

"""
    set_dflux_params(flux_obj, params_dict)

Initialize a dictionary of Daemonflux parameters to those values given in `params_dict`.
Those parameters not mentioned there are set to default values.

# Arguments
- `flux_obj`: A Python object representing a Daemonflux instance.
- `params_dict`: A dictionary with those pairs of parameter names and values which should be modified.
The possible parameter names are K+_158G, K+_2P, K+_31G, K-_158G, K-_2P, K-_31G, n_158G, n_2P, p_158G, p_2P, 
pi+_158G, pi+_20T, pi+_2P, pi+_31G, pi-_158G, pi-_20T, pi-_2P, pi-_31G, GSF_1, GSF_2, GSF_3, GSF_4, GSF_5, GSF_6.
The get_dflux_param_value function helps find specific values.

# Returns
- `Dict{String, Float64}` 
"""
function set_dflux_params(flux_obj, params_dict::Union{Nothing, Dict{String, Any}})

    # Set default Daemonflux parameters
    full_params = set_dflux_params(flux_obj)

    # If params_dict is nothing, just return the defaults without trying to merge
    if isnothing(params_dict)
        return full_params
    end

    # Otherwise, update with pairs from params_dict
    merge!(full_params, params_dict)

    return full_params

end

"""
    get_dflux_param_value(flux_obj, param_name, signed_sigma; max_iterations=10000, max_param=10., tol=0.01)

Find the value of the parameter `param_name` which is `signed_sigma` away from the default value.

# Arguments
- `flux_obj`: A Python object representing a Daemonflux instance.
- `param_name`: The name of the parameter (see the set_dflux_params function for possible parameter names).
- `signed_sigma`
- `max_iterations`: The maximum number of iterations in the Bisection method.
- `max_param`: The maximum search limit for the absolute value of the parameter value.
- `tol`: The relative tolerance on the target distance from the default value.
"""
function get_dflux_param_value(flux_obj, param_name, signed_sigma; max_iterations=10000, max_param=10., tol=0.01)
    
    # Determine the sign of the parameter based on if it should be above or below the default value
    p_sign = sign(signed_sigma)
    # flux_obj.chi2 is symmetric and positive, the search is performed on the positive axis
    sigma = abs(signed_sigma)

    # Define the starting search bounds
    p_low  = 0.
    p_high = max_param

    # Evaluate the function at the upper boundary to ensure a root exists
    chi2_high = pyconvert(Float64, flux_obj.chi2(params=Dict(param_name => p_high)))
    diff_high = sqrt(chi2_high) - sigma
    # If the target is out of bounds, fail immediately 
    if diff_high < 0.0
        error("Target sigma ($signed_sigma) is outside the search bounds. At p = $max_param, chi2 is only $chi2_high. Increase max_p.")
    end

    # Run Bisection search
    for i in 1:max_iterations

        p_mid = (p_low + p_high) / 2
        chi2_mid = pyconvert(Float64, flux_obj.chi2(params=Dict(param_name => p_mid)))
        sigma_diff = sqrt(chi2_mid) - sigma

        # Check if we are within the relative tolerance threshold
        if abs(sigma_diff) <= tol * sigma
            return p_sign * p_mid  
        end

        # Adjust search bounds
        if sigma_diff < 0.0 
            p_low = p_mid
        else
            p_high = p_mid
        end

    end

    error("Bisection failed to converge within $max_iterations iterations. Check your parameters.")

end

function set_hflux_params(nuflux_params::Dict{String, Any})

    loc        = Symbol(nuflux_params["location_hf"])
    season     = Symbol(nuflux_params["season_hf"])
    angles_sep = Symbol(nuflux_params["angles_separation_hf"])
    mountain   = Symbol(nuflux_params["mountain_overburden_hf"])
    solar      = Symbol(nuflux_params["solar_activity_hf"])

    loc_str    = LOCATION_MAPPING[loc]        
    season_str = SEASON_MAPPING[season]
    angles_str = ANGLES_MAPPING[angles_sep]
    mtn_str    = MOUNTAIN_MAPPING[mountain]
    solar_str  = SOLAR_MAPPING[solar]

    filename = "http://www-rccn.icrr.u-tokyo.ac.jp/mhonda/public/nflx2014/$(loc_str)-$(season_str)-$(angles_str)$(mtn_str)-$(solar_str).d.gz"
    return filename
    
end

"""
    produce_neutrino_flux(bin_centers_arrays; model, flux_mode)

Produce the atmospheric neutrino flux for all relevant neutrino types following a given `model`.

If the `model` `Daemonflux` is chosen, compute the atmospheric neutrino fluxes at `IceCube` 
using the Daemonflux model (arXiv:2303.00022).
`Icecube` is, for the moment, the only location for which upgoing neutrinos are computed (July 2026).
`use_calibration` is set to true so that the Daemonflux calibration is used.
If the `model` `Honda` is chosen, fetch the atmospheric neutrino fluxes from 
<http://www-rccn.icrr.u-tokyo.ac.jp/mhonda/public/nflx2014/index.html> (arXiv:1502.03916).
Neither model provides tau-neutrino or -antineutrino fluxes since they are negligible.

# Arguments
- `model::String`: `Daemonflux` (default) or `Honda`.
- `bin_centers_arrays::Tuple{AbstractArray, AbstractArray}`: a 2-tuple containing the bin centers 
for energy and cosθ. The bin centers for cosθ have to be given in increasing order.
- `params_df::Union{Nothing, Dict{String, Any}}`: Daemonflux parameters. By default all Daemonflux parameters
are set to their default values. Only those provided in the `params_df` dictionary are updated to the desired values.
- `flux_mode::String`: `conventional` or `total` (default). `total` means conventional + prompt (only the 
conventional componennt is calibrated).
- `return_uncertainties::Bool`: if true, return neutrino flux uncertainties. Default is false. 
- `location_hf::String`: location for honda flux, options are `kamioka`, `gran_sasso`, `sudbury`, `frejus`, `ino`, 
`south_pole`, `pythasalmi`, `homestake`, `juno`
- `season_hf::String`: season for honda flux, options are `all_year`, `march_may`, `june_aug`, `sept_nov`, `dec_feb`
- `angles_separation_hf::String`: the way the angles are organized in the honda flux, options are 
`variable_ϕ`, `averaged_ϕ`, `averaged_ϕθ`
- `mountain_overburden_hf::String`: `with` or `without`
- `solar_activity_hf::String`: `minimum` or `maximum`

# Returns 
a 4-tuple (or 8-tuple, if `return_uncertainties::Bool` is true) with the following matrices 
in the shape (nθ, nE):
- `NuMu_flux` 
- `Nue_flux`
- `AntiNuMu_flux` 
- `AntiNue_flux`
(followed by the corresponding uncertanties if `return_uncertainties::Bool` is true, in the same order).
"""
function produce_neutrino_flux(nuflux_params::Dict{String, Any})
  
    model = Symbol(nuflux_params["model"])

    if (model === :Daemonflux)
        results = produce_daemonflux_neutrino_flux(nuflux_params)
    elseif (model === :Honda)
        results = produce_honda_neutrino_flux(nuflux_params)
    else 
        error("model $model is not supported. Use :Daemonflux or :Honda.")
    end

    return results

end

"""
    produce_daemonflux_neutrino_flux(nuflux_params::Dict{String, Any})

Produce the atmospheric neutrino flux for all relevant neutrino types using the Daemonflux model.

See produce_neutrino_flux for more information.
"""
function produce_daemonflux_neutrino_flux(nuflux_params::Dict{String, Any})
    
    # Read some inputs
    Ebin_centers, cosθbin_centers = nuflux_params["bin_centers_arrays"]
    flux_mode            = Symbol(nuflux_params["flux_mode"])
    return_uncertainties = Bool(nuflux_params["return_uncertainties"])

    # Import and instantiate Python daemonflux class
    daemonflux = pyimport("daemonflux")
    flux = daemonflux.Flux(location="IceCube", 
                           use_calibration=true, 
                           debug=1)

    # Convert to degrees, as needed by daemonflux
    θbin_centers = string.(rad2deg.(acos.(cosθbin_centers)))

    # Pre-allocate output arrays
    nE = length(Ebin_centers)
    nθ = length(θbin_centers)
    NuMu_flux     = zeros(Float64, nθ, nE)
    Nue_flux      = zeros(Float64, nθ, nE)
    AntiNuMu_flux = zeros(Float64, nθ, nE)
    AntiNue_flux  = zeros(Float64, nθ, nE)
    if return_uncertainties
        NuMu_flux_err     = zeros(Float64, nθ, nE)
        Nue_flux_err      = zeros(Float64, nθ, nE)
        AntiNuMu_flux_err = zeros(Float64, nθ, nE)
        AntiNue_flux_err  = zeros(Float64, nθ, nE)
    end

    # Check mode and define key names depending on the chosen flux_mode 
    if flux_mode !== :total && flux_mode !== :conventional
        error("flux_mode $flux_mode is not a valid option. Use :total or :conventional.")
    end
    prefix = flux_mode === :total ? "total_" : ""
    keys = (
        numu    = prefix * "numu",
        nue     = prefix * "nue",
        anumu   = prefix * "antinumu",
        anue    = prefix * "antinue" 
    )

    # Set daemonflux parameters
    raw_params_df = get(nuflux_params, "params_df", nothing)
    params = (isnothing(raw_params_df) || isempty(raw_params_df)) ? 
             set_dflux_params(flux, Dict{String, Any}()) : 
             set_dflux_params(flux, raw_params_df)
            
    # Energy scaling
    E_scaling = 1.e4 ./ (Ebin_centers .^ 3)
    # Fetch raw daemonflux values
    @views for (i, θ) in enumerate(θbin_centers)
        # flux values
        NuMu_flux[i, :]     .= pyconvert(Vector{Float64}, flux.flux(Ebin_centers, θ, keys.numu,  params=params)) .* E_scaling
        Nue_flux[i, :]      .= pyconvert(Vector{Float64}, flux.flux(Ebin_centers, θ, keys.nue,   params=params)) .* E_scaling
        AntiNuMu_flux[i, :] .= pyconvert(Vector{Float64}, flux.flux(Ebin_centers, θ, keys.anumu, params=params)) .* E_scaling
        AntiNue_flux[i, :]  .= pyconvert(Vector{Float64}, flux.flux(Ebin_centers, θ, keys.anue,  params=params)) .* E_scaling
        if return_uncertainties
            # flux uncertainties
            NuMu_flux_err[i, :]     .= pyconvert(Vector{Float64}, flux.error(Ebin_centers, θ, keys.numu))
            Nue_flux_err[i, :]      .= pyconvert(Vector{Float64}, flux.error(Ebin_centers, θ, keys.nue))
            AntiNuMu_flux_err[i, :] .= pyconvert(Vector{Float64}, flux.error(Ebin_centers, θ, keys.anumu))
            AntiNue_flux_err[i, :]  .= pyconvert(Vector{Float64}, flux.error(Ebin_centers, θ, keys.anue))
        end
    end

    # Prepare dictionary for output
    results = Dict{String, Matrix{Float64}}(
        "NuMu_flux"     => NuMu_flux,
        "Nue_flux"      => Nue_flux,
        "AntiNuMu_flux" => AntiNuMu_flux,
        "AntiNue_flux"  => AntiNue_flux
    )
    if return_uncertainties
        results["NuMu_flux_err"]     = NuMu_flux_err
        results["Nue_flux_err"]      = Nue_flux_err
        results["AntiNuMu_flux_err"] = AntiNuMu_flux_err
        results["AntiNue_flux_err"]  = AntiNue_flux_err
    end

    return results

end

"""
    produce_honda_neutrino_flux(nuflux_params::Dict{String, Any})

Produce the atmospheric neutrino flux for all relevant neutrino types using the honda model.

See produce_neutrino_flux for more information.
"""
function produce_honda_neutrino_flux(nuflux_params::Dict{String, Any})

    # Read parameters and define url from where to fetch the data
    url_honda  = set_hflux_params(nuflux_params)
    angles_sep = Symbol(nuflux_params["angles_separation_hf"])
    
    # Download and decompression
    io_buffer = IOBuffer()
    try
        Downloads.download(url_honda, io_buffer)
    catch e 
        if e isa Downloads.RequestError
            error("Failed to download Honda flux. The URL does not exist or the server is down.\nURL attempted: $url_honda\nDetails: $(e.message)")
        else
            rethrow(e) 
        end
    end
    seekstart(io_buffer)
    
    # Decompress
    decompressed_stream = GzipDecompressorStream(io_buffer)
    
    # Read the text grid into a Matrix{Float64}, omit any lines with text 
    lines = readlines(decompressed_stream)
    close(decompressed_stream)
    numeric_lines = filter(line -> occursin(r"^\s*[0-9.-]", line), lines)
    clean_text_block = join(numeric_lines, "\n")
    matrix_2d = readdlm(IOBuffer(clean_text_block), Float64)
    
    # Reshaping based on angular separation mode
    if angles_sep === :averaged_ϕ

        # Grid layout: 5 columns x 101 energies x 20 cosθ steps
        n_E    = 101
        n_cosθ = 20
        flux_reshaped = reshape(matrix_2d', 5, n_E, n_cosθ)

        # Extract matrices and transpose using optimization views ('.') 
        # to match your required shape without unnecessary memory allocation.
        results = Dict{String, Any}(
            "NuMu_flux"     => collect(flux_reshaped[2, :, :]'),
            "Nue_flux"      => collect(flux_reshaped[4, :, :]'),
            "AntiNuMu_flux" => collect(flux_reshaped[3, :, :]'),
            "AntiNue_flux"  => collect(flux_reshaped[5, :, :]'),
            "Energies"      => flux_reshaped[1, :, 1],
            "cosθs"         => range(-1.0, 1.0, length=20)
        )
        return results

    elseif angles_sep === :variable_ϕ
        # TODO: Implement the reshaping logic for the azimuth-dependent grid
        error("angles_separation_hf mode :variable_ϕ is not yet implemented.")

    elseif angles_sep === :averaged_ϕθ
        # TODO: Implement the reshaping logic for the fully integrated/averaged grid
        error("angles_separation_hf mode :averaged_ϕθ is not yet implemented.")

    else
        error("Invalid angles_separation_hf option: $angles_sep. Use :variable_ϕ, :averaged_ϕ, or :averaged_ϕθ.")
    end

end