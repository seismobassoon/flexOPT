
"""
    CubicSpline

Holds the knots (k) and coefficient matrix (c) for a SciPy-generated cubic spline.
"""
struct CubicSpline
    k::Vector{Float64}  # spline knots 
    c::Matrix{Float64}  # spline polynomial coefficients
end

"""
    read_neutrino_cross_sections_info(filename::String)

Reads cross-section data from a file in the `../data` directory and 
returns a tuple with the information for all 6 interesting neutrino interaction channels.
"""
function read_neutrino_cross_sections_info(filename::String)

    # construct path
    data_dir  = joinpath(@__DIR__, "..", "data")
    data_path = joinpath(data_dir, "$(filename)")

    # read the data
    if endswith(filename, ".json")
        data = JSON.parsefile(data_path)
        # parse JSON dictionary 
        cross_section_splines = Dict{String, CubicSpline}()
        for (key, val) in data
            k = Float64.(val["k"])
            c_list = val["c"]                
            c = Float64.(hcat(c_list...))'   
            cross_section_splines[key] = CubicSpline(k, c)
        end
        # extract specific neutrino interaction channels
        νe_CC  = cross_section_splines["nue_CC_logE"]
        νμ_CC  = cross_section_splines["num_CC_logE"]
        ντ_CC  = cross_section_splines["nut_CC_logE"]
        antiνe_CC = cross_section_splines["anue_CC_logE"]
        antiνμ_CC = cross_section_splines["anum_CC_logE"]
        antiντ_CC = cross_section_splines["anut_CC_logE"]
    else
        error("Unsupported file format: '$filename'.")
    end

    return νe_CC, νμ_CC, ντ_CC, antiνe_CC, antiνμ_CC, antiντ_CC

end 

"""
    evaluate_cubicspline(s::CubicSpline, x::Float64)

Evaluate the cubic spline `s` at a given coordinate `x`.

# Arguments
- `s`: The struct containing the spline knots and coefficients.
- `x`: The value at which to evaluate the spline.

# Returns
- `Float64`: The interpolated value.
"""
function evaluate_cubicspline(s::CubicSpline, x::Float64)

    k = s.k
    c = s.c

    # search in which interval 'x' falls into
    i = searchsortedlast(k, x)

    # if x is beyond the last knot, use the final valid polynomial segment
    if i >= length(k)
        i = length(k) - 1
    # if x is below the first knot, use the first valid polynomial segment
    elseif i < 1
        i = 1
    end

    # distance from the left knot of the current interval
    dx = x - k[i]

    # evaluate the cubic polynomial using Horner's method
    return ((c[1, i] * dx + c[2, i]) * dx + c[3, i]) * dx + c[4, i]
end



