module stat_utils
export energy_resolution, angular_resolution

# Function that gives energy resolution RMS(E)
# EnuT [GeV]
# E₁, E₂ [%/100]
energy_resolution(EnuT, E₁, E₂) = EnuT .* (E₁ .+ E₂ ./ sqrt.(EnuT))

# Function that gives angular resolution RMS(zenith)
# EnuT [GeV]
# θ₁, θ₂ [rad]
angular_resolution(EnuT, θ₁, θ₂) = θ₁ .+ θ₂ ./ sqrt.(EnuT)

end





