# here we define the number of flavours and inter-angles between them

DEFAULT_NEUTRINO_FLAVOURS = Ref(:NormalThree)

function set_default_flavour!(name::Symbol;roundingU=0.01,roundingH=0.00001)
    # if the user does not round U and H, s/he should put 0.e0

    DEFAULT_NEUTRINO_FLAVOURS[]=name
    return makingFlavourConfiguration(name,roundingU,roundingH)
end


function makingFlavourConfiguration(name::Symbol,roundingU,roundingH)
    
    if name == :NormalThree
        osc = OscillationParameters(3)
        setθ!(osc, 1=>2, 0.59)
        setθ!(osc, 1=>3, 0.15)
        setθ!(osc, 2=>3, 0.84)
        setδ!(osc, 1=>3, 3.86)
        setΔm²!(osc, 2=>3, -2.523e-3)
        setΔm²!(osc, 1=>2, -7.39e-5)
        U = PMNSMatrix(osc)
        H = Hamiltonian(osc)
    else
        error("Unknown flavour configuration: $name")
    end

    if roundingU !== 0.e0
        U = roundExt.(U, 0.01)
    end

    if roundingH !== 0.e0
        H = roundExt.(H, 0.00001)
    end


    return U, H
end


U, H = set_default_flavour!(DEFAULT_NEUTRINO_FLAVOURS[])