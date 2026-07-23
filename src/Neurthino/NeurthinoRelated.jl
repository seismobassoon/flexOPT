#= these functions here are the beautiful fruits of Estelle Salomé's 2-month internship in 2025 =#




function solveQuadraticEquation(a,b,c)
    Δ = b^2 - 4*a*c

    if Δ>0
        x1 = ((-b - sqrt(Δ))/(2*a))
        x2 = ((-b + sqrt(Δ))/(2*a))
    else
        x1 = -b/(2 *a)
        x2 = x1
    end
    return x1, x2

end

function posOrNeg(cos_θ, sign = :positive)
    #to choose if we want a positive or a negative angle
    if sign == :positive
        θ = acos.(cos_θ)
    else
        θ = .- acos.(cos_θ)
    end
    return θ
end


function roundExt(x,step)
    #to round the values of H and U
    if x isa Complex
        xreal = round(real(x)/step)*step
        ximag = round(imag(x)/step)*step
        return complex(xreal, ximag)
    else
        return round(x/step)*step
    end
end