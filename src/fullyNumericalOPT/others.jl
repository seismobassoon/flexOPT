function getModelPoints(model,pointsInTime,timeMarching)
    
    fakeNt = 1
    if timeMarching
        fakeNt = pointsInTime+1
        modelPoints = (size(model)...,fakeNt) # Nx, Ny etc thing. Nt is also mentioned and it should be the last element!
    else
        modelPoints = (size(model)...,1)
    end
    return modelPoints 
end

