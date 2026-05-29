function getModelPoints(model,pointsInTime,timeMarching)
    fakeNt = 1
    if timeMarching
        fakeNt = pointsInTime+1
         # Nx, Ny etc thing. Nt is also mentioned and it should be the last element!
    end
    modelPoints = (size(model)...,fakeNt)
    return modelPoints 
end

