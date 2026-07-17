
# some constants which can be useful for global Earth

correlationLengthDefault=(20e3,20e3,20e2) # not yet fully understood this for DIV interpolation
epsilon2Default =1.;


giveMaskPX(Xs,Ys)=giveMaskPX(Xs,Ys,nothing)

function giveMaskPX(Xs,Ys,Zs)

    minX, maxX, nX = @unpack Xs
    minY, maxY, nY = @unpack Ys

    if Zs === nothing

        mask,(pm,pn),(xi,yi) = DIVAnd_rectdom(range(minX,stop=maxX,length=nX),range(minY,stop=maxY,length=nY))
        return (mask=mask,Ps=(pm,pn),XIs=(xi,yi))
                                            
    else
        minZ, maxZ, nZ = @unpack Zs
        mask,(pm,pn,po),(xi,yi,zi) = DIVAnd_rectdom(range(minX,stop=maxX,length=nX),range(minY,stop=maxY,length=nY),range(minZ,stop=maxZ,length=nZ))
        return (mask=mask,Ps=(pm,pn,po),XIs=(xi,yi,zi))
    end
end

interpolateField(Xs,Ys,Xnode,Ynode;correlationLength=correlationLengthDefault,epsilon2=epsilon2Default) = interpolateField(Xs,Ys,nothing,Xnode,Ynode,nothing;correlationLength=correlationLength,epsilon2=epsilon2)

function interpolateField(Xs,Ys,Zs,XXnodes;correlationLength=correlationLengthDefault,epsilon2=epsilon2Default) 
    masks = giveMaskPX(Xs,Ys,Zs)
    if Znode === nothing
        XXnodes = (Xnode,Ynode)
    else
        XXnodes = (Xnode,Ynode,Znode)
    end
    fieldCartesien,_ =DIVAndrun(masks.mask,masks.Ps,masks.XIs,XXnodes,rawField,correlationLength,epsilon2)
    return fieldCartesien
end