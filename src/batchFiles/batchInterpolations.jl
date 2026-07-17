
# some constants which can be useful for global Earth

correlationLengthDefault=(20e3,20e3,20e2) # not yet fully understood this for DIV interpolation
epsilon2Default =1.;


giveMaskPX(Xs,Ys)=giveMaskPX(Xs,Ys,nothing)

function giveMaskPX(Xs,Ys,Zs)
    if Zs === nothing

        minZ=0.0
        maxZ=0.0
        tmpX=correlationLength[1]
        tmpY=correlationLength[2]

        #correlationLength=(tmpX,tmpY)

        return mask,(pm,pn),(xi,yi) = DIVAnd_rectdom(range(minX,stop=maxX,length=nX),
                                                range(minY,stop=maxY,length=nY));
    else

        return mask,(pm,pn,po),(xi,yi,zi) = DIVAnd_rectdom(range(minX,stop=maxX,length=nX),
                                                range(minY,stop=maxY,length=nY),
                                                range(minZ,stop=maxZ,length=nZ));
    end
end


function interpolate2D3DField(;correlationLength=correlationLengthDefault,epsilon2=epsilon2Default) 
    mask, Ps, XIs = giveMaskPX()

fieldCartesien,_ =DIVAndrun(mask,(pm,pn),(xi,yi),(Xnode,Ynode),diffρField,correlationLength,epsilon2)
end