using FileIO, Images,ColorSchemes

Base.retry_load_extensions()

function rescalingImage(image,new_height,new_width)
    new_size = (new_height,new_width)
    image_rescaled = imresize(image, new_size);
    return image_rescaled
end

function getColorPalette(cmap::String)
    newcmap=getproperty(ColorSchemes,Symbol(cmap))
    return newcmap
end

function regenerataionColorMap(colorbar,values,Ncolor)
    # regeneration of colorbar and values
    nSegments=length(colorbar)-1
    nColors, indices=reinterpolateArrayMembers(Ncolor,nSegments)

    newColorBar = Array{RGB}(undef,Ncolor)
    newValues = Array{Number}(undef,Ncolor)
    nColors[end]=nColors[end]-1
    for i in 1:nSegments
        Δcolor = (colorbar[i+1]-colorbar[i])/nColors[i]
        Δvalue = (values[i+1]-values[i])/nColors[i]
        for j in indices[1,i]:indices[2,i]
            newColorBar[j]=colorbar[i]+Δcolor*(j-indices[1,i])
            newValues[j]=values[i]+Δvalue*(j-indices[1,i])

        end
    end
    return newColorBar, newValues
end

function color2Float(image,colorbar,values)
    #region Introduction and small vector analysis explanation
    # here we would like to compute the nearest color that is on the colorbar defined 
    # by several points in 3D RGB space

    # Let us define \mathbf{a} and \mathbf{b} the end points of a segment and
    # assume that the color bar is defined linearly between the two points
    # 
    # we have \mathbf{y} as the given pixel color and we need to find the nearest color \mathbf{x}
    #
    # \mathbf{x} = \mathbf{a} + \alpha \mathbf{p}
    # with
    # \mathbf{p} = \mathbf{b} - \mathbf{a}
    # and 
    # \alpha \in [0,1] for the intermediate segments otherwise <0 for the beginning of the colorbar and >1 for the end
    #
    # now we have the orthogonality 
    # [\mathbf{y}-\mathbf{x}]^T  \mathbf{p} = 0
    # which is true for the case \mathbf{y}=\mathbf{x}
    # 
    # yielding
    #
    # \alpha = \frac{[\mathbf{y}-\mathbf{a}]^T \mathbf{p}}{|\mathbf{p}|^2}
    #
    #endregion
    height, width = size(image)[1:2]
    nSegments = length(colorbar)-1
    α = Array{Float64,3}(undef, height, width, nSegments)
    l² = Array{Float64,4}(undef, 3,height, width, nSegments)
    imageScalar = Array{Float64,2}(undef, height, width)


    for i in 1:nSegments
        a=[red(colorbar[i]) ; green(colorbar[i]) ; blue(colorbar[i])]        
        b=[red(colorbar[i+1]) ; green(colorbar[i+1]) ; blue(colorbar[i+1])]
        p = b-a
        p²=transpose(p)*p
        one_p² = 1/p²


        for j in 1:width
            for k in 1:height
                y = [red(image[k,j]); green(image[k,j]); blue(image[k,j])]
                y_a= y-a
                y_b= y-b
                α[k,j,i] =  transpose(y_a)*p*one_p²
                l²[1,k,j,i] = transpose(y_a)*y_a
                l²[2,k,j,i] = transpose(y_b)*y_b
                l²[3,k,j,i] = l²[1,k,j,i] + l²[2,k,j,i] 
            end
        end
    end
  
    
    for j in 1:width
        for k in 1:height
            αtmp=nothing

            i = argmin(l²[3,k,j,:])
         
            valueA=values[i]
            valueB=values[i+1]
            valueB_A = valueB-valueA

            αtmp = α[k,j,i]
            if αtmp <= 1.0 && αtmp >= 0.0
                imageScalar[k,j] = valueA+valueB_A*αtmp
            elseif l²[1,k,j,i] <= l²[2,k,j,i]
                imageScalar[k,j] = valueA
            else
                imageScalar[k,j] = valueB
            end
            
        end
    end
    return imageScalar

end

# This takes super long!
function color_to_value_old(color::RGB,colorbar,values)
    # Find the closest color in the colorbar
    distances = [colordiff(color, cb) for cb in colorbar]
    idx = argmin(distances)  # Index of the closest color
    return values[idx] # Map index to interpolated value
end

function myReadPhoto(file)
    data=float(load(file))
    return data
end
