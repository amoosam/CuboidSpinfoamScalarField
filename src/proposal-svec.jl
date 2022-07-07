using Base: setindex

function next_step_small_single(lat::Lattice,lenx::T,leny::T,
    lenz::T,lent::T,ind::Union{Int,Nothing} = nothing) where {T<:SArray{Tuple{lat_size}, Float64,1,lat_size}}
    shift = 0.5
    if isnothing(ind)
        index = rand(1:Size^4)
    else
        index = ind
    end

    (i,j,k,l) = lat.node2ind[ind]
    lx = setindex(lenx, lx[i] + shift*(2*rand() - 1),i)
    ly = setindex(leny, ly[j] + shift*(2*rand() - 1),j)
    lz = setindex(lenz, lz[k] + shift*(2*rand() - 1),k)
    lt = setindex(lent, lt[l] + shift*(2*rand() - 1),l)

    return (lx,ly,lz,lt)
end

function next_step_small(lenx::T,leny::T,lenz::T,lent::T) where {T<:SArray{Tuple{lat_size}, Float64,1,lat_size}}
    shift = 1.0


    lx = SVector{lat_size}(lenx .+ shift.*((2 * shift) .* rand(lat_size) .- ones(lat_size)))
    ly = SVector{lat_size}(leny .+ shift.*((2 * shift) .* rand(lat_size) .- ones(lat_size)))
    lz = SVector{lat_size}(lenz .+ shift.*((2 * shift) .* rand(lat_size) .- ones(lat_size)))
    lt = SVector{lat_size}(lent .+ shift.*((2 * shift) .* rand(lat_size) .- ones(lat_size)))
    
    return (lx,ly,lz,lt)
end

function getfactor(s::Int=1)
    r = rand(1:2)
    if s == 1
        if r == 1
            return 0.5 + 0.5*rand()
        else
            return 1.0 + rand() 
        end
    else 
        if r == 1
            return 0.5 .+ 0.5 .* rand(s)
        else
            return 1.0 .+ rand(s) 
        end
    end
end


function next_step_small_rescale(lenx::T,leny::T,lenz::T,lent::T) where {T<:SArray{Tuple{lat_size}, Float64,1,lat_size}}

    lenx_new = getfactor().* lenx
    leny_new = getfactor().* leny
    lenz_new = getfactor().* lenz
    lent_new = getfactor().* lent

    return (lenx_new, leny_new, lenz_new, lent_new)
end

function next_step_rescale(lenx::T,leny::T,lenz::T,lent::T) where {T<:SArray{Tuple{lat_size}, Float64,1,lat_size}}
#rescales all the lengths in all directions by different scalings between a and b.

    lenx_new = getfactor(lat_size).* lenx
    leny_new = getfactor(lat_size).* leny
    lenz_new = getfactor(lat_size).* lenz
    lent_new = getfactor(lat_size).* lent
    return (lenx_new, leny_new, lenz_new, lent_new)
end

function next_step_rescale_single_same(lat::Lattice,lenx::T,leny::T,lenz::T,lent::T,ind::Union{Int,Nothing}=nothing) where {T<:SArray{Tuple{lat_size}, Float64,1,lat_size}}

    s = getfactor(4)

    if isnothing(ind)
        index = rand(1:lat_size^4)
    else 
        index = ind
    end

    (i,j,k,l) = lat.node2ind[index]
    
    lenx_new = setindex(lenx, lenx[i]*s[1], i) 
    leny_new = setindex(leny, leny[j]*s[2], j) 
    lenz_new = setindex(lenz, lenz[k]*s[3], k) 
    lent_new = setindex(lent, lent[l]*s[4], l) 

    return (lenx_new, leny_new, lenz_new, lent_new)
end

function next_step_rescale_single_one(lat::Lattice,lenx::Array{Float64,1},leny::Array{Float64,1},lenz::Array{Float64,1},lent::Array{Float64,1},ind::Union{Int,Nothing}=nothing)
    s = getfactor()

    if isnothing(ind)
        index = rand(1:lat_size^4)
    else
    index = ind
    end

    indices = lat.node2ind[index]
    len = [copy(lenx) copy(leny) copy(lenz) copy(lent)]
    d = rand(1:4)
    len[indices[d],d] *= s
    return (len[:,1], len[:,2], len[:,3], len[:,4])
end

function next_step_rescale_n_any(lenx::T, leny::T, lenz::T, lent::T, n::Int=2) where {T<:SArray{Tuple{lat_size}, Float64,1,lat_size}}
    s = getfactor(n)
    lenx_new, leny_new, lenz_new , lent_new = lenx, leny, lenz, lent
    indices = rand(1:4*lat_size, n)
    for i in 1:n
        ind = mod(indices[i],lat_size) + 1
        if div(indices[i],lat_size) == 0
            lenx_new = setindex(lenx_new, lenx_new[ind]*s[i], ind) 
        elseif div(indices[i], lat_size) == 1
            leny_new = setindex(leny_new, leny_new[ind]*s[i], ind) 
        elseif div(indices[i], lat_size) == 2
            lenz_new = setindex(lenz_new, lenz_new[ind]*s[i], ind) 
        else
            lent_new = setindex(lent_new, lent_new[ind]*s[i], ind) 
        end
    end
    return (lenx_new, leny_new, lenz_new, lent_new)
end

function next_step_field_small(vals::Array{Float64,1},d::Float64=del)
    res = similar(vals)
    res =  @. vals + (-d + 2 * d * rand(Size^4))
    res

end
function next_step_field_small_single(vals::Array{Float64,1},ind::Union{Int,Nothing}=nothing,d::Float64=del)

    if isnothing(ind)
        index = rand(1:Size^4)
    else
    index = ind
    end
    #del = 1.0
    new_vals = copy(vals)
    r = rand(1:2)
    if r == 1
        new_vals[index] = vals[index] + (-d + 2 * d* rand())
    else
        new_vals[index] = getsign(new_vals[index])
    end
    new_vals

end

function getsign(x)
    if rand() <= 0.5
        return -x
    else
        return x
    end
end

function  next_step_field_rescale(vals::Array{Float64,1})
    r = rand(1:2)
    if rand == 1
        vals_new = getfactor(Size^4).*vals
    else
        vals_new = getsign.(vals)
    end
    vals_new
end

function next_step_field_small_rescale(vals::Array{Float64,1})
    r = rand(1:2)
    if r == 1
        vals_new = getfactor().*vals
    else
        vals_new = getsign.(vals)
    end
    vals_new
end


function next_step_field_rescale_single(vals::Array{Float64,1},ind::Union{Int,Nothing}=nothing)
    if isnothing(ind)
        index = rand(1:Size^4)
    else
        index = ind
    end
    new_vals = copy(vals)
    r = rand(1:2)
    if r == 1
        new_vals[index] = getfactor()*vals[index]
    else
        new_vals[index] = getsign(vals[index])
    end 
    new_vals
end
