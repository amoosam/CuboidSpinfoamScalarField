#using StaticArrays
#using LoopVectorization

struct Lattice
    node2ind::Dict{Int, NTuple{4,Int}}
    ind2node::Dict{NTuple{4,Int},Int}
    adjacency::Array{Int,2}
end

function per_bdry(i::Int)

    res = i

    if i == 0

        res = lat_size

    elseif i == lat_size + 1

        res = 1

    end

    res

end


function comb_distance(i::Int,j::Int)

    count1 = 0
    count2 = 0

    i_temp1 = i
    i_temp2 = i

    while i_temp1 != j

        i_temp1 = per_bdry(i_temp1+1)
        count1 += 1

    end

    while i_temp2 != j

        i_temp2 = per_bdry(i_temp2-1)
        count2 += 1

    end

    return min(count1,count2)

end

function phys_distance(i::Int, j::Int, len::SArray{Tuple{lat_size},Float64,1,lat_size})
    dist1 = zero(Float64)
    dist2 = zero(Float64)
    temp = i

    while temp != j
        dist1 += len[temp]
        temp = per_bdry(temp+1)
    end

    temp = i
    while temp!=j
        temp = per_bdry(temp-1)
        dist2 += len[temp]
    end

    return min(dist1, dist2)
end

function nodes_to_indices()
    counter::Int = 1
    d = Dict{Int,NTuple{4,Int}}()
    for i in 1:lat_size
        for j in 1:lat_size
            for k in 1:lat_size
                for l in 1:lat_size
                    d[counter] = (i,j,k,l)
                    counter += 1
                end
            end
        end
    end
    return d             
end

function indices_to_nodes()
    counter::Int = 1
    d = Dict{NTuple{4,Int},Int}()
    for i in 1:lat_size
        for j in 1:lat_size
            for k in 1:lat_size
                for l in 1:lat_size
                    d[(i,j,k,l)] = counter
                    counter += 1
                end
            end
        end
    end
    return d
end

function duallength(length::SArray{Tuple{lat_size},Float64,1,lat_size})
    dual = @SVector [1.0/2.0 * (length[i] + length[per_bdry(i-1)]) for i in 1:lat_size]
    return dual
end

function dual3vol(len1::Real, len2::Real, len3::Real)
    res = 0.
    res = len1 * len2 * len3
    return res
end


function dual4vol(len1::Real, len2::Real, len3::Real, len4::Real)
    res = 0.
    res =  len1 * len2 * len3 * len4 
    return res
end

function define_adjacency_matrix()
    res = zeros(Int,lat_size^4,8)
    for ind in 1:lat_size^4
        (i,j,k,l) = all_indices(ind)
        res[ind,1] = index_length(i,j,k,per_bdry(l+1))
        res[ind,2] = index_length(i,j,k,per_bdry(l-1))
        res[ind,3] = index_length(i,j,per_bdry(k+1),l)
        res[ind,4] = index_length(i,j,per_bdry(k-1),l)
        res[ind,5] = index_length(i,per_bdry(j+1),k,l)
        res[ind,6] = index_length(i,per_bdry(j-1),k,l)
        res[ind,7] = index_length(per_bdry(i+1),j,k,l)
        res[ind,8] = index_length(per_bdry(i-1),j,k,l)
    end
    res
end

function define_field_matrix(lat::Lattice,mass::Float64,len_x::T,len_y::T,len_z::T, len_t::T, dual_x::T, dual_y::T, dual_z::T, dual_t::T) where {T<:SArray{Tuple{lat_size},Float64,1, lat_size}}

    res = zeros((lat_size^4, lat_size^4))
@inbounds Threads.@threads for ind in 1:lat_size^4
        (i,j,k,l) = lat.node2ind[ind]
        res[ind,ind] = (dual3vol(dual_y[j],dual_z[k],dual_t[l])*(1. / len_x[i] + 1. / len_x[per_bdry(i-1)]) +
                         dual3vol(dual_x[i],dual_z[k],dual_t[l])*(1. / len_y[j] + 1. / len_y[per_bdry(j-1)]) +
                         dual3vol(dual_x[i],dual_y[j],dual_t[l])*(1. / len_z[k] + 1. / len_z[per_bdry(k-1)]) +
                         dual3vol(dual_x[i],dual_y[j],dual_z[k])*(1. / len_t[l] + 1. / len_t[per_bdry(l-1)]) +
                         mass^2 * dual4vol(dual_x[i],dual_y[j], dual_z[k], dual_t[l]))
        
        res[ind,lat.ind2node[(i,j,k,per_bdry(l-1))]] = -dual3vol(dual_x[i],dual_y[j],dual_z[k]) / len_t[per_bdry(l-1)]
        res[ind,lat.ind2node[(i,j,k,per_bdry(l+1))]] = -dual3vol(dual_x[i],dual_y[j],dual_z[k]) / len_t[l]
        res[ind,lat.ind2node[(i,j,per_bdry(k-1),l)]] = -dual3vol(dual_x[i],dual_y[j],dual_t[l]) / len_z[per_bdry(k-1)]
        res[ind,lat.ind2node[(i,j,per_bdry(k+1),l)]] = -dual3vol(dual_x[i],dual_y[j],dual_t[l]) / len_z[k]
        res[ind,lat.ind2node[(i,per_bdry(j-1),k,l)]] = -dual3vol(dual_x[i],dual_z[k],dual_t[l]) / len_y[per_bdry(j-1)]
        res[ind,lat.ind2node[(i,per_bdry(j+1),k,l)]] = -dual3vol(dual_x[i],dual_z[k],dual_t[l]) / len_y[j]
        res[ind,lat.ind2node[(per_bdry(i-1),j,k,l)]] = -dual3vol(dual_y[j],dual_z[k],dual_t[l]) / len_x[per_bdry(i-1)]
        res[ind,lat.ind2node[(per_bdry(i+1),j,k,l)]] = -dual3vol(dual_y[j],dual_z[k],dual_t[l]) / len_x[i]
    
    end
    return res
end

function index_length(i::Int,j::Int,k::Int,l::Int)

    res = lat_size^3 * (i-1) + lat_size^2 * (j-1) + lat_size * (k-1) + l

    res

end


function all_indices(index::Int)

    i = 0
    j = 0
    k = 0
    l = 0

    temp = index

    while temp > 0

        i += 1

        temp -= lat_size^3

    end

    temp += lat_size^3

    while temp > 0

        j += 1

        temp -= lat_size^2

    end

    temp += lat_size^2

    while temp > 0

        k += 1

        temp -= lat_size

    end

    temp += lat_size

    l = temp

    return i,j,k,l

end

function define_corr_length(i::Int,j::Int,k::Int,l::Int,inver_mat::Array)

    res = zeros(lat_size^4,2)

    ref_index = index_length(i,j,k,l)

    for index2 in 1:lat_size^4

        (i2,j2,k2,l2) = all_indices(index2)

        #res[index2,1] = sqrt(comb_distance(i,i2)^2 + comb_distance(j,j2)^2 + comb_distance(k,k2)^2 +
        #                comb_distance(l,l2)^2)

        res[index2,1] = comb_distance(i,i2) + comb_distance(j,j2) + comb_distance(k,k2) +
                        comb_distance(l,l2)

        res[index2,2] = inver_mat[ref_index,index2]

    end

    return res

end

function define_phys_corr_length(i::Int,j::Int,k::Int,l::Int,lenx::Array{Float64,1}, leny::Array{Float64,1}, lenz::Array{Float64,1}, lent::Array{Float64,1}, inver_mat::Array)
    res = zeros(lat_size^4,2)

    ref_index = index_length(i,j,k,l)

    for index2 in 1:lat_size^4

        (i2,j2,k2,l2) = all_indices(index2)

        #res[index2,1] = sqrt(comb_distance(i,i2)^2 + comb_distance(j,j2)^2 + comb_distance(k,k2)^2 +
        #                comb_distance(l,l2)^2)

        res[index2,1] = phys_distance(i,i2,lenx) + phys_distance(j,j2,leny) + phys_distance(k,k2,lenz) +
                        phys_distance(l,l2,lent)

        res[index2,2] = inver_mat[ref_index,index2]

    end

    return res
end