function DE(j1, j2, j3, j4, j5, j6)
    res = (-2) * (j1^2 * (j2 + j4) + j2 * j4 * (j2 + j4) + j1 * (j2^2 + (1. + 1im) * j2 * j4 + j4^2))
    res *= (j1^2 * (j3 + j5) + j3 * j5 * (j3 + j5) + j1 * (j3^2 + (1. + 1im) * j3 * j5 + j5^2))
    res *= (j3 * j4 * j5 + j2 * (j4 * j5 + j3 * (j4 + j5)))
    res *= (j2^2 * (j3 + j6) + j3 * j6 * (j3 + j6) + j2 * (j3^2 + (1. + 1im) * j3 * j6 + j6^2))
    res *= (j4^2 * (j5 + j6) + j5 * j6 * (j5 + j6) + j4 * (j5^2 + (1. + 1im) * j5 * j6 + j6^2))
    res *= (j3 * j4 * j6 + j1 * (j4 * j6 + j3 * (j4 + j6)))
    res *= (j2 * j5 * j6 + j1 * (j5 * j6 + j2 * (j5 + j6)))

    res
end


function DF(j1,j2,j3,j4,j5,j6)
    res = (-2) * (j1^2 * (j2 + j4) + j2 * j4 * (j2 + j4) + j1 * (j2^2 + (1. - 1im) * j2 * j4 + j4^2))
    res *= (j1^2 * (j3 + j5) + j3 * j5 * (j3 + j5) + j1 * (j3^2 + (1. - 1im) * j3 * j5 + j5^2))
    res *= (j3 * j4 * j5 + j2 * (j4 * j5 + j3 * (j4 + j5)))
    res *= (j2^2 * (j3 + j6) + j3 * j6 * (j3 + j6) + j2 * (j3^2 + (1. - 1im) * j3 * j6 + j6^2))
    res *= (j4^2 * (j5 + j6) + j5 * j6 * (j5 + j6) + j4 * (j5^2 + (1. - 1im) * j5 * j6 + j6^2))
    res *= (j3 * j4 * j6 + j1 * (j4 * j6 + j3 * (j4 + j6)))
    res *= (j2 * j5 * j6 + j1 * (j5 * j6 + j2 * (j5 + j6)))

    res
end


function Ampl(alpha,j1,j2,j3,j4,j5,j6)

    res = 1/((2. * pi)^3) * 2. ^(24)
    res *= (1/(sqrt(-DE(j1,j2,j3,j4,j5,j6))) + 1/(sqrt(-DF(j1,j2,j3,j4,j5,j6))))
    res *= (1/(sqrt(-DE(j1,j2,j3,j4,j5,j6))) + 1/(sqrt(-DF(j1,j2,j3,j4,j5,j6))))
    res *= (j1 + j2) * (j1 + j3) * (j1 + j4) * (j1 + j5) * (j2 + j3) * (j2 + j4)
    res *= (j2 + j6) * (j3 + j5) * (j3 + j6) * (j4 + j5) * (j4 + j6) * (j5 + j6)
    res *= (j1 * j2 * j3 * j4 * j5 * j6)^(2 * alpha)

    res
end


function FP(x, y, z, t)

    result = 2.0 * sqrt(x^2 * y^2 * z^2 * (x^2 + y^2 + z^2) +
    t^4 * (y^2 * z^2 + x^2 * (y^2 + z^2)) +
    t^2 * (x^4 * (y^2 + z^2) + y^2 * z^2 * (y^2 + z^2) + x^2 * (y^4 + z^4)))

    result

end


function Am(a,x,y,z,t)

    result = Ampl(a, y*z, x*y, x*z, y*t, z*t, x*t)
    result *= FP(x, y, z, t)

    convert(BigFloat,real(result))

end


# function full_sf_amplitude(a,lenx::Array,leny::Array,
#         lenz::Array,lent::Array)
#     res = 1.

#     for i in 1:lat_size, j in 1:lat_size, k in 1:lat_size, l in 1:lat_size

#         res *= Am(a,lenx[i],leny[j],lenz[k],lent[l]) #/ Am(a,1.,1.,1.,1.)

#     end

#     real(res)

# end

function full_sf_amplitude(a,lat::Lattice, lenx::T, leny::T, lenz::T, lent::T ) where {T<:SArray{Tuple{lat_size}, Float64,1,lat_size}}

    res = BigFloat(1.)
    @inbounds Threads.@threads for node in 1:lat_size^4
        (i,j,k,l) = lat.node2ind[node]
        res*= Am(a,lenx[i],leny[j],lenz[k],lent[l]) 
    end
    return res
end