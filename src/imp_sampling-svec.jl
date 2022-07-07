using Serialization
using HDF5
using DelimitedFiles
using Random
using LinearAlgebra
using StaticArrays

# structure to define and save the state of the simulation
struct State
    lx::SArray{Tuple{lat_size}, Float64,1,lat_size}
    ly::SArray{Tuple{lat_size}, Float64,1,lat_size}
    lz::SArray{Tuple{lat_size}, Float64,1,lat_size}
    lt::SArray{Tuple{lat_size}, Float64,1,lat_size}
    d::Float64
    steps::Int64
    acc::Float32
    rng:: Random.MersenneTwister
end

# functions to load data from file for next thermalization block or sampling. Checks if thermalization file exists and generates new data if it doesn't if called with sample=false, otherwise prints error.
# when version is 'a', thermalization file is overwritten.
function load_lengths(sample::Bool = false,filename::Union{String,Nothing}=nothing)
    if sample == true
        fname = string(path,"thermstate-",filename,ver,".dat")
        if isfile(fname)
            #temp = Ref{state}()
            #open(io -> read!(io,temp),file,"r")
            thst = deserialize(fname)
            copy!(Random.default_rng(),thst.rng)
            return (thst.lx,thst.ly,thst.lz,thst.lt,thst.d,thst.acc)
        else
            println("Error: Thermalization runs not done.")
            return nothing
        end
    else        
        if ver != 'a'
            preVer = Char(Int(ver)-1)
            fname = string(path,"thermstate-",filename,preVer,".dat")
        else 
            fname = string(path,"thermstate-",filename,ver,".dat")
        end
            
        if isfile(fname)
            println("loading lengths from stored thermalization state")
            thst = deserialize(fname)
            copy!(Random.default_rng(),thst.rng)
            lenx, leny, lenz, lent = thst.lx, thst.ly, thst.lz, thst.lt
        else
            println("Thermalization file doesn't exist. Generating new length values")
            lenx = SVector{lat_size}(init_lower_cut .+ (init - init_lower_cut) .* rand(lat_size))
            leny = SVector{lat_size}(init_lower_cut .+ (init - init_lower_cut) .* rand(lat_size))
            lenz = SVector{lat_size}(init_lower_cut .+ (init - init_lower_cut) .* rand(lat_size))
            lent = SVector{lat_size}(init_lower_cut .+ (init - init_lower_cut) .* rand(lat_size))
        end
        return (lenx, leny, lenz, lent)
    end
end

function load_fields(sample::Bool = false,filename::Union{String,Nothing}=nothing)::Union{Array{Float64,1}, Nothing}
    field = Array{Float64,1}
    if sample == true
        fname = string(path,"thermfieldstate-",filename,ver,".txt")
        if isfile(fname)
            field = vec(readdlm(fname, '\t', Float64, '\n'))
            #open(io -> read!(io,field),file,"r")
            return field
        else
            println("Thermalization runs not done. Fields not found.")
            return field
        end
    else
        if ver != 'a'
            preVer = Char(Int(ver)-1)
            fname = string(path,"thermfieldstate-",filename,preVer,".txt")
        else
            fname = string(path,"thermfieldstate-",filename,ver,".txt")
        end

        if isfile(fname)
            println("loading fields from stored thermalization state")
            field = vec(readdlm(fname))
            #open(io -> read!(io,field),fname,"r")
            return field
        else
            println("Thermalization file doesn't exist. Generating new field values")
            field = -init_field .+ (2 * init_field) .* rand(lat_size^4)
            return field
        end
    end
end

#calculates the full amplitude with lengths and scalar fields
function full_amplitude(lat::Lattice, len_x::T,len_y::T,len_z::T,len_t::T,
    dual_x::T,dual_y::T,dual_z::T,dual_t::T,
    field_vals::Array{Float64,1}, alpha::Float64, mass::Float64) where {T<:SArray{Tuple{lat_size}, Float64,1,lat_size}}
    
    Sact  = BigFloat(0.)
    Sfamp = BigFloat(1.)

    Mat = define_field_matrix(lat,mass,len_x,len_y,len_z,len_t, dual_x,dual_y,dual_z,dual_t)

    Sfamp *= prod([Am(alpha,len_x[v[1]],len_y[v[2]],len_z[v[3]],len_t[v[4]]) for v in values(lat.node2ind)])
    
    Sact += (Transpose(field_vals) * Mat * field_vals)
    res::BigFloat = Sfamp * exp(-Sact/2.0)

    return res, Mat

end

function partial_amplitude_fields(lat::Lattice, Mat::Array{Float64,2}, field_vals_old::Array{Float64,1}, field_vals_new::Array{Float64,1}, ind::Int)
    deltaS = BigFloat(0.)
    delta = field_vals_new[ind] - field_vals_old[ind]
    neighbors = lat.adjacency[ind,:]
    deltaS += delta^2 * (Mat[ind,ind]) + delta * (2.0 * Mat[ind,ind] * field_vals_old[ind] + 2.0 * sum([field_vals_old[j]*Mat[ind,j] for j in neighbors]))

    res::BigFloat = exp(-deltaS/2.0)
    return res
end

# function for thermalization
function importance_geom_matter_therm(lat::Lattice, thermalization::Int, alpha::Float64, mass::Float64, del::Float64)

    filename = string(convert(Int,round(alpha*1e4)), "_", convert(Int,round(mass*1e4)))
    f = x -> (x > cut) || (x < lower_cut)
    steps = convert(Int,round(1e5))
    counter = 1
    d = del
    
    amp_therm = zeros(BigFloat, convert(Int,thermalization / steps))
    vol_therm = zeros(BigFloat,convert(Int,thermalization / steps))
    
    acc = zeros(Float64,convert(Int,thermalization / steps))
    acc_therm = 0.0

    #Pick random initial data  or load from previous thermalization run
    lengths = load_lengths(false,filename)
    duals = (duallength(lengths[1]), duallength(lengths[2]),duallength(lengths[3]),duallength(lengths[4]))

    field_vals = load_fields(false, filename)
    amp, Mat = full_amplitude(lat,lengths..., duals..., field_vals, alpha, mass)

    println("Initial amplitude = ", amp)
    
    for i in 1:thermalization
        ind = rand(1:lat_size^4)
        choice = rand(1:2)
        if choice == 1
            # multiple rescale functions available in proposal-svec.jl
            lengths_new =  next_step_rescale_single_same(lat,lengths..., ind)
            if sum(count.(f,lengths_new)) == 0
                duals_new = (duallength(lengths_new[1]), duallength(lengths_new[2]),duallength(lengths_new[3]),duallength(lengths_new[4]))
                new_amp, new_Mat = full_amplitude(lat,lengths_new..., duals_new..., field_vals, alpha, mass)
                if rand() < real(new_amp / amp)
                    lengths = lengths_new
                    duals = duals_new
                    amp = new_amp
                    Mat = new_Mat
                    @inbounds acc[counter] += 1
                end
            end
        else
            field_vals_new = next_step_field_small_single(field_vals,ind,d)
            delta_amp = partial_amplitude_fields(lat,Mat, field_vals, field_vals_new, ind)
            if rand() < delta_amp 
                field_vals = field_vals_new
                amp *= delta_amp
                @inbounds acc[counter] += 1
            end

        end   

        if mod(i,steps) == 0
            @inbounds amp_therm[counter] = amp
            @inbounds vol_therm[counter] = sum(convert.(BigFloat,lengths[1])) * sum(convert.(BigFloat,lengths[2])) *
                                    sum(convert.(BigFloat,lengths[3])) * sum(convert.(BigFloat,lengths[4]))

            acc_therm = acc[counter] / steps
            if acc_therm < 0.3 
                d *= 0.75
                println("Shift parameter reduced by 25%. shift =$d")
            elseif acc_therm > 0.6
                d *= 1.25
                println("Shift parameter increased by 25%. shift = $d")
            end
            counter += 1
        end
    end
    acc_therm = sum(acc) / thermalization
    state_chk = State(lengths..., d, thermalization, acc_therm, copy(Random.default_rng()))
    acc = acc ./ steps


    fname = open(string(path,"thermstate-",filename,ver,".dat"),"w")
    serialize(fname, state_chk)
    close(fname)
    open(io -> writedlm(io,field_vals),string(path,"thermfieldstate-",filename,ver,".txt"),"w")

    println("Thermalization acceptance rate after $thermalization steps = $acc_therm")
    println("Amplitude at $thermalization step = $amp")
    println("=================================")
    plot_therm_acc_results(alpha, mass, acc,amp_therm, vol_therm,filename)
end

function importance_geom_matter_sample(lat::Lattice, samples::Int, steps_between::Int, alpha::Float64, mass::Float64)
    filename = string(convert(Int,round(alpha*1e4)), "_", convert(Int,round(mass*1e4)))
    f = x -> (x > cut) || (x < lower_cut)
    #observables
    amplitude = zeros(BigFloat,samples)
    acceptance_rate = zeros(Float64, samples)
    sample_data = zeros(Float64,(lat_size,samples,4))
    field_samples = zeros(Float64, lat_size^4,samples)
    acc = 0

    len_x,len_y,len_z,len_t, shift, acc_therm = load_lengths(true,filename)

    if acc_therm == 0
        println("Thermalization runs unsuccessful. Acceptance rate = $acc_therm. Tweak initial configuration.")
        return nothing
    end

    lengths = (len_x,len_y,len_z,len_t)
    field_vals = load_fields(true,filename)

    duals =  (duallength(lengths[1]), duallength(lengths[2]),duallength(lengths[3]),duallength(lengths[4]))
    amp, Mat = full_amplitude(lat,lengths..., duals..., field_vals,alpha, mass)

    for i in 1:samples
        @inbounds for _ in 1:steps_between
            choice = rand(1:2)
            ind = rand(1:lat_size^4)
            if choice == 1
                lengths_new =  next_step_rescale_single_same(lat,lengths..., ind)
                #lengths_new = next_step_rescale_n_any(lengths...,2)
                if sum(count.(f,lengths_new)) == 0
                    duals_new = (duallength(lengths_new[1]), duallength(lengths_new[2]),duallength(lengths_new[3]),duallength(lengths_new[4]))
                    new_amp, new_Mat = full_amplitude(lat,lengths_new..., duals_new..., field_vals, alpha, mass)
                    if rand() < real(new_amp / amp)
                        lengths = lengths_new
                        duals = duals_new
                        amp = new_amp
                        Mat = new_Mat
                        acc += 1
                    end
                end
            else
                field_vals_new = next_step_field_small_single(field_vals,ind,shift)
                delta_amp = partial_amplitude_fields(lat, Mat, field_vals, field_vals_new, ind)
                if rand() < delta_amp 
                    field_vals = field_vals_new
                    amp *= delta_amp
                    acc += 1
                end
    
            end
        end

        @inbounds acceptance_rate[i] = acc /steps_between
        @inbounds amplitude[i] =  amp
        @inbounds sample_data[:,i,1] = lengths[1]
        @inbounds sample_data[:,i,2] = lengths[2]
        @inbounds sample_data[:,i,3] = lengths[3]
        @inbounds sample_data[:,i,4] = lengths[4]
        @inbounds field_samples[:,i] = field_vals
        acc = 0
    end
    println("acceptance rate during sampling: ", sum(acceptance_rate)/samples)
    fname = string(path,"samples-",filename,ver,".h5")
    fid = h5open(fname,"w")
    fid["lengths"] = convert.(Float64,sample_data)
    fid["fields"] = field_samples
    close(fid)
    plot_results(alpha, mass, acceptance_rate,amplitude, filename)

end
