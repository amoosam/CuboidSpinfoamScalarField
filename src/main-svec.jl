using LinearAlgebra
using StaticArrays
using Dates
using TOML

cfg = TOML.parsefile("../config.toml")

BLAS.set_num_threads(cfg["machine"]["blasThreads"])

# Simulation Parameters
lat_params = cfg["lattice"]
const lat_size = lat_params["latticeSize"]
const cut = lat_params["upperCutOff"]
const lower_cut = lat_params["lowerCutOff"]
const init = lat_params["upperInitBound"]
const init_lower_cut = lat_params["lowerInitBound"]
const init_field = lat_params["fieldInit"]


include("lattice-svec.jl")

const del_dict = Dict([
  (5500, 0.1), (5600,0.1), (5700, 0.1), (5800, 0.1), (5900,0.1), 
  (6000,0.1), (6100,0.1), (6200,0.1), (6300,0.1), (6400, 0.1), (6500,0.1),
  (6600,0.1), (6700,0.1), (6800,0.1), (6900,0.1), (6900, 0.1), (7000,0.1),
  (7100,0.5), (7200,0.5), (7300,0.2), (7400,0.2), (7500, 0.2), (7600,0.2)
])


folder = today()

const path = string(mkpath(string("../Data/RescaleShift/", replace(string(folder), ['-'] => ""),"_equal")),"/")
println("Data directory: $path")
#const path = string("/home/mali/Data/SF_Matter/RescaleSingle_Update/", replace(string(today()), ['-'] => ""),"/alpha_doublefine/")
const ver = 'a'

include("spin_foam.jl")
include("imp_sampling-svec.jl")
include("plots.jl")
include("proposal-svec.jl")

#setprecision(importance_geom_matter_therm, BigFloat, 128)  # precision to be used for BigFloats within specified function

function main_importance(alpha::Float64, mass::Float64, del::Float64)
    thermalization  = convert(Int, 1e6)
    samples = 100
    steps_between = 10
    lat = Lattice(nodes_to_indices(),indices_to_nodes(), define_adjacency_matrix())
    importance_geom_matter_therm(lat,thermalization, alpha, mass, del)
    println("Thermalization completed.")
    println("=================================")
    println("Starting sampling...")
    importance_geom_matter_sample(lat,samples,steps_between, alpha, mass)
    println("Sampling completed. Data written to disk")
    println("=================================")

end

main_importance(0.6,0.5,0.1)
