using Dates
using StaticArrays
using LinearAlgebra
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

sim_params = cfg["MCsimulation"]
const thermalization  = convert(Int,sim_params["thermalizationSteps"])
const samples = sim_params["samples"]
const steps_between = sim_params["samplesToSkip"]

const ind = parse(Int64,ARGS[1])
const alpha = 0.61 + (ind * 0.001)
const mass = parse(Float64,ARGS[2])

const del_dict = Dict([
  (5500, 0.1), (5600,0.1), (5700, 0.1), (5800, 0.2), (5900,0.2), 
  (6000,3.0), (6010,3.0), (6020,3.0), (6030,3.0), (6040,3.0), (6050,3.0), (6060,3.0), (6070,3.0), (6080,3.0), (6090,3.0),
  (6100,3.0), (6110,3.0), (6120,3.0), (6130,3.0), (6140,3.0), (6150,3.0), (6160,3.0), (6170,3.0), (6180,3.0), (6190,3.0),
  (6200,3.0), (6210,3.0), (6220,3.0), (6230,3.0), (6240,3.0), (6250,3.0), (6260,3.0), (6270,3.0), (6280,3.0), (6290,3.0),
  (6300,3.0), (6310,3.0), (6320,3.0), (6330,3.0), (6340,3.0), (6350,3.0), (6360,3.0), (6370,3.0), (6380,3.0), (6390,3.0),
  (6400, 3.0), (6410, 3.0), (6420, 3.0), (6430, 3.0), (6440, 3.0), (6450, 3.0), (6460, 3.0), (6470, 3.0), (6480, 3.0), (6490, 3.0),
  (6500,2.0), (6510,2.0), (6520,2.0), (6530,2.0), (6540,2.0), (6550,2.0), (6560,2.0), (6570,2.0), (6580,2.0), (6590,2.0),
  (6600,1.2), (6610,1.2), (6620,1.2), (6630,1.2), (6640,1.2), (6650,1.2), (6660,1.2), (6670,1.2), (6680,1.2), (6690,1.2), 
  (6700,1.2), (6710,1.2), (6720,1.2), (6730,1.2), (6740,1.2), (6750,1.2), (6760,1.2), (6770,1.2) ,(6780,1.2), (6790,1.2),
  (6800,0.5), (6810,0.5), (6820,0.5), (6830,0.5), (6840,0.5), (6850,0.5), (6860,0.5), (6870,0.5), (6880,0.5), (6890,0.5),
  (6900,0.5), (6910,0.5), (6920,0.5), (6930,0.5), (6940,0.5), (6950,0.5), (6960,0.5), (6970,0.5), (6980,0.5), (6990,0.5),
  (7000,0.5), (7010,0.5), (7020,0.5), (7030,0.5), (7040,0.5), (7050,0.5), (7060,0.5), (7070,0.5), (7080,0.5), (7090,0.5),
  (7100,0.5), (7200,0.5), (7300,0.2), (7400,0.2), (7500, 0.2), (7600,0.2), (7700,0.2), (7800,0.2), (7900,0.2), (8000,0.2)
])

const del = del_dict[convert(Int64, round(alpha * 1e4))] 

#folder = today()
folder = Date(2022,04,02)
const path = string(mkpath(string(cfg["data"]["outputDir"], replace(string(folder), ['-'] => ""),"_equal")),"/")
println("Data directory: $path")
const ver = 'a'

function write_simulation_details()
    fname = string(path,"array_simulation_details.txt")
    open(fname, "w") do io
        write(io, "size = $lat_size \n")
        write(io, "thermalization steps = $thermalization \n")
        write(io, "# of samples = $samples \n")
        write(io, "steps_between = $steps_between \n")
        write(io, "upper_cutoff = $cut \n")
        write(io, "lower_cutoff = $lower_cut \n")
        write(io, "init_lower_cut = $init_lower_cut \n")
        write(io, "init_field = $init_field \n")
        write(io, "proposal scheme: either length or field with probability 1/2. when length, change all lengths at a node, when field rescale field or flip sign.\n")
  end
end


include("lattice-svec.jl")
include("spin_foam.jl")
include("imp_sampling-svec.jl")
include("plots.jl")
include("proposal-svec.jl")

#setprecision(importance_geom_matter_therm, BigFloat, 128)  # precision to be used for BigFloats within specified function

function importance_sampling()

    println("alpha:$alpha , mass:$mass")
    lat = Lattice(nodes_to_indices(),indices_to_nodes(), define_adjacency_matrix())
    importance_geom_matter_therm(lat, thermalization, alpha, mass, del)
    println("Thermalization completed.")
    println("=================================")
    println("Starting sampling...")
    importance_geom_matter_sample(lat, samples, steps_between, alpha, mass)
    println("Sampling completed. Data written to disk")
    println("=================================")

end

write_simulation_details()
importance_sampling()
