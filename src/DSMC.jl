module DSMC

using Distributions

include("typedefs.jl")
include("triangles.jl")
include("gas.jl")
include("octree.jl")
include("io.jl")
include("raytrace.jl")
include("user.jl")

export SphericalBody,
       Block,
       MeshBody,
       Particle,
       Particles,
       Triangle,
       collect_blocks!,
       compute_macroscopic_params,
       initialize_domain,
       insert_new_particles,
       load_ply_file,
       mySettings,
       refine_domain,
       time_step,
       save2vtk,
       save_particles,
       send_particles_to_cpu,
       lostParticles,
       lostIDs,
       myPoint,
       particle_buffer,
       coords,
       domain,
       sayHi


#globals
const MyID = myid()
const nowpos = zeros(Float64, 3)
const nextpos = zeros(Float64, 3)
const vStartStop = zeros(Float64, 3)
const pI = zeros(Float64, 3)
const u = zeros(Float64, 3)
const v = zeros(Float64, 3)
const w = zeros(Float64, 3)

lostParticles = Particles(200)
lostIDs = Int64[]
myPoint = zeros(Float64, 3)
particle_buffer = Array(Particle, 400)
coords = zeros(Float64, 3)

domain = initialize_domain(mySettings)

end
