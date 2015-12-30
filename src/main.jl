using Triangles
using Types
@everywhere include("gas.jl")
@everywhere include("octree.jl")
@everywhere include("Physical.jl")
@everywhere include("io.jl")
@everywhere include("user.jl")

@everywhere const global MyID = myid()

@everywhere lostParticles = Particles(200)
@everywhere lostIDs = Int64[]
@everywhere myPoint = zeros(Float64, 3)
@everywhere particle_buffer = Array(Particle, 400)

################################################################################
# load shape model
################################################################################
@everywhere nTriangles, allTriangles, surfaceArea = load_ply_file(mySettings.meshFileName)
@everywhere flux = mySettings.nNewParticlesPerIteration / surfaceArea
@everywhere body = MeshBody(allTriangles, 300.0, 2e24, Float64[flux],
                            Float64[18.0*amu],  Float64[1.e5])
################################################################################
# initialize simulation domain
################################################################################
@everywhere oct = initialize_domain(mySettings)

################################################################################
# refine domain and compute volume of cut cells
################################################################################
@everywhere refine_domain(oct, allTriangles, mySettings)
@everywhere allBlocks = Block[]
@everywhere collect_blocks!(oct, allBlocks)

################################################################################
# main loop
################################################################################
@everywhere const nParticles = mySettings.nNewParticlesPerIteration
@everywhere const f = nParticles / surfaceArea
@time begin
  for iteration = 1:mySettings.nIterations
    @everywhere insert_new_particles(oct, body, myPoint)
    @everywhere compute_macroscopic_params(oct)
    @everywhere time_step(oct, lostParticles)
    @everywhere send_particles_to_cpu(lostParticles)
    @everywhere lostParticles.nParticles = 0
    @show(iteration)
    readline(STDIN)
  end
end

################################################################################
# save results
################################################################################
#remotecall_fetch(1, save2vtk, oct)
for iProc in workers()
  remotecall_fetch(iProc, save_particles, "../output/particles_" *string(iProc)* ".csv")
end
