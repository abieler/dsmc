using Triangles
using Types
@everywhere include("gas.jl")
@everywhere include("octree.jl")
@everywhere include("Physical.jl")
@everywhere include("io.jl")
@everywhere include("user.jl")

@everywhere const global MyID = myid()

@everywhere lostParticles = Particle[]
@everywhere lostIDs = Int64[]
@everywhere myPoint = zeros(Float64, 3)

@everywhere particle_buffer = Array(Particle, 400)

################################################################################
# initialize source, time step, particle, weight, number of rep particles???
################################################################################
#=
bodyradius = TitanRadius+1430.0*km
bodymass = TitanMass
density = 4.48e13
temperature = 161.0
Source = UniformSource(bodyradius,bodymass,density,temperature)
DELTA=time_step(Source.SourceTemperature,28.0*amu)
println(" - SourceTemperature: ", Source.SourceTemperature)
println(" - TimeStep: ", DELTA)
w_factor = constant_weight(DELTA,Source,28.0*amu)
println(w_factor)
=#

################################################################################
# load shape model
################################################################################
@everywhere nTriangles, allTriangles, surfaceArea = load_ply_file(mySettings.meshFileName)
@everywhere flux = mySettings.nNewParticlesPerIteration / surfaceArea
@everywhere body = MeshBody(allTriangles, 300.0, Float64[flux], Float64[18.0*amu],
                Float64[1.e5] )
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
#@everywhere blocks2proc!(allBlocks)
#@everywhere rr = distribute(allBlocks)

################################################################################
# main loop
################################################################################
@everywhere const nParticles = mySettings.nNewParticlesPerIteration
@everywhere const f = nParticles / surfaceArea

for iteration = 1:mySettings.nIterations
  tic()
  println("iteration: ", iteration)
  if iteration >= 1
    @everywhere insert_new_particles(oct, body, myPoint)
  end
  @everywhere compute_macroscopic_params(oct)
  @everywhere time_step(oct, lostParticles, particle_buffer)
  @everywhere send_lost_particles(lostParticles, lostIDs)
  toc()
end

################################################################################
# save results
################################################################################
remotecall_fetch(1, save2vtk, oct)
for iProc in workers()
  remotecall_fetch(iProc, save_particles, "../output/particles_" *string(iProc)* ".csv")
end
