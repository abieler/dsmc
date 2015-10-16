using Octree
using Triangles
using Gas
using Types

include("Physical.jl")
include("io.jl")
include("user.jl")

lostParticles = Particle[]
myPoint = zeros(Float64,3)

################################################################################
# initialize source, time step, particle, weight, number of rep particles???
################################################################################
bodyradius = TitanRadius+1430.0*km
bodymass = TitanMass
density = 4.48e13
temperature = 161.0
Source = UniformSource(bodyradius,bodymass,density,temperature)
DELTA=time_step(Source.SourceTemperature,28.0*amu)
println("SourceTemperature: ", Source.SourceTemperature)
println("TimeStep: ", DELTA)
w_factor = constant_weight(DELTA,Source,28.0*amu)
println(w_factor)

################################################################################
# load shape model
################################################################################
nTriangles, allTriangles, surfaceArea = load_ply_file(mySettings.meshFileName)

################################################################################
# initialize simulation domain
################################################################################
oct = initialize_domain(mySettings)

################################################################################
# refine domain and compute volume of cut cells
################################################################################
refine_domain(oct, allTriangles, mySettings)

################################################################################
# main loop
################################################################################
const nParticles = mySettings.nNewParticlesPerIteration
const f = nParticles / surfaceArea

for iteration = 1:mySettings.nIterations
  println("iteration: ", iteration)
  if iteration >= 1
    #insert_new_particles(oct, nParticles, myPoint, Source, DELTA)
    insert_new_particles(oct, allTriangles, f, myPoint, Source)
  end
  if iteration % 20 == 1
    save_particles(oct, "../output/particles_" *string(iteration)* ".csv")
  end
  compute_macroscopic_params(oct)
  @time time_step(oct, lostParticles)
  @time assign_particles!(oct, lostParticles, myPoint)
  lostParticles = Particle[]
end

################################################################################
# save results
################################################################################
#save2vtk(oct)
