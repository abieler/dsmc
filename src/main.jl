using Octree
using Triangles
using Gas
using Types

include("Physical.jl")
include("io.jl")

lostParticles = Particle[]
myPoint = zeros(Float64,3)

nTriangles, allTriangles = load_ply_file("../input/sphere-Titan-4005.ply")

################################################################################
# initialize simulation domain
################################################################################
origin = zeros(Float64, 3)
halfSize = ones(Float64, 3) * 20.0 #5000.0e3
nCellsX = 5
nCellsY = 5
nCellsZ = 5
isLeaf = true
refLevel = 0
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

oct = Block(origin, halfSize, isLeaf, Array(Block,8), Cell[], refLevel,
            nCellsX, nCellsY, nCellsZ)

#split octree into 8 children
split_block(oct)

################################################################################
# refine domain and compute volume of cut cells
################################################################################
assign_triangles!(oct, allTriangles)
for i=1:4
  refine_tree(oct, 1)
  assign_triangles!(oct, allTriangles)
end
pStart = [5.0, 0.0, 0.0]
@time cut_cell_volume!(oct, pStart, 2000)

################################################################################
# main loop
################################################################################
nParticles = 500
for iteration = 1:4
  println("iteration: ", iteration)
  #insert_new_particles(oct, nParticles, myPoint)
  if iteration == 1
    #insert_new_particles_sphere(oct, nParticles, myPoint,Source,DELTA)
    @time insert_new_particles_body(oct, allTriangles, nParticles, myPoint)
  end
  save_particles(oct, "../output/particles_init_" *string(iteration)* ".csv")
  println("compute mac params")
  compute_macroscopic_params(oct)
  println("mac p done")
  time_step(oct, lostParticles)
  println("timestep done")
  assign_particles!(oct, lostParticles, myPoint)
  println("assign done")
  lostParticles = Particle[]
  println(iteration)
  println()
end

################################################################################
# save results
################################################################################
save2vtk(oct)
