using Octree
using Triangles
using Gas
using Types

include("io.jl")

lostParticles = Particle[]
myPoint = zeros(Float64,3)

nTriangles, allTriangles = load_ply_file("../input/cow.ply")

################################################################################
# initialize simulation domain
################################################################################
origin = zeros(Float64, 3)
halfSize = ones(Float64, 3) * 20
nCellsX = 5
nCellsY = 5
nCellsZ = 5
isLeaf = true
refLevel = 0

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
  #insert_new_particles(oct, nParticles, myPoint)
  if iteration == 1
    @time insert_new_particles_body(oct, allTriangles, nParticles, myPoint)
  end
  save_particles(oct, "../output/particles_init_" *string(iteration)* ".csv")
  compute_macroscopic_params(oct)
  @time time_step(oct, lostParticles)
  @time assign_particles!(oct, lostParticles, myPoint)
  lostParticles = Particle[]
  println(iteration)
  println()
end

################################################################################
# save results
################################################################################
save2vtk(oct)
