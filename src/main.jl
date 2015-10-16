using Octree
using Triangles
using Gas
using Types

include("io.jl")

lostParticles = Particle[]
myPoint = zeros(Float64,3)

nTriangles, allTriangles, surfaceArea = load_ply_file("../input/cow.ply")

################################################################################
# initialize simulation domain
################################################################################
origin = zeros(Float64, 3)
halfSize = ones(Float64, 3) * 10
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
for i=1:3
  refine_tree(oct, 1)
  assign_triangles!(oct, allTriangles)
end
pStart = [150.0, 0.0, 0.0]
@time cut_cell_volume!(oct, pStart, 2000)

################################################################################
# give every cell a unique ID
################################################################################
label_cells!(oct)

################################################################################
# main loop
################################################################################
const nParticles = 20000
const f = nParticles / surfaceArea

println("f: ", f)
println("surfaceArea: ", surfaceArea)

for iteration = 1:200
  if iteration == 1
    #insert_new_particles_sphere(oct, nParticles, myPoint)
    insert_new_particles_body(oct, allTriangles, f, myPoint)
  end
  if iteration % 20 == 1
    save_particles(oct, "../output/particles_init_" *string(iteration)* ".csv")
  end
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
#save2vtk(oct)
