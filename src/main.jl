using Octree
using Triangles
using Gas
using Types

include("io.jl")

lostParticles = Particle[]
myPoint = zeros(Float64,3)

nTriangles, allTriangles = load_ply_file("../input/cow.ply")

origin = zeros(Float64, 3)
halfSize = ones(Float64, 3)*2*10
nCellsX = 5
nCellsY = 5
nCellsZ = 5
isLeaf = true 
refLevel = 0

# initialize octree with no cells.
oct = Block(origin, halfSize, isLeaf, Array(Block,8), Cell[], refLevel, nCellsX, nCellsY, nCellsZ)
#split octree into 8 children
split_block(oct)
assign_triangles!(oct, allTriangles)
for i=1:5
  refine_tree(oct)
  assign_triangles!(oct, allTriangles)
end
pStart = [18.0, 0.0, 0.0]
@time cut_cell_volume!(oct, pStart, 2000)

for i=1:1
  insert_new_particles(oct, myPoint)
  compute_macroscopic_params(oct)
  time_step(oct, lostParticles)
  assign_particles!(oct, lostParticles, myPoint)
  lostParticles = Particle[]
  #println(i)
end
save2vtk(oct)
