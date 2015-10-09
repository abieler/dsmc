using Octree
using Triangles
using Gas

include("io.jl")

lostParticles = Particle[]
myPoint = zeros(Float64,3)

nTriangles, allTriangles = load_ply_file("../input/cow.ply")

origin = zeros(Float64, 3)
halfSize = ones(Float64, 3)*2*10
nCellsX = 5
nCellsY = 5
nCellsZ = 5
isLeaf = 1
refLevel = 0

# initialize octree with no cells.
oct = Block(origin, halfSize, isLeaf, Array(Block,8), Cell[], refLevel, nCellsX, nCellsY, nCellsZ)
#split octree into 8 children
split_block(oct)
for i=1:6
  assign_triangles!(oct, allTriangles)
  refine_tree(oct)
end

for i=1:1000
  insert_new_particles(oct, myPoint)
  compute_macroscopic_params(oct)
  @time time_step(oct, lostParticles)
  @time assign_particles!(oct, lostParticles, myPoint)
  lostParticles = Particle[]
  println(i)
end
save2vtk(oct)
