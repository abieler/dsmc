using Octree
using Triangles

include("io.jl")

#nTriangles, nodeCoords, triIndices, triangles, n_hat, triCenters, triAreas
nTriangles, allTriangles = load_ply_file("../input/plane.ply")

origin = zeros(Float64, 3)
halfSize = ones(Float64, 3)*2*1000
nCellsX = 5
nCellsY = 5
nCellsZ = 5
isLeaf = 1
refLevel = 0

# initialize octree with no cells.
oct = Block(origin, halfSize, isLeaf, Array(Block,8), Cell[], refLevel, nCellsX, nCellsY, nCellsZ)
#split octree into 8 children
split_block(oct)
for i=1:5
  assign_triangles!(oct, allTriangles)
  refine_tree(oct)
end

save2vtk(oct)
