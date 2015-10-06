using Octree
using Triangles

include("io.jl")

#nTriangles, nodeCoords, triIndices, triangles, n_hat, triCenters, triAreas
nTriangles, allTriangles = load_ply_file("../input/SHAP5_100k.ply")

origin = zeros(Float64, 3)
halfSize = ones(Float64, 3)*5.0
nCellsX = 2
nCellsY = 2
nCellsZ = 2
isLeaf = 1
refLevel = 0

# initialize octree with no cells.
oct = Block(origin, halfSize, isLeaf, Array(Block,8), Cell[], refLevel, nCellsX, nCellsY, nCellsZ)
insert_cells(oct)

#split octree into 8 children
#split_block(oct)

save2vtk(oct)
