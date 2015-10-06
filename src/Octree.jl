module Octree
using Triangles
using PyCall

@pyimport abieler.shapeUtils as sU

export Cell,
       Block,
       split_block,
       insert_cells,
       blockContainingPoint,
       cellContainingPoint,
       populate_blocks,
       out_of_bounds,
       count_cells,
       refine_tree,
       octree_slice!,
       assign_triangles!,
       save2vtk

type Cell
  origin::Vector{Float64}
  halfSize::Vector{Float64}
  nodes::Array{Float64,2}
  volume::Float64
  data::Vector{Float64}
  triangles::Vector{Triangle}
  hasTriangles::Bool
end

type Block
  origin::Vector{Float64}
  halfSize::Vector{Float64}
  isLeaf::Int64
  children::Vector{Block}
  cells::Vector{Cell}
  nestingLevel::Int64
  nx::Int64
  ny::Int64
  nz::Int64
end

type Point3D
  x::Float64
  y::Float64
  z::Float64
end

function save2vtk(oct)

  indexTransform = Dict{Int64, Int64}()
  indexTransform[1] = 1
  indexTransform[2] = 5
  indexTransform[3] = 3
  indexTransform[4] = 7
  indexTransform[5] = 2
  indexTransform[6] = 6
  indexTransform[7] = 4
  indexTransform[8] = 8

  allCells = Cell[]
  all_cells!(oct, allCells)
  nCells = length(allCells)
  println("nCells: ", nCells)
  epsilon = 1e-10
  coord = zeros(Float64,3)

  #fill in first coordinate. --> needs to be a non empty list to later be
  # able to loop over
  nodeCoords = Point3D[]
  cell = allCells[1]
  for i=1:3
    coord[i] = cell.nodes[i,1]
  end
  push!(nodeCoords, Point3D(coord[1], coord[2], coord[3]))

  for cell in allCells
    for nNode = 1:8
      for i=1:3
        coord[i] = cell.nodes[i,nNode]
      end

      sameScore = 0
      for i=1:length(nodeCoords)
        cond1 = (coord[1]-epsilon <= nodeCoords[i].x <= coord[1]+epsilon)
        cond2 = (coord[2]-epsilon <= nodeCoords[i].y <= coord[2]+epsilon)
        cond3 = (coord[3]-epsilon <= nodeCoords[i].z <= coord[3]+epsilon)
        if (cond1 & cond2 & cond3)
          sameScore += 1
          break
        end
      end
      if sameScore == 0
        push!(nodeCoords, Point3D(coord[1], coord[2], coord[3]))
      end
    end
  end
  nUniqueCoords = length(nodeCoords)
  println("nUniqueCoords: ", nUniqueCoords)

  allIndexes = zeros(Int64, 8, nCells)
  allIndexesVTK = zeros(Int64, 8, nCells)
  jj = 1
  for cell in allCells
    for nNode = 1:8
      kk = 0
      for p in nodeCoords
        if ((cell.nodes[1,nNode] == p.x) & (cell.nodes[2,nNode] == p.y) & (cell.nodes[3,nNode] == p.z))
          allIndexes[nNode, jj] = kk
          break
        end
        kk += 1
      end
    end
    jj += 1
  end

  for i=1:nCells
    for k=1:8
      allIndexesVTK[k,i] = allIndexes[indexTransform[k],i]
    end
  end

  oFile = open("../output/firstTry.vtk", "w")
  write(oFile , "# vtk DataFile Version 3.0\n")
  write(oFile, "some mesh\n")
  write(oFile, "ASCII\n")
  write(oFile, "\n")
  write(oFile, "DATASET UNSTRUCTURED_GRID\n")
  write(oFile, "POINTS " * string(nUniqueCoords) * " float\n")
  nodeCoords_array = zeros(Float64, 3, length(nodeCoords))
  i=1
  for p in nodeCoords
    write(oFile, string(p.x), " ", string(p.y), " ", string(p.z), "\n")
    nodeCoords_array[1, i] = p.x
    nodeCoords_array[2, i] = p.y
    nodeCoords_array[3, i] = p.z
    i+=1
  end
  write(oFile, "\n")

  write(oFile, "CELLS " * string(nCells) * " " * string(nCells*9) * "\n")
  for i=1:size(allIndexes,2)
    write(oFile, "8 ")
    for k = 1:7
      write(oFile, string(allIndexesVTK[k,i]) * " ")
    end
    write(oFile, string(allIndexesVTK[8,i]) * "\n")
  end

  write(oFile, "\n")

  write(oFile, "CELL_TYPES " *string(nCells) * "\n")
  for i=1:nCells
    write(oFile, "11\n")
  end
  close(oFile)
end

function assign_triangles!(oct, allTriangles)
  for tri in allTriangles
      cell = cellContainingPoint(oct, tri.center)
      push!(cell.triangles, tri)
  end
end

function refine_tree(oct)
  for child in oct.children
    if child.isLeaf == 1
      nCellsMax = 20
      for cell in child.cells
        if length(cell.triangles) > nCellsMax
          split_block(child)
        end
      end
    else
      refine_tree(child)
    end
  end
end

function octree_slice_blocks!(oct, coords)
  for child in oct.children
    if child.isLeaf == 1
        for i=1:8
          if -0.001 < child.nodes[1,i] < 0.001
            push!(coords, cell.nodes[2,i])
            push!(coords, cell.nodes[3,i])
          end
        end
    else
      octree_slice!(child, coords)
    end
  end
end

function octree_slice!(oct, coords)
  for child in oct.children
    if child.isLeaf == 1
      for cell in child.cells
        for i=1:8
          if -0.001 < cell.nodes[1,i] < 0.001
            push!(coords, cell.nodes[2,i])
            push!(coords, cell.nodes[3,i])
          end
        end
      end
    else
      octree_slice!(child, coords)
    end
  end
end

function populate_blocks(block::Block)
  for child in block.children
    if child.isLeaf == 1
      insert_cells(child)
    else
      populate_blocks(child)
    end
  end
end

function count_cells(block::Block)
  for child in block.children
    if child.isLeaf == 1
      for cell in child.cells
        if length(cell.triangles) > 0
          println("triangles in cell: ", length(cell.triangles))
        end
      end
    else
      count_cells(child)
    end
  end
end

function all_cells!(block::Block, allCells::Vector{Cell})
  if block.isLeaf == 0
    for child in block.children
      if child.isLeaf == 1
        for cell in child.cells
          push!(allCells, cell)
        end
      else
        all_cells!(child, allCells)
      end
    end
  else
    for cell in block.cells
      push!(allCells, cell)
    end
  end
end

function insert_cells(b::Block)
  lx = b.halfSize[1]*2/b.nx
  ly = b.halfSize[2]*2/b.ny
  lz = b.halfSize[3]*2/b.nz
  for iz = 0:b.nz-1
    for iy = 0:b.ny-1
      for ix = 0:b.nx-1
        cell = Cell(zeros(Float64,3), zeros(Float64,3), zeros(Float64,3,8), 0.0,
                    zeros(Float64,8), Triangle[], false)

        cell.origin[1] = 0.5 * lx + ix * lx + b.origin[1] - b.halfSize[1]
        cell.origin[2] = 0.5 * ly + iy * ly + b.origin[2] - b.halfSize[2]
        cell.origin[3] = 0.5 * lz + iz * lz + b.origin[3] - b.halfSize[3]

        cell.halfSize[1] = lx/2
        cell.halfSize[2] = ly/2
        cell.halfSize[3] = lz/2


        cell.nodes[1:3,1] = [cell.origin[1] - lx/2, cell.origin[2] - ly/2,  cell.origin[3] - lz/2]
        cell.nodes[1:3,2] = [cell.origin[1] - lx/2, cell.origin[2] - ly/2,  cell.origin[3] + lz/2]
        cell.nodes[1:3,3] = [cell.origin[1] - lx/2, cell.origin[2] + ly/2,  cell.origin[3] - lz/2]
        cell.nodes[1:3,4] = [cell.origin[1] - lx/2, cell.origin[2] + ly/2,  cell.origin[3] + lz/2]
        cell.nodes[1:3,5] = [cell.origin[1] + lx/2, cell.origin[2] - ly/2,  cell.origin[3] - lz/2]
        cell.nodes[1:3,6] = [cell.origin[1] + lx/2, cell.origin[2] - ly/2,  cell.origin[3] + lz/2]
        cell.nodes[1:3,7] = [cell.origin[1] + lx/2, cell.origin[2] + ly/2,  cell.origin[3] - lz/2]
        cell.nodes[1:3,8] = [cell.origin[1] + lx/2, cell.origin[2] + ly/2,  cell.origin[3] + lz/2]

        push!(b.cells, cell)
      end
    end
  end

end

function split_block(b::Block)
  if b.isLeaf == 0
    b.cells = Cell[]
  end

  xc1 = [b.origin[1] - b.halfSize[1]/2.0, b.origin[2] - b.halfSize[2]/2.0,  b.origin[3] - b.halfSize[3]/2.0]
  xc2 = [b.origin[1] - b.halfSize[1]/2.0, b.origin[2] - b.halfSize[2]/2.0,  b.origin[3] + b.halfSize[3]/2.0]
  xc3 = [b.origin[1] - b.halfSize[1]/2.0, b.origin[2] + b.halfSize[2]/2.0,  b.origin[3] - b.halfSize[3]/2.0]
  xc4 = [b.origin[1] - b.halfSize[1]/2.0, b.origin[2] + b.halfSize[2]/2.0,  b.origin[3] + b.halfSize[3]/2.0]
  xc5 = [b.origin[1] + b.halfSize[1]/2.0, b.origin[2] - b.halfSize[2]/2.0,  b.origin[3] - b.halfSize[3]/2.0]
  xc6 = [b.origin[1] + b.halfSize[1]/2.0, b.origin[2] - b.halfSize[2]/2.0,  b.origin[3] + b.halfSize[3]/2.0]
  xc7 = [b.origin[1] + b.halfSize[1]/2.0, b.origin[2] + b.halfSize[2]/2.0,  b.origin[3] - b.halfSize[3]/2.0]
  xc8 = [b.origin[1] + b.halfSize[1]/2.0, b.origin[2] + b.halfSize[2]/2.0,  b.origin[3] + b.halfSize[3]/2.0]

  b.children[1] = Block(xc1, b.halfSize/2.0, 1, Array(Block, 8), Cell[], b.nestingLevel+1, b.nx, b.ny, b.nz)
  b.children[2] = Block(xc2, b.halfSize/2.0, 1, Array(Block, 8), Cell[], b.nestingLevel+1, b.nx, b.ny, b.nz)
  b.children[3] = Block(xc3, b.halfSize/2.0, 1, Array(Block, 8), Cell[], b.nestingLevel+1, b.nx, b.ny, b.nz)
  b.children[4] = Block(xc4, b.halfSize/2.0, 1, Array(Block, 8), Cell[], b.nestingLevel+1, b.nx, b.ny, b.nz)
  b.children[5] = Block(xc5, b.halfSize/2.0, 1, Array(Block, 8), Cell[], b.nestingLevel+1, b.nx, b.ny, b.nz)
  b.children[6] = Block(xc6, b.halfSize/2.0, 1, Array(Block, 8), Cell[], b.nestingLevel+1, b.nx, b.ny, b.nz)
  b.children[7] = Block(xc7, b.halfSize/2.0, 1, Array(Block, 8), Cell[], b.nestingLevel+1, b.nx, b.ny, b.nz)
  b.children[8] = Block(xc8, b.halfSize/2.0, 1, Array(Block, 8), Cell[], b.nestingLevel+1, b.nx, b.ny, b.nz)
  b.isLeaf = 0

  for child in b.children
    insert_cells(child)
  end

end

function getOctantContainingPoint(point::Array{Float64,1}, block::Block)
    octant::Int64 = 1
    if (point[1] >= block.origin[1])
      octant += 4
    end
    if (point[2] >= block.origin[2])
      octant += 2
    end
    if (point[3] >= block.origin[3])
      octant += 1
    end
  return octant
end

function blockContainingPoint(block::Block, point::Array{Float64,1})
  if (block.isLeaf == 0)
    oct = getOctantContainingPoint(point, block)
    blockContainingPoint(block.children[oct], point)
  elseif (block.isLeaf == 1)
    return block
  end
end

function cellContainingPoint(oct::Block, point::Array{Float64, 1})
  block = blockContainingPoint(oct, point)
  nx = block.nx
  ny = block.ny
  nz = block.nz
  x = point[1] - block.cells[1].nodes[1,1]
  y = point[2] - block.cells[1].nodes[2,1]
  z = point[3] - block.cells[1].nodes[3,1]

  lx = block.halfSize[1] * 2.0 / nx
  ly = block.halfSize[2] * 2.0 / ny
  lz = block.halfSize[3] * 2.0 / nz

  fx = fld(x, lx)
  fy = fld(y, ly)
  fz = fld(z, lz)

  if fx > (nx-1.0)
      fx = nx - 1.0
  end
  if fy > (ny-1.0)
      fy = ny - 1.0
  end
  if fz > (nz-1.0)
      fz = nz - 1.0
  end

  cellIndex = 1 + fx + fy*nx + fz*nx*ny
  return block.cells[cellIndex]

end

function triLinearInterpolation(cell::Cell, point::Array{Float64,1})

  xd = (point[1] - cell.nodes[1,1]) / (cell.nodes[1,2] - cell.nodes[1,1])
  yd = (point[2] - cell.nodes[2,1]) / (cell.nodes[2,3] - cell.nodes[2,1])
  zd = (point[3] - cell.nodes[3,1]) / (cell.nodes[3,5] - cell.nodes[3,1])

  c00 = cell.data[1] * (1-xd) + cell.data[2] * xd
  c10 = cell.data[5] * (1-xd) + cell.data[6] * xd
  c01 = cell.data[4] * (1-xd) + cell.data[3] * xd
  c11 = cell.data[8] * (1-xd) + cell.data[7] * xd

  c0 = c00*(1-zd) + c10*zd
  c1 = c01*(1-zd) + c11*zd

  c = c0*(1-yd) + c1*yd

  return c
end

end
