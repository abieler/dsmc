module Octree
using Triangles
using Types
using RayTrace

export split_block,
       insert_cells,
       block_containing_point,
       cell_containing_point,
       populate_blocks,
       out_of_bounds,
       count_cells,
       refine_tree,
       octree_slice!,
       all_cells!,
       is_out_of_bounds,
       cut_cell_volume!

function cut_cell_volume!(oct, pStart, N)
  pRandom = zeros(Float64, 3)
  vRandom = zeros(Float64, 3)
  r = zeros(Float64, 3)
  for block in oct.children
    if block.isLeaf
      for cell in block.cells
        if cell.hasTriangles
          xMin = cell.origin[1] - cell.halfSize[1]
          xMax = cell.origin[1] + cell.halfSize[1]
          yMin = cell.origin[2] - cell.halfSize[2]
          yMax = cell.origin[2] + cell.halfSize[2]
          zMin = cell.origin[3] - cell.halfSize[3]
          zMax = cell.origin[3] + cell.halfSize[3]
          dx = (xMax - xMin) / 1000.0
          dy = (yMax - yMin) / 1000.0
          dz = (zMax - zMin) / 1000.0

          counter = 0
          x = rand(xMin:dx:xMax,N)
          y = rand(yMin:dy:yMax,N)
          z = rand(zMin:dz:zMax,N)

          nInside = 0
          nOutside = 0
          for i=1:N
            pRandom[1] = x[i]
            pRandom[2] = y[i]
            pRandom[3] = z[i]
            counter = nTrianglesIntersects(cell.triangles, pStart, pRandom, vRandom)
            if (counter % 2) == 1
              nInside += 1
            else
              nOutside += 1
            end
          end
          if nInside !== 0
            cell.volume *= (nOutside/(nOutside+nInside))
          end
        end
      end
    else
      cut_cell_volume!(block, pStart, N)
    end
  end
end

function refine(b::Block, nCellsMax)
  if nCellsMax == 0
    split_block(b)
  else
    for cell in b.cells
      if length(cell.triangles) >= nCellsMax
        split_block(b)
      end
    end
  end
end

function refine_tree(oct, nCellsMax=10)
  for block in oct.children
    if block.isLeaf
      refine(block, nCellsMax)
    else
      refine_tree(block, nCellsMax)
    end
  end
end

function count_cells(block::Block)
  for child in block.children
    if child.isLeaf
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
  if !block.isLeaf
    for child in block.children
      if child.isLeaf
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
  volume = lx*ly*lz
  for iz = 0:b.nz-1
    for iy = 0:b.ny-1
      for ix = 0:b.nx-1
        cell = Cell(zeros(Float64,3), zeros(Float64,3), zeros(Float64,3,8), volume,
                    zeros(Float64,8), Triangle[], false, Particle[])

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
  if !b.isLeaf
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

  b.children[1] = Block(xc1, b.halfSize/2.0, true, Array(Block, 8), Cell[], b.nestingLevel+1, b.nx, b.ny, b.nz)
  b.children[2] = Block(xc2, b.halfSize/2.0, true, Array(Block, 8), Cell[], b.nestingLevel+1, b.nx, b.ny, b.nz)
  b.children[3] = Block(xc3, b.halfSize/2.0, true, Array(Block, 8), Cell[], b.nestingLevel+1, b.nx, b.ny, b.nz)
  b.children[4] = Block(xc4, b.halfSize/2.0, true, Array(Block, 8), Cell[], b.nestingLevel+1, b.nx, b.ny, b.nz)
  b.children[5] = Block(xc5, b.halfSize/2.0, true, Array(Block, 8), Cell[], b.nestingLevel+1, b.nx, b.ny, b.nz)
  b.children[6] = Block(xc6, b.halfSize/2.0, true, Array(Block, 8), Cell[], b.nestingLevel+1, b.nx, b.ny, b.nz)
  b.children[7] = Block(xc7, b.halfSize/2.0, true, Array(Block, 8), Cell[], b.nestingLevel+1, b.nx, b.ny, b.nz)
  b.children[8] = Block(xc8, b.halfSize/2.0, true, Array(Block, 8), Cell[], b.nestingLevel+1, b.nx, b.ny, b.nz)
  b.isLeaf = false

  for child in b.children
    insert_cells(child)
  end

end

function octant_containing_point(point::Array{Float64,1}, block::Block)
  if !is_out_of_bounds(block, point)
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
  else
    return -1
  end
end

function block_containing_point(block::Block, point::Array{Float64,1})
  if !block.isLeaf
    oct = octant_containing_point(point, block)
    if oct == -1
      return false, block
    end
    block_containing_point(block.children[oct], point)
  elseif block.isLeaf
    if !is_out_of_bounds(block, point)
      return true, block
    else
      return false, block
    end
  end
end

function cell_containing_point(oct::Block, point::Array{Float64, 1})
  foundBlock, block = block_containing_point(oct, point)
  if foundBlock
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

    cellIndex = round(Int, 1 + fx + fy*nx + fz*nx*ny)
    return true, block.cells[cellIndex]
  else
    return false, block.cells[1]
  end

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

function is_out_of_bounds(oct, r)
  for i=1:3
    if ((r[i] > (oct.origin[i] + oct.halfSize[i])) | (r[i] < (oct.origin[i]-oct.halfSize[i])))
      return true
    end
  end
  return false
end

end
