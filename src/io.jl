using Triangles
using Octree

function load_ply_file(fileName::String)
  nNodes::Int64 = 0
  nTriangles::Int64 = 0
  iHeader::Int64 = 0
  i::Int64 = 0
  iFile = open(fileName, "r")
  while !eof(iFile)
    line = readline(iFile)
    if contains(line, "element vertex")
      nNodes = parse(Int, split(line, " ")[3])
    elseif contains(line, "element face")
      nTriangles = parse(Int, split(line, " ")[3])
    elseif contains(line, "end_header")
      iHeader = i
      break
    end
    i += 1
  end
  close(iFile)

  println("nNodes: ", nNodes)
  println("nTriangles : ", nTriangles)
  println("iHeader: ", iHeader)

  nodeCoords = zeros(Float64, 3, nNodes)
  triIndices = zeros(Int64, 3, nTriangles)

  i = 0
  iFile = open(fileName, "r")
  while !eof(iFile)
    line = readline(iFile)
    if iHeader < i <= iHeader+nNodes
      xyz = matchall(r"-?\d+(\.\d+)?", line)
      nodeCoords[1,i-iHeader] = float(xyz[1])
      nodeCoords[2,i-iHeader] = float(xyz[2])
      nodeCoords[3,i-iHeader] = float(xyz[3])
    elseif i > iHeader+nNodes
      ijk = matchall(r"(\d+)", line)
      triIndices[1,i-iHeader-nNodes] = parse(Int, ijk[2])+1
      triIndices[2,i-iHeader-nNodes] = parse(Int, ijk[3])+1
      triIndices[3,i-iHeader-nNodes] = parse(Int, ijk[4])+1
    end
    i += 1
  end
  close(iFile)

  triangles = build_triangles(nodeCoords, triIndices, nTriangles)
  n_hat = calculate_surface_normals(nodeCoords, triIndices, nTriangles)
  triCenters = calculate_tri_centers(triangles, nTriangles)
  triAreas = calculate_tri_areas(triangles, nTriangles)
  n_hat = calculate_surface_normals(nodeCoords, triIndices, nTriangles)

  allTriangles = Array(Triangle, nTriangles)
  for i=1:nTriangles
    tri = Triangle(i, triCenters[1:3, i], triangles[1:3,1:3,i],
                   triAreas[i], n_hat[1:3,i])
    allTriangles[i] = tri
  end
  return nTriangles, allTriangles
end

function save2vtk(oct)
  println("saving simulation domain to disk")
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

  oFile = open("../output/domain.vtk", "w")
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


  write(oFile, "\n")
  write(oFile, "CELL_DATA " * string(nCells) * "\n")
  write(oFile, "SCALARS density float\n")
  write(oFile, "LOOKUP_TABLE default\n")

  for i=1:nCells
    write(oFile, string(allCells[i].data[1]) * "\n")
  end
  close(oFile)
end
