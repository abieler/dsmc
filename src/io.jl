using Triangles

function load_ply_file(fileName::String)
  nNodes::Int64 = 0
  nTriangles::Int64 = 0
  iHeader::Int64 = 0
  i::Int64 = 0
  iFile = open(fileName, "r")
  while !eof(iFile)
    line = readline(iFile)
    if contains(line, "element vertex")
      nNodes = int(split(line, " ")[3])
    elseif contains(line, "element face")
      nTriangles = int(split(line, " ")[3])
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
      triIndices[1,i-iHeader-nNodes] = int(ijk[2])+1
      triIndices[2,i-iHeader-nNodes] = int(ijk[3])+1
      triIndices[3,i-iHeader-nNodes] = int(ijk[4])+1
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
