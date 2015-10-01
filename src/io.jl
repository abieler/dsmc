
function angle_between(a,b)
  return acos(dot(a, b)/(norm(a) * norm(b)))
end

function cartesian2spherical(v)
  x = v[1]
  y = v[2]
  z = v[3]
  r = sqrt(x^2 + y^2 + z^2)
  theta = asin(z/r)
  phi = atan2(x,y)
  return r, phi, theta
end


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
  triCenters = calculateTriCenters(triangles, nTriangles)
  triAreas = calculateTriAreas(triangles, nTriangles)

  return nTriangles, nodeCoords, triIndices, triangles, n_hat, triCenters, triAreas
end

function rotateCoordinates(nodeCoords, R)
  nNodes = size(nodeCoords)[2]
  nodeCoordsRot = zeros(Float64, 3, nNodes)
  for i=1:nNodes
    vRot = R * nodeCoords[1:3,i]
    for k=1:3
      nodeCoordsRot[k,i] = vRot[k]
    end
  end
  return nodeCoordsRot
end


function rayTrace(triangles::Array{Triangle,1}, r, p)
  nTriangles = length(triangles)
  i = 0
  j = 0
  k = 0

  a = 0.0
  b = 0.0
  rI = 0.0

  dot_uv = 0.0
  dot_uu = 0.0
  dot_vv = 0.0
  dot_wv = 0.0
  dot_wu = 0.0

  divisor = 0.0
  sI = 0.0
  tI = 0.0

  pI = [0.,0.,0.]
  u = [0.,0.,0.]
  v = [0.,0.,0.]
  w = [0.,0.,0.]

  for i=1:nTriangles
    a = 0.0
    b = 0.0
    @simd for k=1:3
      @inbounds a = a + triangles[i].surfaceNormal[k] * (triangles[i].nodes[k,1] - p[k])
      @inbounds b = b + triangles[i].surfaceNormal[k] * r[k]
    end
    if ((a != 0.0) && (b != 0.0))
      rI = a / b
      if rI >= 0.0
        dot_uv = 0.0
        dot_uu = 0.0
        dot_vv = 0.0
        dot_wu = 0.0
        dot_wv = 0.0

        for k=1:3
          @inbounds pI[k] = p[k] + rI * r[k]
          @inbounds u[k] = triangles[i].nodes[k,2] - triangles[i].nodes[k,1]
          @inbounds v[k] = triangles[i].nodes[k,3] - triangles[i].nodes[k,1]
          @inbounds w[k] = pI[k] - triangles[i].nodes[k,1]

          @inbounds dot_uv = dot_uv + u[k]*v[k]
          @inbounds dot_uu = dot_uu + u[k]*u[k]
          @inbounds dot_vv = dot_vv + v[k]*v[k]
          @inbounds dot_wu = dot_wu + w[k]*u[k]
          @inbounds dot_wv = dot_wv + w[k]*v[k]
        end

        divisor = dot_uv*dot_uv - dot_uu * dot_vv
        sI = (dot_uv*dot_wv - dot_vv*dot_wu) / divisor
        tI = (dot_uv*dot_wu - dot_uu*dot_wv) / divisor

        if ((tI >= 0.0) && (sI >= 0.0) && (sI + tI < 1.0))
          if triangles[i].center != p
            return 0
          end
        end
      end

    end
  end
  return 1

end

function buildOctree(triangles, triCenters, triAreas, n_hat, nTriangles)
  allTriangles = Array(Triangle, nTriangles)
  for i=1:nTriangles
    allTriangles[i] = Triangle(i,
                               triCenters[1:3,i],
                               triCenters[1:3,i],
                               triangles[:,:,i],
                               triangles[:,:,i],
                               triAreas[i],
                               n_hat[:,i],
                               n_hat[:,i],
                               0.0,
                               0.0,
                               1,
                               1)
  end

  const halfSize = [3000.0, 3000.0, 3000.0]
  const root = [0.0, 0.0, 0.0]
  const nCellsX = 25
  const nCellsY = 25
  const nCellsZ = 25
  const lxCell = halfSize[1] * 2.0 / nCellsX
  const lyCell = halfSize[2] * 2.0 / nCellsY
  const lzCell = halfSize[3] * 2.0 / nCellsZ

  const nCells = nCellsX * nCellsY * nCellsZ

  cells = Array(Cell, nCells)
  node = zeros(Float64, 3, 8, nCells)
  ii::Int64 = 1
  for iz = 0:nCellsZ-1
    for iy = 0:nCellsY-1
      for ix =0:nCellsX-1
        x = root[1] - halfSize[1] + lxCell/2.0 + ix * lxCell
        y = root[2] - halfSize[2] + lyCell/2.0 + iy * lyCell
        z = root[3] - halfSize[3] + lzCell/2.0 + iz * lzCell
        node[1:3,1,ii] = [x - lxCell/2, y - lyCell/2, z - lzCell/2]
        node[1:3,2,ii] = [x - lxCell/2, y - lyCell/2, z + lzCell/2]
        node[1:3,3,ii] = [x - lxCell/2, y + lyCell/2, z - lzCell/2]
        node[1:3,4,ii] = [x - lxCell/2, y + lyCell/2, z + lzCell/2]
        node[1:3,5,ii] = [x + lxCell/2, y - lyCell/2, z - lzCell/2]
        node[1:3,6,ii] = [x + lxCell/2, y - lyCell/2, z + lzCell/2]
        node[1:3,7,ii] = [x + lxCell/2, y + lyCell/2, z - lzCell/2]
        node[1:3,8,ii] = [x + lxCell/2, y + lyCell/2, z + lzCell/2]

        cells[ii] = Cell(ii, [x,y,z], [lxCell/2, lyCell/2, lzCell/2], node[:,:,ii], Triangle[])
        ii += 1
      end
    end
  end

  octree = Block(root, halfSize, 1, initChildren, cells, nCellsX, nCellsY, nCellsZ)

  for tri in allTriangles
    id = findCellInBlock(octree, tri.center)
    push!(octree.cells[id].triangles, tri)
  end
 return octree
end


function isVis!(oct, r_hat, nTriangles, dr)
  kkk = 0
  for ii=1:3
    dr[ii] = oct.cells[1].halfSize[ii] * r_hat[ii] * 0.8
  end

  newPoint = zeros(Float64, 3)
  newCell = oct.cells[1]
  currentCellID::Int64 = 0
  newCellID::Int64 = 0
  isVisibleBool::Bool = true
  for cell in oct.cells
    currentCellID = cell.id
    for tri in cell.triangles
      if tri.cos_rza > 0.0
        isVisibleBool = rayTrace(cell.triangles, r_hat, tri.center)
        if isVisibleBool == 0
          tri.isVisible = 0
        else
          jjj = 1
          while true
            for ii = 1:3
              @inbounds newPoint[ii] = tri.center[ii] + jjj * dr[ii]
            end
            newCellID = findCellInBlock(oct, newPoint)
            if (newPoint[1] > 3000) | (newPoint[1] < -3000) | (newPoint[2] > 3000) | (newPoint[2] < -3000) | (newPoint[3] > 3000.0) | (newPoint[3] < -3000.0)
              tri.isVisible = 1
              break
            end
            if newCellID != currentCellID
              newCell = oct.cells[newCellID]
              isVisibleBool = rayTrace(newCell.triangles, r_hat, tri.center)
              if isVisibleBool == 0
                tri.isVisible = 0
            	  break
              end
              currentCellID = newCellID
            end
            jjj += 1
          end
        end
      else
        tri.isVisible = 0
      end
        kkk+=1
    end
  end
  return 0
end


function get_spice_stuff(et)
  frame = "67P/C-G_CK"
  rRosetta_km, lt = spice.spkpos("ROSETTA", et, frame, "NONE", "CHURYUMOV-GERASIMENKO")
  rRosetta = Float64[value*1000.0 for value in rRosetta_km]
  rRosetta_hat = rRosetta/norm(rRosetta)

  rSun_km, lt = spice.spkpos("SUN", et, frame, "NONE", "CHURYUMOV-GERASIMENKO")
  rSun = Float64[value for value in rSun_km]
  rSun_hat = rSun/norm(rSun)

  rCG_km, lt = spice.spkpos("CHURYUMOV-GERASIMENKO", et, "ROS_OSIRIS_NAC", "NONE", "ROSETTA")
  rCG = Float64[value*1000.0 for value in rCG_km]
  rCG_hat = rCG/norm(rCG)

  # get rotation matrix
  R_spice = spice.pxform("ROS_OSIRIS_NAC", frame, et)
  R = zeros(Float64, 3, 3)
  i = 1
  for row in R_spice
    for j=1:3
      R[j,i] = row[j]
    end
    i+=1
  end
  return rSun_hat, rCG_hat, rCG, rRosetta_hat, R
end

function rotate_triangles(oct, R, rCG, nodeCoords, triIndices, nTriangles)
  # rotate node coordinates into osiris nac frame
  nodeCoordsRot = rotateCoordinates(nodeCoords, R)
  nodeCoordsRot[1,:] += rCG[1]
  nodeCoordsRot[2,:] += rCG[2]
  nodeCoordsRot[3,:] += rCG[3]
  trianglesRot = build_triangles(nodeCoordsRot, triIndices, nTriangles)
  n_hatRot = calculate_surface_normals(nodeCoordsRot, triIndices, nTriangles)
  triCentersRot = calculateTriCenters(trianglesRot, nTriangles)

  for cell in oct.cells
    for tri in cell.triangles
      tri.centerRot = R * tri.center
      tri.surfaceNormalRot = R * tri.surfaceNormal
      for i=1:3
	      tri.nodesRot[1:3,i] = R * tri.nodes[1:3,i]
      end
      tri.cos_rza = cos(angle_between(n_hatRot[1:3,tri.id], triCentersRot[1:3,tri.id]))
    end
  end
end
