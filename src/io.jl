type Triangle
  id::Int64
  center::Array{Float64,1}
  centerRot::Array{Float64,1}
  nodes::Array{Float64,2}
  nodesRot::Array{Float64,2}
  area::Float64
  surfaceNormal::Array{Float64,1}
  surfaceNormalRot::Array{Float64,1}
  cos_sza::Float64
  cos_rza::Float64
  isVisible::Int64
  isSunlit::Int64
end

type Cell
  id::Int64
  origin::Array{Float64,1}
  halfSize::Array{Float64,1}
  nodes::Array{Float64,2}
  triangles::Array{Triangle, 1}
end

type Block
  origin::Array{Float64, 1}
  halfSize::Array{Float64, 1}
  isLeaf::Int64
  children::Array{Block, 1}
  cells::Array{Cell,1}
  nCellsX::Int64
  nCellsY::Int64
  nCellsZ::Int64
end

initChildren = Array(Block, 8)
initCells = Array(Cell, 1)

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

function calculate_surface_normals(nodeCoords, triIndices, nTriangles)

  n_hat = zeros(Float64, 3, nTriangles)
  vi = zeros(Float64, 3)
  vj = zeros(Float64, 3)
  vk = zeros(Float64, 3)

  for ii=1:nTriangles
    i = triIndices[1, ii]
    j = triIndices[2, ii]
    k = triIndices[3, ii]

    vi = vec(nodeCoords[1:3, i])
    vj = vec(nodeCoords[1:3, j])
    vk = vec(nodeCoords[1:3, k])
    r = cross(vj-vi, vk-vi)
    r = r/norm(r)
    for kk = 1:3
      n_hat[kk, ii] = r[kk]
    end
  end

  return n_hat

end

function build_triangles(nodeCoords, triIndices, nTriangles)
  triangles = zeros(Float64, 3, 3, nTriangles)
  for i=1:nTriangles
    for j=1:3
      for k=1:3
        triangles[k,j,i] = nodeCoords[k,triIndices[j,i]]
      end
    end
  end
  return triangles
end

function calculateTriCenters(triangles, nTriangles)
  triCenters = zeros(Float64, 3, nTriangles)
  for i=1:nTriangles
    for j=1:3
      triCenters[j,i] = sum(triangles[j,1:3,i])/3.0
    end
  end
  return triCenters
end


function calculateTriAreas(triangles, nTriangles)
  triAreas = zeros(Float64, nTriangles)
  for i=1:nTriangles
    P = vec(triangles[1:3,2,i] - triangles[1:3,1,i])
    Q = vec(triangles[1:3,3,i] - triangles[1:3,1,i])
    S = sqrt(sum(cross(P,Q).^2))
    triAreas[i] = 0.5 * S
  end
  return triAreas
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


function findCellInBlock(block::Block, point::Array{Float64, 1})
  const nx = block.nCellsX
  const ny = block.nCellsY
  const nz = block.nCellsZ

  x = point[1] - block.cells[1].nodes[1,1]
  y = point[2] - block.cells[1].nodes[1,2]
  z = point[3] - block.cells[1].nodes[1,3]
  lx = block.halfSize[1] * 2 / nx
  ly = block.halfSize[2] * 2 / ny
  lz = block.halfSize[3] * 2 / nz
  fx = fld(x, lx)
  fy = fld(y, ly)
  fz = fld(z, lz)

  if fx > (nx-1)
    fx = nx-1
  end
  if fy > (ny-1)
    fy = ny-1
  end
  if fz > (nz-1)
    fz = nz-1
  end

  cellIndex = 1 + fx + fy*nx + fz*nx * ny
  return int(cellIndex)
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

function isSun!(oct, r_hat, nTriangles, dr)
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
      if tri.cos_sza > 0.0
        # check if triangle is visible from the Sun (in shadow or not)
        isVisibleBool = rayTrace(cell.triangles, r_hat, tri.center)
        if isVisibleBool == 0
          tri.isSunlit = 0
        # if not interseption from triangles within own cell, check triangles
        # of the next cell
        else
          jjj = 1
          while true
            for ii = 1:3
              # propagate ray by dr, find new cell and check against those
              # triangles
              @inbounds newPoint[ii] = tri.center[ii] + jjj * dr[ii]
            end
            newCellID = findCellInBlock(oct, newPoint)
            # if the ray leaves the 3 km bounding box no more triangles can
            # intercept the ray --> it is set to visible
            if (newPoint[1] > 3000) | (newPoint[1] < -3000) | (newPoint[2] > 3000) | (newPoint[2] < -3000) | (newPoint[3] > 3000.0) | (newPoint[3] < -3000.0)
              tri.isSunlit = 1
              break
            end
            if newCellID != currentCellID
              newCell = oct.cells[newCellID]
              isVisibleBool = rayTrace(newCell.triangles, r_hat, tri.center)
              if isVisibleBool == 0
                tri.isSunlit = 0
      	        break
              end
              currentCellID = newCellID
            end
            jjj += 1
          end
        end
      else
	      tri.isSunlit = 0
      end
      kkk+=1
    end
  end
  return 0
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

function compute_sza(oct, r_hat)
  for cell in oct.cells
    for tri in cell.triangles
      tri.cos_sza = cos(angle_between(r_hat, tri.surfaceNormal))
    end
  end
end

function compute_rza(oct, r_hat)
  for cell in oct.cells
    for tri in cell.triangles
      tri.cos_rza = cos(angle_between(r_hat, tri.surfaceNormal))
    end
  end
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

function my_pmap(f, lst, oct, nodeCoords, triIndices, nTriangles)
    np = nprocs()  # determine the number of processes available
    n = length(lst)
    results = cell(n)
    i = 1
    # function to produce the next work item from the queue.
    # in this case it's just an index
    nextidx() = (idx=i; i+=1; idx)
    @sync begin
        for p=1:np
            if p != myid() || np == 1
                println(myid())
                @async begin
                    while true
                        idx = nextidx()
                        #println(myid(), "   ")
                        #println("idx     : ", idx)
                        #println("lst[idx]: ", lst[idx])
                        #rSun_hat, rCG_hat, rCG, R = get_spice_stuff(lst[idx])
                        #compute_sza(oct, rSun_hat)
                        #rotate_triangles(oct, R, rCG, nodeCoords, triIndices, nTriangles)
                        if idx > n
                            break
                        end
                        #results[idx] = sum(remotecall_fetch(p, doItAll, oct))
                        results[idx] = bii()
                        println(myid(), " xxx ", results[idx])
                        readline(STDIN)
                    end
                end
            end
        end
    end
    results
end

function bii()
  sleep(1)
  rand(1:100)
end

function parallel_spice(et)
  rSun_hat, rCG_hat, rCG, R = get_spice_stuff(et)
  compute_sza(oct, rSun_hat)
  rotate_triangles(oct, R, rCG, nodeCoords, triIndices, nTriangles)

  return sqrt(rCG[1]^2 + rCG[2]^2 + rCG[3]^3)

end

function save_result(oct)
  fH2O = zeros(Float64, nTriangles)
  for i=1:nTriangles
    r, lon, lat = cartesian2spherical(triCenters[1:3,i])
    fH2O[i] = 2.0 * maximum([0.3, sin(lat)])
  end

  sza = zeros(Float64, nTriangles)
  rza = zeros(Float64, nTriangles)
  isSunlit = zeros(Int64, nTriangles)
  isVisible = zeros(Int64, nTriangles)
  for cell in oct.cells
    for tri in cell.triangles
      sza[tri.id] = tri.cos_sza
      rza[tri.id] = tri.cos_rza
      isSunlit[tri.id] = tri.isSunlit
      isVisible[tri.id] = tri.isVisible
    end
  end
  np.save(string(t) * "_sza.npy", sza)
  np.save(string(t) * "_rza.npy", rza)
  np.save(string(t) * "_isSunlit.npy", isSunlit)
  np.save(string(t) * "_isVisible.npy", isVisible)
  np.save(string(t) * "_fH2O.npy", fH2O)
end
