module RayTrace

using Types

export raytrace,
       nTrianglesIntersects

function raytrace(oct, r, r_hat)
  # raytrace through octree from starting point r (x,y,z) unitl line of
  # sight intersects with a surface element, or the ray gets out of the
  # simulation box
  dr = zeros(Float64, 3)
  foundCell, cell = traverse_domain(oct, r, r_hat, dr)
  doesIntersect = intersect(cell.triangles, r_hat, g)
end

function nTrianglesIntersects(triangles, p, pRandom, vRandom)
  nTriangles = length(triangles)
  #println("start intersect for ", nTriangles, " triangles")
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

  for i=1:3
    vRandom[i] = pRandom[i] - p[i]
  end
  lRandom = norm(vRandom)
  r = vRandom/lRandom
  counter = 0

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
        lIntersect = norm(pI .- p)
        #println("tri.center: ", triangles[i].center)
        #println("pRandom: ", p)
        #println("pI: ", pI)
        #println("lRandom: ", lRandom)
        #println("lIntersect: ", lIntersect)

        if lIntersect > lRandom
          #println("not intersected")
          continue
        end

        divisor = dot_uv*dot_uv - dot_uu * dot_vv
        sI = (dot_uv*dot_wv - dot_vv*dot_wu) / divisor
        tI = (dot_uv*dot_wu - dot_uu*dot_wv) / divisor

        if ((tI >= 0.0) && (sI >= 0.0) && (sI + tI < 1.0))
          counter += 1
          #println("counter += 1")
          continue
        end
        #println()
      end
    end
  end
  return counter
end

function traverse_domain(oct::Block, r, r_hat, dr)
  cell = cellContainingPoint(oct, r)
  dummyCell = deepcopy(cell)

  while true
    if out_of_bounds(oct, r)
      return false, dummyCell
    end
    cell = cellContainingPoint(oct, r)
    if cell.hasTriangles
      return true, cell
    end

    for k=1:3
      dr[k] = cell.halfSize[k] * r_hat[k]
      r[k] = r[k] + dr[k]
    end
  end
end


end
