function raytrace(oct, r, r_hat)
  # raytrace through octree from starting point r (x,y,z) unitl line of
  # sight intersects with a surface element, or the ray gets out of the
  # simulation box
  dr = zeros(Float64, 3)
  foundCell, cell = traverse_domain(oct, r, r_hat, dr)
  doesIntersect = intersect(cell.triangles, r_hat, g)
end

function iTriangleIntersect(triangles, pStart, pStop, vStartStop,
                            pI, u, v, w)
  nTriangles = length(triangles)
  iIntersectedTriangle = -1
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

  lIntersectMin = 0.0
  pIntersect = [0.0,0.0,0.0]

  particleStepSize = 0.0
  for i=1:3
    vStartStop[i] = pStop[i] - pStart[i]
    particleStepSize += vStartStop[i] * vStartStop[i]
  end
  particleStepSize = sqrt(particleStepSize)
  r = vStartStop/particleStepSize

  for i=1:nTriangles
    if dot(triangles[i].surfaceNormal, r) > 0.0
      continue
    end
    a = 0.0
    b = 0.0
    @simd for k=1:3
      @inbounds a = a + triangles[i].surfaceNormal[k] * (triangles[i].nodes[k,1] - pStart[k])
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
          @inbounds pI[k] = pStart[k] + rI * r[k]
          @inbounds u[k] = triangles[i].nodes[k,2] - triangles[i].nodes[k,1]
          @inbounds v[k] = triangles[i].nodes[k,3] - triangles[i].nodes[k,1]
          @inbounds w[k] = pI[k] - triangles[i].nodes[k,1]

          @inbounds dot_uv = dot_uv + u[k] * v[k]
          @inbounds dot_uu = dot_uu + u[k] * u[k]
          @inbounds dot_vv = dot_vv + v[k] * v[k]
          @inbounds dot_wu = dot_wu + w[k] * u[k]
          @inbounds dot_wv = dot_wv + w[k] * v[k]
        end
        lIntersect = 0.0
        for kk = 1:3
          lIntersect += (pI[kk] - pStart[kk]) * (pI[kk] - pStart[kk])
        end
        lIntersect = sqrt(lIntersect)
        if lIntersect > particleStepSize
          continue
        end

        divisor = dot_uv * dot_uv - dot_uu * dot_vv
        sI = (dot_uv * dot_wv - dot_vv * dot_wu) / divisor
        tI = (dot_uv * dot_wu - dot_uu * dot_wv) / divisor

        if ((tI >= 0.0) && (sI >= 0.0) && (sI + tI < 1.0))
          if (lIntersectMin == 0.0) | (lIntersect < lIntersectMin)
            lIntersectMin = lIntersect
            iIntersectedTriangle = i
            for kk =1:3
              pIntersect[kk] = pI[kk]
            end
          end
          continue
        end
      end
    end
  end
  # return index of intersected triangle and coordinates of intersection
  # if no intersection, return iIntersectedTriangle == -1
  return iIntersectedTriangle, pIntersect
end

function nTrianglesIntersects(triangles, pStart, pStop, vStartStop,
                              pI, u, v, w)
  # check how many triangles are crossed going from pStart to pStop.
  # vStartStop is the vector pointing from pStart to pStop.
  # it is given to the function as placeholder that is computed inside
  # this function --> this prevents the function to allocate a new Array
  # every time it is called.
  # assuming pStart is in the flow field, pStop is also in the flow field
  # (outside the body) if an even number of triangles is corssed. if an odd
  # number of triangles is crossed, pStop is inside the body
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

  nnorm = 0.0
  for i=1:3
    vStartStop[i] = pStop[i] - pStart[i]
    nnorm += vStartStop[i] * vStartStop[i]
  end
  lRandom = sqrt(nnorm)
  r = vStartStop/lRandom
  counter = 0

  for i=1:nTriangles
    a = 0.0
    b = 0.0
    @simd for k=1:3
      @inbounds a = a + triangles[i].surfaceNormal[k] * (triangles[i].nodes[k,1] - pStart[k])
      @inbounds b = b + triangles[i].surfaceNormal[k] * r[k]
    end
    if ((a != 0.0) && (b != 0.0))
      rI = a / b
      lIntersect = 0.0
      if rI >= 0.0
        dot_uv = 0.0
        dot_uu = 0.0
        dot_vv = 0.0
        dot_wu = 0.0
        dot_wv = 0.0

        for k=1:3
          @inbounds pI[k] = pStart[k] + rI * r[k]
          @inbounds u[k] = triangles[i].nodes[k,2] - triangles[i].nodes[k,1]
          @inbounds v[k] = triangles[i].nodes[k,3] - triangles[i].nodes[k,1]
          @inbounds w[k] = pI[k] - triangles[i].nodes[k,1]

          @inbounds dot_uv = dot_uv + u[k] * v[k]
          @inbounds dot_uu = dot_uu + u[k] * u[k]
          @inbounds dot_vv = dot_vv + v[k] * v[k]
          @inbounds dot_wu = dot_wu + w[k] * u[k]
          @inbounds dot_wv = dot_wv + w[k] * v[k]

          lIntersect += ((pI[k] - pStart[k]) * (pI[k] - pStart[k]))
        end
        lIntersect = sqrt(lIntersect)
        if lIntersect > lRandom
          continue
        end

        divisor = dot_uv*dot_uv - dot_uu * dot_vv
        sI = (dot_uv*dot_wv - dot_vv*dot_wu) / divisor
        tI = (dot_uv*dot_wu - dot_uu*dot_wv) / divisor

        if ((tI >= 0.0) && (sI >= 0.0) && (sI + tI < 1.0))
          counter += 1
          continue
        end
      end
    end
  end
  return counter
end
