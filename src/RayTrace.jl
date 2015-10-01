module RayTrace

export raytrace

function raytrace(oct, r, r_hat)
  dr = zeros(Float64, 3)
  foundCell, cell = traverse_domain(oct, r, r_hat, dr)
  doesIntersect = intersect(cell.triangles, r_hat, r)

end

function intersect(triangles, r, p)

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
          return true
        end
      end

    end
  end
  return false

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

function out_of_bounds(oct, r)
  for i=1:3
    if (r[i] > oct.halfSize[i]) | (r[i] < -oct.halfSize[i])
      return true
    end
  end
  return false
end

end
