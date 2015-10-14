
module Gas
using Distributions
using Octree
using Types

export move!,
       insert_new_particles,
       assign_particles!,
       compute_macroscopic_params,
       time_step


distY = Normal(0,25)
distX = Normal()
distZ = Normal()


function move!(p::Particle, dt)
  p.x = p.x + dt * p.vx
  p.y = p.y + dt * p.vy
  p.z = p.z + dt * p.vz
end

function next_pos!(p::Particle, dt, pos)
  pos[1] = p.x + dt * p.vx
  pos[2] = p.y + dt * p.vy
  pos[3] = p.z + dt * p.vz
end


function gas_surface_collisions!(block)
    for child in block.children
      if child.isLef

    counter = nTrianglesIntersects(cell.triangles, r, pStart, pRandom, vRandom)

end

function insert_new_particles(oct, coords)
  N = 100
  xMin = 19.5
  xMax = 19.9
  yMin = -5
  yMax = 5
  zMin = -5
  zMax = 5

  xInit = rand(xMin:0.01:xMax, N)
  yInit = rand(yMin:0.01:yMax, N)
  zInit = rand(zMin:0.01:zMax, N)

  vxInit = -ones(Float64, N)
  vyInit = zeros(Float64, N)
  vzInit = zeros(Float64, N)
  newParticles = Array(Particle, N)
  for i=1:N
    newParticles[i] = Particle(xInit[i], yInit[i], zInit[i],
                 vxInit[i], vyInit[i], vzInit[i],
                 18.0, 1.0)
  end

  assign_particles!(oct, newParticles, coords)
end

function compute_macroscopic_params(oct)
  for block in oct.children
    if block.isLeaf == 1
      compute_params(block)
    else
      compute_macroscopic_params(block)
    end
  end
end

function compute_params(block)
  for cell in block.cells
    density = length(cell.particles)/cell.volume
    if isnan(density)
      density = 0.0
    end
    cell.data[1] = density
  end
end

function time_step(oct, lostParticles)
  for block in oct.children
    if block.isLeaf == 1
      perform_time_step(block, lostParticles)
    else
      time_step(block, lostParticles)
    end
  end
end

function perform_time_step(b::Block, lostParticles)
  dt = 0.01
  coords = zeros(Float64, 3)
  pos = zeros(Float64, 3)
  for cell in b.cells
    nParticles = length(cell.particles)
    if nParticles > 0
      for p in copy(cell.particles)
        if cell.hasTriangles
          next_pos!(p, dt, pos)
        end
        move!(p, dt)
        wasAssigned = assign_particle!(b, p, coords)
        if !wasAssigned
          push!(lostParticles, p)
        end
      end
      splice!(cell.particles, 1:nParticles)
    end
  end
end

function assign_particles!(oct, particles, coords)
  for p in particles
    coords[1] = p.x
    coords[2] = p.y
    coords[3] = p.z
    if !is_out_of_bounds(oct, coords)
      foundCell, cell = cellContainingPoint(oct, coords)
      if foundCell
        push!(cell.particles, p)
      end
    end

  end
  return 0
end

function assign_particle!(oct, p, coords)
    coords[1] = p.x
    coords[2] = p.y
    coords[3] = p.z
    foundCell, cell = cellContainingPoint(oct, coords)
    if foundCell
      push!(cell.particles, p)
      return true
    else
      return false
    end
end

end
