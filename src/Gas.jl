module Gas
using Distributions
using Octree
using Types
using Triangles

include("Physical.jl")

export move!,
       insert_new_particles,
       assign_particles!,
       compute_macroscopic_params,
       time_step,
       time_step2,
       constant_weight

function accelerate!(pos, p::Particle, accl, S)
  r = sqrt(pos[1]^2 +pos[2]^2 + pos[3]^2)
  accl[1] = G * S.SourceMass * p.mass * pos[1] / r^3.0
  accl[2] = G * S.SourceMass * p.mass * pos[2] / r^3.0
  accl[3] = G * S.SourceMass * p.mass * pos[3] / r^3.0
end

function move_RK2!(p::Particle, dt, S)
  pos[1] = p.x
  pos[2] = p.y
  pos[3] = p.z

  acclerate! = (pos, p, a, S)

  rkPos[1] = p.x + p.vx * dt / 2.0
  rkPos[2] = p.y + p.vy * dt / 2.0
  rkPos[3] = p.z + p.vz * dt / 2.0

  rkVel[1] = p.vx + a[1] * dt / 2.0
  rkVel[2] = p.vy + a[2] * dt / 2.0
  rkVel[3] = p.vz + a[3] * dt / 2.0

  acclerate! = (rkPos, p, a, S)

  p.x = p.x + rkVel[1] * dt
  p.y = p.y + rkVel[2] * dt
  p.z = p.z + rkVel[3] * dt

  p.vx = p.vx + a[1] * dt
  p.vy = p.vy + a[2] * dt
  p.vz = p.vz + a[3] * dt
end

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
   end
end

function insert_new_particles_body(oct, allTriangles, f, coords, S)
  particleMass = 18.0 * amu
  w_factor = 1.0
  cellID = 0
  for tri in allTriangles
    N = round(Int, tri.area * f)
    newParticles = Array(Particle, N)
    for i=1:N
      pick_point!(tri, coords)
      x = coords[1]
      y = coords[2]
      z = coords[3]
      #vx = tri.surfaceNormal[1]
      #vy = tri.surfaceNormal[2]
      #vz = tri.surfaceNormal[3]
 	    vx, vy, vz = maxwell_boltzmann_flux_v(S.SourceTemperature, particleMass)
      vx, vy, vz = rotate_vec_to_pos(vx, vy, vz, x, y, z)
      newParticles[i] = Particle(cellID, x, y, z, vx, vy, vz, particleMass, w_factor)
    end
    assign_particles!(oct, newParticles, coords)
  end
end

function insert_new_particles_sphere(oct, nParticles, coords, S, dt)
# #(! function input parameters are changed)
   newParticles = Array(Particle, nParticles)

   particleMass = 28.0*amu
   w_factor = constant_weight(dt,S,particleMass)

  cellID = 0
   for i=1:nParticles
     theta = 2.0*pi*rand()
     phi = acos(2.0*rand()-1.0)
     x = S.SourceRadius*cos(theta)*sin(phi)
 	   y = S.SourceRadius*sin(theta)*sin(phi)
 	   z = S.SourceRadius*cos(phi)
 	   vx, vy, vz = maxwell_boltzmann_flux_v(S.SourceTemperature,particleMass)
     vx, vy, vz = rotate_vec_to_pos(vx, vy, vz, x, y, z)
     newParticles[i] = Particle(cellID, x, y, z, vx, vy, vz, particleMass, w_factor)
  end

  assign_particles!(oct, newParticles, coords)
end

# ############O.J.10-13-15###############################################
# #Maxwwell Boltzmann flux velocity
# ############################
 function maxwell_boltzmann_flux_v(temperature,mass)
  velmax = 5000.0
  beta::Float64 = mass / 2.0 / k_boltz / temperature
  prb::Float64 = 0.0
  r = 1.0

  vel = 0.0
  C = (1.5/beta)^(1.5)*exp(-1.5)
  while r > prb
     vel = rand() * velmax
     a = vel * vel * beta
     prb = vel^3.0 * exp(-a) / C
     r = rand()
   end
   theta = 2.0*pi*rand()
   #polar angle determined from cosine distribution
   phi = asin(sqrt(rand()))
   vx = vel*cos(theta)*sin(phi)
   vy = vel*sin(theta)*sin(phi)
   vz = vel*cos(phi)
#   #Need to rotate vector to particle position
   return vx,vy,vz
 end

function rotate_vec_to_pos(vecx,vecy,vecz,posx,posy,posz)
	r = sqrt(posx*posx+posy*posy+posz*posz)
	cosphi = posz/r
	sinphi = sqrt(posx*posx+posy*posy)/r
	costheta = posx/sinphi/r
	sintheta = posy/sinphi/r

	rotated_vectorx = vecx*costheta*cosphi-vecy*sintheta+vecz*costheta*sinphi
	rotated_vectory = vecx*sintheta*cosphi+vecy*costheta+vecz*sintheta*sinphi
	rotated_vectorz = -vecx*sinphi+vecz*cosphi
    return rotated_vectorx, rotated_vectory, rotated_vectorz
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

function time_step(oct::Block, lostParticles, particle_buffer)
  for block in oct.children
    if block.isLeaf == 1
      perform_time_step(block, lostParticles, particle_buffer)
    else
      time_step(block, lostParticles, particle_buffer)
    end
  end
end

function perform_time_step(b::Block, lostParticles, particle_buffer)
  dt = 0.005
  coords = zeros(Float64, 3)
  pos = zeros(Float64, 3)
  for cell in b.cells
    nParticles = length(cell.particles)
    if nParticles > 0
      if nParticles > length(particle_buffer)
        particle_buffer = Array(Particle, nParticles)
      end
      for i = 1:nParticles
        particle_buffer[i] = cell.particles[i]
      end
      for p in particle_buffer[1:nParticles]
        if p.cellID == cell.ID
          move!(p, dt)
        end
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
      foundCell, cell = cell_containing_point(oct, coords)
      if foundCell
        p.cellID = cell.ID
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
    foundCell, cell = cell_containing_point(oct, coords)
    if foundCell
      p.cellID = cell.ID
      push!(cell.particles, p)
      return true
    else
      return false
    end
end

function constant_weight(dt,S::UniformSource,mass)
  nParticles = 50
  vth = sqrt(8.0*k_boltz*S.SourceTemperature/pi/mass)
  Flux = pi*S.SourceRadius^2*S.SourceDensity*vth
  return Flux*dt/nParticles
end

function time_step(temperature,mass)
  ####for test purpose using path lenth of 500 should be based on cell length
  return sqrt(8.0*k_boltz*temperature/pi/mass)/500.0
end

insert_new_particles = insert_new_particles_body
#insert_new_particles = insert_new_particles_sphere

end
