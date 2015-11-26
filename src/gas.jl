using Distributions
using Types
using Triangles
include("octree.jl")
include("Physical.jl")

function insert_new_particles(oct::Block, body::MeshBody, coords)
 cellID = 0
 nSpecies = length(body.particleMass)
 for tri in body.triangles
   for iSpecies in nSpecies
     N = round(Int, tri.area * body.particleFlux[iSpecies])
     newParticles = Array(Particle, N)
     particleMass = body.particleMass[iSpecies]
     w = body.particleWeight[iSpecies]
     for i=1:N
       pick_point!(tri, coords)
       x = coords[1]
       y = coords[2]
       z = coords[3]
  	    vx, vy, vz = maxwell_boltzmann_flux_v(body.temperature, particleMass)
       vx, vy, vz = rotate_vec_to_pos(vx, vy, vz, x, y, z)
       newParticles[i] = Particle(cellID, x, y, z, vx, vy, vz, particleMass, w)
     end
     assign_particles!(oct, newParticles, coords)
   end
 end
end

function insert_new_particles(oct::Block, body::SphericalBody, coords)
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

function maxwell_boltzmann_flux_v(temperature,mass)
  velmax = 3000.0
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
   theta = 2.0 * pi * rand()
   #polar angle determined from cosine distribution
   phi = asin(sqrt(rand()))
   vx = vel * cos(theta) * sin(phi)
   vy = vel * sin(theta) * sin(phi)
   vz = vel * cos(phi)
#   #Need to rotate vector to particle position
   return vx, vy, vz
end

function rotate_vec_to_pos(vecx,vecy,vecz,posx,posy,posz)
  r = sqrt(posx * posx + posy * posy + posz * posz)
  cosphi = posz / r
  sinphi = sqrt(posx * posx + posy * posy) / r
  costheta = posx / sinphi / r
  sintheta = posy / sinphi / r

  rotated_vectorx = vecx * costheta * cosphi - vecy * sintheta + vecz * costheta * sinphi
  rotated_vectory = vecx * sintheta * cosphi + vecy*costheta + vecz * sintheta * sinphi
  rotated_vectorz = -vecx * sinphi + vecz * cosphi
  return rotated_vectorx, rotated_vectory, rotated_vectorz
end

function compute_macroscopic_params(oct)
  for block in oct.children
    if (block.isLeaf == 1)
      if block.procID == MyID
        compute_params(block)
      end
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

function decompose_particles(particles)
  nParticles = length(particles)
  p_arr = zeros(Float64, 9, nParticles)
  for i=1:nParticles
    @inbounds p_arr[1,i] = particles[i].x
    @inbounds p_arr[2,i] = particles[i].y
    @inbounds p_arr[3,i] = particles[i].z
    @inbounds p_arr[4,i] = particles[i].vx
    @inbounds p_arr[5,i] = particles[i].vy
    @inbounds p_arr[6,i] = particles[i].vz
    @inbounds p_arr[7,i] = particles[i].cellID
    @inbounds p_arr[8,i] = particles[i].mass
    @inbounds p_arr[9,i] = particles[i].weight
  end
  return p_arr
end

function rebuild_particles(p_arr)
  nParticles = size(p_arr, 2)
  particles = Array(Particle, nParticles)
  for i=1:nParticles
    @inbounds particles[i] = Particle(Int(p_arr[7,i]), p_arr[1,i], p_arr[2,i], p_arr[3,i],
                            p_arr[4,i], p_arr[5,i], p_arr[6,i], p_arr[8,i],
                            p_arr[9,i])
  end
  return particles
end

function send_lost_particles(lostParticles, lostIDs)
  @sync begin
    for iProc in unique(lostIDs)
      @async begin
        particles_to_send = Particle[]
        for k = 1:length(lostIDs)
          if lostIDs[k] == iProc
            push!(particles_to_send, lostParticles[k])
          end
        end
        p_arr = decompose_particles(particles_to_send)
        remotecall_fetch(iProc, assign_particles_rem!, p_arr, MyID)
      end
    end
  end
  lostParticles = Particle[]
  lostIDs = Int64[]
end


function time_step(oct::Block, lostParticles, particle_buffer)
  for block in oct.children
    if block.isLeaf == 1
      if block.procID == MyID
        perform_time_step(block, lostParticles, particle_buffer, lostIDs)
      end
    else
      time_step(block, lostParticles, particle_buffer)
    end
  end
  nothing
end

function perform_time_step(b::Block, lostParticles, particle_buffer, lostIDs)
  dt = get_time_step(b)
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
        else
          p.cellID = cell.ID
        end
        wasAssigned, iProc = assign_particle!(p, coords)
        if iProc != MyID
          push!(lostParticles, p)
          push!(lostIDs, iProc)
        end
      end
      splice!(cell.particles, 1:nParticles)
    end
  end
  nothing
end


function assign_particles!(oct, particles, coords)
  for p in particles
    coords[1] = p.x
    coords[2] = p.y
    coords[3] = p.z
    if !is_out_of_bounds(oct, coords)
      foundCell, cell, iProc = cell_containing_point(oct, coords)
      if foundCell && (iProc == MyID)
        p.cellID = cell.ID
        push!(cell.particles, p)
      end
    end
  end
  nothing
end

function assign_particles_rem!(p_arr, senderID)
  coords = zeros(Float64, 3)
  particles = rebuild_particles(p_arr)
  for p in particles
    coords[1] = p.x
    coords[2] = p.y
    coords[3] = p.z
    if !is_out_of_bounds(oct, coords)
      foundCell, cell, iProc = cell_containing_point(oct, coords)
      p.cellID = cell.ID
      push!(cell.particles, p)
    end
  end
  nothing
end

function assign_particle!(p, coords)
    coords[1] = p.x
    coords[2] = p.y
    coords[3] = p.z
    foundCell, cell, iProc = cell_containing_point(oct, coords)
    if foundCell && (iProc == MyID)
      push!(cell.particles, p)
      return true, iProc
    else
      return false, iProc
    end
end

function assign_particle!(block::Block, p, coords)
    coords[1] = p.x
    coords[2] = p.y
    coords[3] = p.z
    foundCell, cell, iProc = cell_containing_point(block, coords)
    if foundCell
      p.cellID = cell.ID
      push!(cell.particles, p)
      return true, iProc
    else
      return false, iProc
    end
end

function constant_weight(dt, S, mass)
  nParticles = 50
  vth = sqrt(8.0 * k_boltz * S.SourceTemperature/pi/mass)
  Flux = pi * S.SourceRadius^2 * S.SourceDensity*vth
  return Flux * dt / nParticles
end

function time_step(temperature,mass)
  ####for test purpose using path lenth of 500 should be based on cell length
  return sqrt(8.0*k_boltz*temperature/pi/mass)/500.0
end

function get_time_step(b::Block)
  return 0.0001
end
