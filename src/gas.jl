
#=
function insert_new_particles(domain::Block, coords)
 N = 500
 R = 10.0
 procID = 0
 w = 1.0
 particleMass = 18
 v = zeros(Float64,3)
 pos = zeros(Float64, 3)
 for i=1:N
   theta = 2.0 * pi * rand()
   phi = acos(2.0 * rand() - 1.0)
   x = 800.0 + R * cos(theta) * sin(phi)
   y = R * sin(theta) * sin(phi)
   z = 20.0 + R * cos(phi)
   vx = -5000.0
   vy = 0.0
   vz = 0.0
   assign_particle!(domain, procID, x, y, z, vx, vy, vz,
                    particleMass, w, coords)
 end
end
=#


function insert_new_particles(domain::Block, coords)
  k = 0
  for j in 1:5
    for i in 1:100
        procID = 0
        w = 1.0
        particleMass = 18
        x = 400.0 + 20 * randn() 
        y = 0.0 +  3 * randn()
        z = 0.0
        vx = -5000.0
        vy = 0.0
        vz = 0.0
        assign_particle!(domain, procID, x, y, z, vx, vy, vz,
                         particleMass, w, coords)
        k += 1
    end
  end
end

function insert_new_particles(domain::Block, body::MeshBody, coords)
 procID = 0
 nSpecies = length(body.particleMass)
 for tri in body.triangles
   for iSpecies in nSpecies
     N = round(Int, tri.area * body.particleFlux[iSpecies])
     particleMass = body.particleMass[iSpecies]
     w = body.particleWeight[iSpecies]
     for i=1:N
        pick_point!(tri, coords)
        x = coords[1]
        y = coords[2]
        z = coords[3]
  	    vx, vy, vz = maxwell_boltzmann_flux_v(body.temperature,
                                              particleMass)
        vx, vy, vz = rotate_vec_to_pos(vx, vy, vz, x, y, z)
        assign_particle!(domain, procID, x, y, z, vx, vy, vz,
                         particleMass, w, coords)
     end
   end
 end
end

function insert_new_particles(domain::Block, body::SphericalBody, coords)
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
    #newParticles[i] = Particle(cellID, x, y, z, vx, vy, vz, particleMass, w_factor)
    newParticles[i] = Particle(cellID, x, y, z, -vx, -vy, -vz, particleMass, w_factor)
  end

  assign_particles!(domain, newParticles, coords)
end
function accelerate!(pos, p::Particle, accl, body)
  r = sqrt(pos[1]^2 +pos[2]^2 + pos[3]^2)
  accl[1] = G * body.mass * p.mass * pos[1] / r^3.0
  accl[2] = G * body.nass * p.mass * pos[2] / r^3.0
  accl[3] = G * body.mass * p.mass * pos[3] / r^3.0
end

function move_RK2!(p::Particle, dt, body)
  pos = zeros(Float64, 3)
  accl = zeros(Float64, 3)
  rkPos = zeros(Float64, 3)
  rkVel = zeros(Float64, 3)

  pos[1] = p.x
  pos[2] = p.y
  pos[3] = p.z

  acclerate! = (pos, p, accl, body)

  rkPos[1] = p.x + p.vx * dt / 2.0
  rkPos[2] = p.y + p.vy * dt / 2.0
  rkPos[3] = p.z + p.vz * dt / 2.0

  rkVel[1] = p.vx + accl[1] * dt / 2.0
  rkVel[2] = p.vy + accl[2] * dt / 2.0
  rkVel[3] = p.vz + accl[3] * dt / 2.0

  acclerate! = (rkPos, p, accl, body)

  p.x = p.x + rkVel[1] * dt
  p.y = p.y + rkVel[2] * dt
  p.z = p.z + rkVel[3] * dt

  p.vx = p.vx + a[1] * dt
  p.vy = p.vy + a[2] * dt
  p.vz = p.vz + a[3] * dt
end

function move!(p::Particles, dt)
  for i=1:p.nParticles
    @inbounds p.x[i] = p.x[i] + dt * p.vx[i]
    @inbounds p.y[i] = p.y[i] + dt * p.vy[i]
    @inbounds p.z[i] = p.z[i] + dt * p.vz[i]
  end
  #fid = open("../output/traj.csv", "a")
  #write(fid, string(p.x[1]),",", string(p.y[1]), ",", string(p.z[1]), "\n")
  #close(fid)
end

function stop!(p::Particles)
  for i=1:p.nParticles
    @inbounds p.vx[i] = 0.0
    @inbounds p.vy[i] = 0.0
    @inbounds p.vz[i] = 0.0
  end
end

function move!(cell::Cell, dt)
  if length(cell.triangles) > 0
    #gas_surface_collisions!(cell, dt)
    #println("checked for coll")
    move!(cell.particles, dt)
  else
    move!(cell.particles, dt)
  end
end

function gas_surface_collisions!(cell::Cell, dt)
  for i in cell.particles.nParticles
    for k =1:3
      nowpos[k] = 0.0
      nextpos[k] = 0.0
      vStartStop[k] = 0.0
      pI[k] = 0.0
      u[k] = 0.0
      v[k] = 0.0
      w[k] = 0.0
    end
    nextpos[1] = cell.particles.x[i] + dt * cell.particles.vx[i]
    nextpos[2] = cell.particles.y[i] + dt * cell.particles.vy[i]
    nextpos[3] = cell.particles.z[i] + dt * cell.particles.vz[i]
    nowpos[1] = cell.particles.x[i]
    nowpos[2] = cell.particles.y[i]
    nowpos[3] = cell.particles.z[i]
    iTri, pIntersect = iTriangleIntersect(cell.triangles, nowpos, nextpos,
                                          vStartStop, pI, u, v, w)
    if iTri != -1
      println("r_old       : ", sqrt(nowpos[1]^2 + nowpos[2]^2 + nowpos[3]^2))
      println("collision! r: ", sqrt(pIntersect[1]^2 + pIntersect[2]^2 + pIntersect[3]^2))
      println("  ")

      cell.particles.x[i] = pIntersect[1]
      cell.particles.y[i] = pIntersect[2]
      cell.particles.z[i] = pIntersect[3]

      cell.particles.vx[i] = 0.0
      cell.particles.vy[i] = 0.0
      cell.particles.vz[i] = 0.0
      data2CSV(cell.particles, i)
    else
      cell.particles.x[i] = cell.particles.x[i] + dt * cell.particles.vx[i]
      cell.particles.y[i] = cell.particles.y[i] + dt * cell.particles.vy[i]
      cell.particles.z[i] = cell.particles.z[i] + dt * cell.particles.vz[i]
    end
    #stop!(cell.particles)
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

function compute_macroscopic_params(domain)
  for block in domain.children
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
    density = length(cell.particles.nParticles)/cell.volume
    if isnan(density)
      density = 0.0
    end
    cell.data[1] = density
  end
end

function decompose_particles(p::Particles)
  nParticles = p.nParticles
  p_arr = zeros(Float64, 9, nParticles)
  for i=1:nParticles
    @inbounds p_arr[1,i] = p.x[i]
    @inbounds p_arr[2,i] = p.y[i]
    @inbounds p_arr[3,i] = p.z[i]
    @inbounds p_arr[4,i] = p.vx[i]
    @inbounds p_arr[5,i] = p.vy[i]
    @inbounds p_arr[6,i] = p.vz[i]
    @inbounds p_arr[7,i] = p.procID[i]
    @inbounds p_arr[8,i] = p.mass[i]
    @inbounds p_arr[9,i] = p.weight[i]
  end
  return p_arr
end

function rebuild_particles(p_arr)
  nParticles = size(p_arr, 2)
  p = Particles(nParticles)
  for i=1:nParticles
    @inbounds p.x[i] = p_arr[1,i]
    @inbounds p.y[i] = p_arr[2,i]
    @inbounds p.z[i] = p_arr[3,i]
    @inbounds p.vx[i] = p_arr[4,i]
    @inbounds p.vy[i] = p_arr[5,i]
    @inbounds p.vz[i] = p_arr[6,i]
    @inbounds p.procID[i] = p_arr[7,i]
    @inbounds p.mass[i] = p_arr[8,i]
    @inbounds p.weight[i] = p_arr[9,i]
  end
  p.nParticles = nParticles
  return p
end

function send_particles_to_cpu(lostParticles)
  if MyID > 1
    doSend = (MyID in workers()) & (lostParticles.nParticles > 0)
    if doSend == true
      lostIDs = lostParticles.procID
      uniqueIDs = unique(lostIDs)
      splice!(uniqueIDs, findfirst(uniqueIDs, 0))
      #println("sending ", lostParticles.nParticles, " to iprocs: ", uniqueIDs)
      @sync begin
        for iProc in uniqueIDs
          @async begin
            particles_to_send = Particles(count(i->i==iProc, lostIDs))
            for k = 1:length(lostIDs)
              if lostIDs[k] == iProc
                add!(lostParticles, particles_to_send, k)
              end
            end
            p_arr = decompose_particles(particles_to_send)
            remotecall_fetch(iProc, assign_particles_rem!, p_arr, MyID)
          end
        end
      end
    end
  end
  nothing
end

function time_step(domain::Block, lostParticles, coords)
  for block in domain.children
    if block.isLeaf == 1
      if block.procID == MyID
        perform_time_step(block, lostParticles, coords)
      end
    else
      time_step(block, lostParticles, coords)
    end
  end
  nothing
end

function perform_time_step(b::Block, lostParticles, coords)
  dt = get_time_step(b)
  for cell in b.cells
    if cell.particles.nParticles > 0
      move!(cell, dt)
      assign_particles!(domain, cell.particles, coords, lostParticles)
    end
  end
  nothing
end

function assign_particles_rem!(p_arr, senderID)
  coords = zeros(Float64, 3)
  #println("got ", size(p_arr, 2), " particles from ", senderID)
  p = rebuild_particles(p_arr)
  for i in 1:p.nParticles
    coords[1] = p.x[i]
    coords[2] = p.y[i]
    coords[3] = p.z[i]
    foundCell, cell, iProc = cell_containing_point(domain, coords)
    add!(p, cell.particles, i)
  end
  nothing
end

function remove!(p::Particles, i::Int64)
  splice!(p.x, i)
  splice!(p.y, i)
  splice!(p.z, i)
  splice!(p.vx, i)
  splice!(p.vy, i)
  splice!(p.vz, i)
  splice!(p.procID, i)
  splice!(p.mass, i)
  splice!(p.weight, i)
  p.nParticles -= 1
end

function add!(from::Particles, to::Particles, i::Int64)
  n = to.nParticles + 1
  nMax = length(to.x)
  if n > nMax
    allocate_particles!(to)
  end
  to.nParticles += 1
  to.x[n] = from.x[i]
  to.y[n] = from.y[i]
  to.z[n] = from.z[i]
  to.vx[n] = from.vx[i]
  to.vy[n] = from.vy[i]
  to.vz[n] = from.vz[i]
  to.procID[n] = from.procID[i]
  to.mass[n] = from.mass[i]
  to.weight[n] = from.weight[i]
end

function assign_particles!(domain::Block, p::Particles, coords, lostParticles)
  rmIndexes = Int64[]
  for i=1:p.nParticles
    coords[1] = p.x[i]
    coords[2] = p.y[i]
    coords[3] = p.z[i]
    if !is_out_of_bounds(domain, coords)
      foundCell, cell, iProc = cell_containing_point(domain, coords)
      p.procID[i] = iProc
      if iProc != MyID
        add!(p, lostParticles, i)
        push!(rmIndexes, i)
      else
        add!(p, cell.particles, i)
        push!(rmIndexes, i)
      end
    else
      push!(rmIndexes, i)
    end
  end

  for i in sort!(rmIndexes, rev=true)
    remove!(p, i)
  end
end

function allocate_particles!(p::Particles; N=25)
  append!(p.procID, zeros(Int64, N))
  append!(p.x, zeros(Float64, N))
  append!(p.y, zeros(Float64, N))
  append!(p.z, zeros(Float64, N))
  append!(p.vx, zeros(Float64, N))
  append!(p.vy, zeros(Float64, N))
  append!(p.vz, zeros(Float64, N))
  append!(p.mass, zeros(Float64, N))
  append!(p.weight, zeros(Float64, N))
  nothing
end

function assign_particle!(domain, procID, x, y, z, vx, vy, vz,
                          particleMass, w, coords)
  coords[1] = x
  coords[2] = y
  coords[3] = z
  foundCell, cell, iProc = cell_containing_point(domain, coords)

  if (foundCell) & (iProc == MyID)
    maxParticles = length(cell.particles.x)
    iParticle = cell.particles.nParticles + 1
    cell.particles.nParticles += 1
    if iParticle > maxParticles
      allocate_particles!(cell.particles)
    end
    cell.particles.procID[iParticle] = MyID
    cell.particles.x[iParticle] = x
    cell.particles.y[iParticle] = y
    cell.particles.z[iParticle] = z
    cell.particles.vx[iParticle] = vx
    cell.particles.vy[iParticle] = vy
    cell.particles.vz[iParticle] = vz
    cell.particles.mass[iParticle] = particleMass
    cell.particles.weight[iParticle] = w
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
