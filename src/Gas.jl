
module Gas
using Distributions
using Octree
using Types
using Physical

include("constants.jl")

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

function insert_new_particles(oct, nParticles, coords)

  particleMass = 18.0 * amu
  w_factor = 1.0
  xMin = oct.origin[1] + oct.halfSize[1] * 0.95
  xMax = oct.origin[1] + oct.halfSize[1] * 0.99
  yMin = oct.origin[2] - oct.halfSize[2] * 0.1
  yMax = oct.origin[2] + oct.halfSize[2] * 0.1
  zMin = oct.origin[3] - oct.halfSize[3] * 0.1
  zMax = oct.origin[3] + oct.halfSize[3] * 0.1

  dx = (xMax - xMin) / 1000.0
  dy = (yMax - yMin) / 1000.0
  dz = (zMax - zMin) / 1000.0


  xInit = rand(xMin:dx:xMax, N)
  yInit = rand(yMin:dy:yMax, N)
  zInit = rand(zMin:dz:zMax, N)

  vxInit = -ones(Float64, N)
  vyInit = zeros(Float64, N)
  vzInit = zeros(Float64, N)

  newParticles = Array(Particle, N)
  for i=1:N
    newParticles[i] = Particle(xInit[i], yInit[i], zInit[i],
                 vxInit[i], vyInit[i], vzInit[i],
                 particleMass, w_factor)
  end
  assign_particles!(oct, newParticles, coords)
end


##########O.J 10-13-15#################################################
#insert_new_particles_sphere
#insert particles uniformly about sphere
#possible update insert particles from each spherical surface
#####################
function insert_new_particles_sphere(oct, N, coords)
# need radius of body it will be used in multiple places
# (! function input parameters are changed)
# coords is a 3-Vector passed to assign_particles unmodified
  body_radius = 3.0
  mass_N2 = 28.0 * amu
  w_factor = 1.0
  newParticles = Array(Particle, N)
  for i=1:N
    theta = 2.0 * pi * rand()
    phi = acos(2.0 * rand() - 1.0)

    xInit = body_radius * cos(theta) * sin(phi)
	  yInit = body_radius * sin(theta) * sin(phi)
	  zInit = body_radius * cos(phi)

	  vxInit, vyInit, vzInit = maxwell_boltzmann_flux_v(theta, phi)

    newParticles[i] = Particle(xInit, yInit, zInit, vxInit, vyInit, vzInit,
                               mass_N2, w_factor)
  end
  assign_particles!(oct, newParticles, coords)
end

############O.J.10-13-15###############################################
#Maxwwell Boltzmann flux velocity
############################
function maxwell_boltmann_flux_v(thetaPos, phiPos)
  velmax = 3000.0e3
  temperature = source_temperature()
  beta = mass_N2 / 2.0 / k_boltz / temperature
  while r > prb
    vel = rand() * velmax
    a = vel * vel * beta
    prb = vel^3.0 * exp(-a) / ((1.5 / beta)^(1.5) * exp(-1.5))
  end
  theta = 2.0 * pi * rand()
  #polar angle determined from cosine distribution
  phi = asin(sqrt(rand()))
  vx = vel * cos(theta) * sin(phi)
  vy = vel * sin(theta) * sin(phi)
  vz = vel * cos(phi)
  #Need to rotate vector to particle position
  return vx, vy, vz
end

############O.J.10-13-15###############################################
#Source Temperature
#Change as need for distribution about source based on particle coordinate on surface node face
############################
function source_temperature()
  return 150.0
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
  coords = zeros(Float64, 3)
  for cell in b.cells
    nParticles = length(cell.particles)
    if nParticles > 0
      for p in copy(cell.particles)
        move!(p, 0.01)
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
