module Gas
using Distributions
using Octree

export Particle,
       move!,
       insert_new_particles


distY = Normal(0,25)
distX = Normal()
distZ = Normal()

type Particle
  x::Float64
  y::Float64
  z::Float64
  vx::Float64
  vy::Float64
  vz::Float64
  mass::Float64
  weight::Float64
end

function move!(p::Particle, dt)
  p.x = p.x + dt * p.vx
  p.y = p.y + dt * p.vy
  p.z = p.z + dt * p.vz
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

end
