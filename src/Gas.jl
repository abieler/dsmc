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

function insert_new_particles(oct)
  N = 200
  xInit = rand(oct.origin[1]-oct.halfSize[1]:0.01:oct.origin[1]+oct.halfSize[1], N)
  yInit = rand(oct.origin[2]-oct.halfSize[2]:0.01:oct.origin[2]+oct.halfSize[2], N)
  zInit = rand(oct.origin[3]-oct.halfSize[3]:0.01:oct.origin[3]+oct.halfSize[3], N)

  vxInit = rand(distX, N) * 4
  vyInit = rand(distY, N) .+ 300
  vzInit = rand(distZ, N) * 4
  newParticles = Array(Particle, N)
  for i=1:N
    newParticles[i] = Particle(xInit[i], yInit[i], zInit[i],
                 vxInit[i], vyInit[i], vzInit[i],
                 18.0, 1.0)
  end

  assign_particles!(oct, newParticles)
end

end
