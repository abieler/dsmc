module Types

export UserSettings,
       Cell,
       Block,
       Point3D,
       Particle,
       Triangle,
       MeshBody,
       SphericalBody,
       _x,
       _y,
       _z

const _x = 1
const _y = 2
const _z = 3


type UserSettings
  nIterations::Int64
  meshFileName::ASCIIString
  domainSizeX::Float64
  domainSizeY::Float64
  domainSizeZ::Float64
  nCellsPerBlockX::Int64
  nCellsPerBlockY::Int64
  nCellsPerBlockZ::Int64
  nMaxRefinementLevel::Int64
  nMaxTrianglesPerCell::Int64
  nNewParticlesPerIteration::Int64
end

# default user settings
UserSettings() = UserSettings(200, "../input/sphere2.ply", 200.0, 200.0, 200.0,
                              5, 5, 5, 3, 1, 5000)

type Triangle
  id::Int64
  center::Array{Float64,1}
  nodes::Array{Float64,2}
  area::Float64
  surfaceNormal::Array{Float64,1}
end

type Particle
  cellID::Int64
  x::Float64
  y::Float64
  z::Float64
  vx::Float64
  vy::Float64
  vz::Float64
  mass::Float64
  weight::Float64
end

type Cell
    ID::Int64
    origin::Vector{Float64}
    halfSize::Vector{Float64}
    nodes::Array{Float64,2}
    volume::Float64
    data::Vector{Float64}
    triangles::Vector{Triangle}
    hasTriangles::Bool
    particles::Vector{Particle}
end

type Block
    origin::Vector{Float64}
    halfSize::Vector{Float64}
    isLeaf::Bool
    children::Vector{Block}
    cells::Vector{Cell}
    nestingLevel::Int64
    nx::Int64
    ny::Int64
    nz::Int64
end

type Point3D
    x::Float64
    y::Float64
    z::Float64
end

type MeshBody
  triangles::Vector{Triangle}
  temperature::Float64
  particleFlux::Vector{Float64}
  particleMass::Vector{Float64}
  particleWeight::Vector{Float64}

end

type SphericalBody
  radius::Float64
  mass::Float64
  density::Float64
  temperature::Float64
  nNewParticles::Vector{Int64}
  particleMass::Vector{Int64}
  particleWeight::Vector{Float64}
end




end
