@everywhere include("constants.jl")
using DSMC
@everywhere using DSMC

function sendtosimple(iProc, nm, val)
  @spawnat(iProc, (Main, Expr(:(=), $nm, $val)))
end

macro sendto(iProc, nm, val)
  return :(sendtosimple($iProc, $nm, $var))
end

@everywhere S = SphericalBody(128.0, 1e24, 3000.0, 250.0,[500],[18],[1e20])
################################################################################
# load shape model
################################################################################
@everywhere nTriangles, allTriangles, surfaceArea = load_ply_file(mySettings.meshFileName)
@everywhere flux = mySettings.nNewParticlesPerIteration / surfaceArea
@everywhere body = MeshBody(allTriangles, 300.0, 2e24, Float64[flux],
                            Float64[18.0*amu],  Float64[1.e5])
################################################################################
# initialize simulation domain
################################################################################
#domain = initialize_domain(mySettings)

################################################################################
# refine domain and compute volume of cut cells
################################################################################
@everywhere refine_domain(domain, allTriangles, mySettings)
@everywhere allBlocks = Block[]
@everywhere collect_blocks!(domain, allBlocks)
#=
println("distribute it")
for iProc in workers()
  @show(iProc)
  @spawnat(iProc, (Main, Expr(:(=), 'domain', domain)))
  @spawnat(iProc, (Main, Expr(:(=), 'allBlocks', allBlocks)))
end
println("distribute it done")
=#

@everywhere const nParticles = mySettings.nNewParticlesPerIteration
@everywhere const f = nParticles / surfaceArea
#=
for kk in 1:100
  @everywhere x1 = -200.0 + randn() * 200.
  @everywhere x2 = -200.0 + randn() * 200.
  @everywhere x3 = -200.0 + randn() * 200.
  @everywhere coords = [x1,x2,x3]
  @everywhere cell_containing_point(domain, coords)
end
=#

################################################################################
# main loop
################################################################################
#@everywhere const nParticles = mySettings.nNewParticlesPerIteration
#@everywhere const f = nParticles / surfaceArea
@time begin
  for iteration = 0:mySettings.nIterations
    if iteration < 1
      @everywhere insert_new_particles(domain, myPoint)
    end
    if iteration % 5 == 0
      for iProc in workers()
        remotecall_fetch(iProc, save_particles, "../output/particles_iProc_" *string(iProc) * "_" * string(iteration) * ".csv")
      end
    end
    @everywhere compute_macroscopic_params(domain)
    @everywhere time_step(domain, lostParticles, coords)
    @everywhere send_particles_to_cpu(lostParticles)
    @everywhere lostParticles.nParticles = 0
    @show(iteration)
    #readline(STDIN)
  end
end

################################################################################
# save results
################################################################################
save2vtk(domain)
