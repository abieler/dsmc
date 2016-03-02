using DSMC
@everywhere using DSMC
@everywhere include("constants.jl")


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

@everywhere sayHi(domain)
readline(STDIN)
################################################################################
# main loop
################################################################################
@everywhere const nParticles = mySettings.nNewParticlesPerIteration
@everywhere const f = nParticles / surfaceArea
@time begin
  for iteration = 0:mySettings.nIterations
    #@everywhere insert_new_particles(domain, body, myPoint)
    if iteration < 1
      #@everywhere insert_new_particles(domain, myPoint)
      @everywhere insert_new_particles(domain, myPoint)
    end
    @everywhere compute_macroscopic_params(domain)
    @everywhere time_step(domain, lostParticles, coords)
    @everywhere send_particles_to_cpu(lostParticles)
    @everywhere lostParticles.nParticles = 0
    if iteration % 10 == 0
      for iProc in workers()
        remotecall_fetch(iProc, save_particles, "../output/particles_iProc_" *string(iProc) * "_" * string(iteration) * ".csv")
      end
    end
    @show(iteration)
  end
end

################################################################################
# save results
################################################################################
@everywhere save2vtk(domain)
