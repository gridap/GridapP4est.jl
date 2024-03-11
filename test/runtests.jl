module GridapP4estTests

using GridapP4est
using MPI
using ArgParse
using Test

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--image-file", "-i"
        help = "Path to the image file that one can use in order to accelerate MPI tests"
        arg_type = String
        default="GridapDistributed.so"
    end
    return parse_args(s)
end

"""
  run_tests(testdir)
"""
function run_tests(testdir)
    repodir = joinpath(@__DIR__,"..")

    parsed_args = parse_commandline()
    image_file_path=parsed_args["image-file"]
    image_file_exists=isfile(image_file_path)

    nprocs_str = get(ENV, "JULIA_GRIDAP_P4EST_TEST_NPROCS","")
    nprocs = nprocs_str == "" ? clamp(Sys.CPU_THREADS, 2, 4) : parse(Int, nprocs_str)
    #mpiexec_args = Base.shell_split("--allow-run-as-root --tag-output") #Base.shell_split(get(ENV, "JULIA_MPIEXEC_TEST_ARGS", ""))
    istest(f) = endswith(f, ".jl") && !(f=="runtests.jl")
    testfiles = sort(filter(istest, readdir(testdir)))
    @time @testset "$f" for f in testfiles
      MPI.mpiexec() do cmd
        if f in ["PoissonUniformlyRefinedOctreeModelsTests.jl"]
          return
          np = [4]
          extra_args = "-s 2 2 -r 2"
        elseif f in ["OctreeDistributedDiscreteModelsTests.jl",
                     "OctreeDistributedDiscreteModelsNoEnvTests.jl",
                     "AdaptivityFlagsMarkingStrategiesTests.jl"]
          return
          np = [4]
          extra_args = ""
        elseif f in ["DarcyNonConformingOctreeModelsTests.jl"]
          return
          np = [1,4]
          extra_args = ""
        elseif f in ["PoissonNonConformingOctreeModelsTests.jl"]
          return
          np = [1,2,4]
          extra_args = ""
        elseif f in ["PoissonAnisotropicOctreeModelsTests.jl"] 
          np = [1,4]
          extra_args = ""
        else
          np = [nprocs]
          extra_args = ""
        end
        for ip in np
          if MPI.MPI_LIBRARY == "OpenMPI" || (isdefined(MPI, :OpenMPI) && MPI.MPI_LIBRARY == MPI.OpenMPI)
            if ! image_file_exists
              cmd2 = `$cmd -n $(ip) --allow-run-as-root --oversubscribe $(Base.julia_cmd()) --project=$(repodir) $(joinpath(testdir, f)) $(split(extra_args))`
            else
              cmd2 = `$cmd -n $(ip) --allow-run-as-root --oversubscribe $(Base.julia_cmd()) --project=$(repodir) -J$(image_file_path) --project=. $(joinpath(testdir, f)) $(split(extra_args))`
            end
          else
            if ! image_file_exists
              cmd2 = `$cmd -n $(ip) $(Base.julia_cmd()) --project=$(repodir) $(joinpath(testdir, f)) $(split(extra_args))`
            else
              cmd2 = `$cmd -n $(ip) $(Base.julia_cmd()) -J$(image_file_path) --project=$(repodir) $(joinpath(testdir, f)) $(split(extra_args))`
            end
          end 
          @show cmd2
          run(cmd2)
        end 
     end 
      @test true
      end
    end  

run_tests(joinpath(@__DIR__,"mpi"))

end # module
