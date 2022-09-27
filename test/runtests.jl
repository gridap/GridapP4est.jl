module GridapP4estTests

using GridapP4est
using MPI
using ArgParse
using Test

@time @testset "RichardsonSmoothers" begin include("GMG/seq/RichardsonSmoothersTests.jl") end
@time @testset "PatchLinearSolverTests" begin include("GMG/seq/PatchLinearSolverTests.jl") end

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
        if f in ["UniformlyRefinedForestOfOctreesDiscreteModelsTests.jl"]
          np = 4
          extra_args = "-s 2 2 -r 2"
        elseif f in ["MeshHierarchiesTests.jl",
                     "RedistributeToolsTests.jl",
                     "GMGLinearSolversPoissonTests.jl",
                     "OctreeDistributedDiscreteModelsTests.jl"]
          np = 6
          extra_args = ""
        elseif f in ["GMGLinearSolversHDivRTTests.jl",
                     "GMGLinearSolversLaplacianTests.jl",
                     "GMGLinearSolversVectorLaplacianTests.jl"]
          np = 1
          extra_args = ""
        else
          np = nprocs
          extra_args = ""
        end
        if ! image_file_exists
          cmd = `$cmd -n $(np) --allow-run-as-root --oversubscribe $(Base.julia_cmd()) --project=. $(joinpath(testdir, f)) $(split(extra_args))`
        else
          cmd = `$cmd -n $(np) --allow-run-as-root --oversubscribe $(Base.julia_cmd()) -J$(image_file_path) --project=. $(joinpath(testdir, f)) $(split(extra_args))`
        end
        @show cmd
        run(cmd)
        @test true
      end
    end
end

run_tests(@__DIR__)
run_tests(joinpath(@__DIR__, "GMG/mpi"))

end # module
