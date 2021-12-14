using Test
using Sunny
using Random

Random.seed!(1111)

# Idea taken from StaticArrays.jl
enabled_tests = lowercase.(ARGS)
function addtests(fname)
    key = lowercase(splitext(fname)[1])
    if isempty(enabled_tests) || key in enabled_tests
        include(fname)
    end
end

addtests("test_lattice.jl")
addtests("test_interactions.jl")
addtests("test_units.jl")
addtests("test_ewald.jl")
addtests("test_symmetry.jl")
addtests("test_metropolis.jl")
addtests("test_fourier.jl")
addtests("test_dynamics.jl")
