# PositPR.jl test suite
#
# Run with:
#     julia --project -e 'using Pkg; Pkg.test()'
# or, from the project folder:
#     julia --project test/runtests.jl
#
# The tests compare the library against an independent reference implementation
# (test/reference.jl), EXHAUSTIVELY for 8-bit posits in most cases (all 256
# patterns and all 65 025 operand pairs) and by sampling for 16 and 32 bits.
#
# Each arithmetic operation is checked at increasing levels of strictness:
#   1. it does not throw;
#   2. when the exact result is representable, it is returned exactly;
#   3. correct rounding (to nearest, ties to even), as required by the Posit
#      standard;
#   4. error of at most 1 ULP ("faithful" rounding).
# This makes it clear at which level a failure happens.
#
# Functions with `while` loops (eps, nextPosit, exp, log, sin, ...) run in a
# separate process with a time limit (test/watchdog.jl): an infinite loop becomes
# a test failure instead of hanging the whole suite.

using Test
using Random

include("reference.jl")

@testset "PositPR.jl" begin
    include("test_api.jl")
    include("test_encoding_decoding.jl")
    include("test_conversions.jl")
    include("test_comparisons.jl")
    include("test_operations.jl")
    include("test_functions.jl")
    include("test_random.jl")
    include("test_quire.jl")
end

stop_watchdog()
