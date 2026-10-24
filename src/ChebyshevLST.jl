# Entry point for calling the solver from other code, such as the XCALibre.jl cases:
#     Pkg.add(url="https://github.com/James-FK/ChebyshevLST")   # or Pkg.develop(path=...) for a local copy
#     using ChebyshevLST
# The module keeps the solver's helpers (build_FS_profile, solve_fs, …) apart from a case's own.
module ChebyshevLST

using LinearAlgebra
include("base_flow.jl")
include("build_chebyshev_operators.jl")
include("eigen_solve_OSS.jl")
include("chebyshev_interpolate.jl")
include("inlet_eigenfunctions.jl")

export inlet_eigenfunctions

end
