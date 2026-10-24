using LinearAlgebra, Plots
include("base_flow.jl")
include("build_chebyshev_operators.jl")
include("eigen_solve_OSS.jl")
include("build_oblique_pair.jl")

N      = 120
R      = 998/1.7208        # R = √(U∞x/ν): Jordinson's Re_δ* = 998, with δ* = 1.7208ℓ for Blasius
ω_TS   = 0.1122/1.7208     # frequency of the fundamental 2D TS wave (Jordinson's ω_δ* = 0.1122)
ω      = ω_TS/2            # subharmonic: the oblique pair runs at half the TS frequency
β      = 0.12              # spanwise wavenumber of the pair ±β (not the Hartree β in base_flow.jl)
y_half = 8.6
y_max  = 100.0


y, Dy1, Dy2, Dy4, Dy2_dirichlet = build_chebyshev_operators(N,y_half,y_max)

U, Upp, Up = build_FS_profile(y[2:end-1])

α,û,v̂,ŵ = eigen_solve_OSS(U,Up,Upp,Dy1,Dy2,Dy4,Dy2_dirichlet,R,ω,β)

z = range(0, 2π/β, length=65)                               # one spanwise wavelength
u′,v′,w′ = build_oblique_pair(û,v̂,ŵ,α,β,ω,0.0,z,0.0)        # velocity of the pair at x = 0, t = 0


p_mode = plot([abs.(û) abs.(v̂) abs.(ŵ)], y, ylims=(0, 14), ylabel="y/ℓ", label=["|û|" "|v̂|" "|ŵ|"])
p_pair = contourf(z, reverse(y), reverse(u′, dims=1), ylims=(0, 17), xlabel="z/ℓ", title="u′ at x = 0, t = 0")
display(plot(p_mode, p_pair, layout=(1, 2)))

savefig(plot(p_mode, p_pair, layout=(1, 2)),"3d_spanwise_contour.svg")