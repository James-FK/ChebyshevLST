using LinearAlgebra, Plots
include("base_flow.jl")
include("build_chebyshev_operators.jl")
include("eigen_solve_OSE.jl")

N      = 120
R      = 998/1.7208        # R = √(U∞x/ν): Jordinson's Re_δ* = 998, with δ* = 1.7208ℓ for Blasius
ω      = 0.1122/1.7208     # Jordinson's ω_δ* = 0.1122
y_half = 8.6
y_max  = 100.0


y, Dy1, Dy2, Dy4 = build_chebyshev_operators(N,y_half,y_max)

U, Upp = build_FS_profile(y[2:end-1])

α,û,v̂ =  eigen_solve_OSE(U,Upp,Dy1,Dy2,Dy4,R,ω)


# p_φ = plot(abs.(φ), y, ylims=(0, 17), ylabel="y/ℓ", xlabel="|φ|", label="R=$R, ω=$ω")
# p_u = plot(abs.(û), y, ylims=(0, 17), ylabel="y/ℓ", xlabel="|u'|",label="")
# display(plot(p_φ, p_u, layout=(1, 2)))
