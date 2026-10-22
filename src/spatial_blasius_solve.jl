using LinearAlgebra, Plots
include("base_flow.jl")
include("build_chebyshev_operators.jl")
include("eigen_solve_OSE.jl")

N      = 120
Re_δ   = 998
ω      = 0.1122
y_half = 5.0
y_max  = 60.0


y, Dy1, Dy2, Dy4 = build_chebyshev_operators(N,y_half,y_max)

U, Upp = build_FS_profile(y[2:end-1])

α,û,v̂ =  eigen_solve_OSE(U,Upp,Dy1,Dy2,Dy4,Re_δ,ω)


# p_φ = plot(abs.(φ), y, ylims=(0, 10), ylabel="y/δ*", xlabel="|φ|", label="Re=$Re_δ, ω=$ω")
# p_u = plot(abs.(û), y, ylims=(0, 10), ylabel="y/δ*", xlabel="|u'|",label="")
# display(plot(p_φ, p_u, layout=(1, 2)))
