# Interface to the CFD: solve for the spatial TS wave at frequency ω and spanwise wavenumber β (β = 0 gives
# the 2D wave) and return α with û, v̂, ŵ interpolated onto the inlet cell-centre heights y_inlet [m] at
# the inlet station x_inlet [m]. R = √(U∞x/ν), and α, ω, β are in units of ℓ = √(νx/U∞) = x/R.
# β_H is the Hartree parameter of the Falkner–Skan base flow.
function inlet_eigenfunctions(N,R,ω,β,x_inlet,y_inlet; β_H=0.0, y_half=8.6, y_max=100.0)
    y, Dy1, Dy2, Dy4, Dy2_dirichlet = build_chebyshev_operators(N,y_half,y_max)
    U, Upp, Up = build_FS_profile(y[2:end-1]; β=β_H)
    α,û,v̂,ŵ = eigen_solve_OSS(U,Up,Upp,Dy1,Dy2,Dy4,Dy2_dirichlet,R,ω,β)

    ℓ       = x_inlet/R
    ξ_inlet = y_to_ξ(min.(y_inlet ./ ℓ, y_max),y_half,y_max)   # above y_max the modes have decayed to 0
    on_inlet(f) = chebyshev_interpolate(f, ξ_inlet)
    return α, on_inlet(û), on_inlet(v̂), on_inlet(ŵ)
end
