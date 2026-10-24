# Spatial Orr–Sommerfeld–Squire solve for an oblique wave q̂(y)exp(i(αx + βz − ωt)), with ω and β real
# and α the complex eigenvalue. Squire is forced by v̂ but does not feed back into Orr–Sommerfeld (the
# coupled system is block lower-triangular), so α and v̂ come from Orr–Sommerfeld alone and η̂ is the
# Squire response to that v̂.
function eigen_solve_OSS(U,Up,Upp,Dy1,Dy2,Dy4,Dy2_dirichlet,R,ω,β)
    n  = size(Dy1,1)
    Id = I(n)
    Z  = zeros(ComplexF64, n, n)

    # Orr–Sommerfeld for v̂: the 2D coefficients with D² → D² − β²
    Dβ2 = Dy2 - β^2*Id
    Dβ4 = Dy4 - 2β^2*Dy2 + β^4*Id              # (D² − β²)², built from Dy4 so it keeps the clamped BCs
    C3  = im*R*Diagonal(U)
    C2  = -(im*ω*R*Id + 2Dβ2)
    C1  = im*R*Diagonal(Upp) - im*R*Diagonal(U)*Dβ2
    C0  = Dβ4 + im*R*ω*Dβ2
    L_companion = [-C3 -C2 -C1 -C0;  Id Z Z Z;  Z Id Z Z;  Z Z Id Z]
    αs, V = eigen(L_companion)

    idx  = findall(α -> 0 < real(α) < 2 && abs(imag(α)) < 0.5 && ω/real(α) < 0.9, αs)
    k_ts = idx[argmin(imag.(αs[idx]))]

    α    = αs[k_ts]
    println("α = $α")

    # Squire for η̂, forced by v̂: (α² + αS1 + S0)η̂ = −Kv̂
    v̂int = V[3n+1:4n, k_ts]
    S1   = im*R*Diagonal(U)
    S0   = -Dy2_dirichlet + (β^2 - im*ω*R)*Id
    K    = im*β*R*Diagonal(Up)
    η̂int = -(α^2*Id + α*S1 + S0) \ (K*v̂int)

    peak(v) = v[argmax(abs.(v))]       # complex value where |v| is largest

    Dv̂   = Dy1*v̂int
    k²   = α^2 + β^2
    û    = [0; im*(α*Dv̂ - β*η̂int)/k²; 0]    # from continuity and η = ∂u/∂z − ∂w/∂x
    v̂    = [0; v̂int; 0]
    ŵ    = [0; im*(β*Dv̂ + α*η̂int)/k²; 0]
    s    = peak(û)                     # same complex factor for all three
    û  ./= s
    v̂  ./= s
    ŵ  ./= s
    return α, û, v̂, ŵ
end
