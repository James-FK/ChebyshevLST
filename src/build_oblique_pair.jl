# Physical velocity of the oblique pair Re{q̂₊e^{i(αx+βz−ωt)} + q̂₋e^{i(αx−βz−ωt)}} on the (y, z) grid
# at a given x and t. The −β wave has the same α, û and v̂ but ŵ₋ = −ŵ₊, so the pair is a standing
# wave in z: u and v vary as cos(βz), w as sin(βz). Rows are y, columns are z.
function build_oblique_pair(û,v̂,ŵ,α,β,ω,x,z,t)
    e = exp(im*(α*x - ω*t))
    u = 2real(û*e)    .* cos.(β*z')
    v = 2real(v̂*e)    .* cos.(β*z')
    w = 2real(im*ŵ*e) .* sin.(β*z')
    return u, v, w
end
