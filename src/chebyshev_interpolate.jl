# Barycentric interpolation (Berrut & Trefethen, SIAM Review 2004) of values f at the Chebyshev points
# ξ = cos(πj/N) onto the points ξq. It evaluates the polynomial through the grid values, so it keeps
# the spectral accuracy of the collocation solution, unlike a spline.
function chebyshev_interpolate(f,ξq)
    N = length(f) - 1
    _, ξ = cheb(N)
    w = [0.5; ones(N-1); 0.5] .* (-1).^(0:N)
    return map(ξq) do x
        d = x .- ξ
        j = findfirst(iszero, d)                   # x is a grid point: return its value
        isnothing(j) ? sum(w .* f ./ d)/sum(w ./ d) : f[j]
    end
end
