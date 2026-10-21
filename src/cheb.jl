# compute chebyshev differentiation matrix in Julia, this is a port of the MATLAB code 
# from the excellent textbook by Trefethen: Spectral Methods in Matlab. https://doi.org/10.1137/1.9780898719598
using LinearAlgebra

function cheb(N)
    if N == 0
        return zeros(1, 1), [1.0]
    end
    x = cos.(π * (0:N) / N)
    c = [2; ones(N-1); 2] .* (-1).^(0:N)
    dX = x .- x'
    D = (c ./ c') ./ (dX + I )
    D = D - Diagonal(vec(sum(D, dims=2)))
    return D, x
end