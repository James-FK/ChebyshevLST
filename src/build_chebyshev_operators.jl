# compute chebyshev differentiation matrix in Julia, this is a port of the MATLAB code 
# from the excellent textbook by Trefethen: Spectral Methods in Matlab. https://doi.org/10.1137/1.9780898719598

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



function build_chebyshev_operators(N,y_half,y_max)

    D, ξ  = cheb(N)
    inner = 2:N                        # interior points (walls excluded)

    a = y_half*y_max/(y_max - 2y_half)
    b = 1 + 2a/y_max
    c = a*(b+1)
    y = a*(1 .+ ξ) ./ (b .- ξ)

    #This section is the Trefethen clamped boundary condition trick introduced in Chapter 14 
    S = Diagonal([0; 1 ./ (1 .- ξ[inner].^2); 0])
    W = Diagonal(1 .- ξ.^2)
    X = Diagonal(ξ)
    D2, D3, D4 = D^2, D^3, D^4
    Dξ1 = ((W*D  - 2X)*S)[inner, inner]
    Dξ2 = ((W*D2 - 4X*D  - 2I)*S)[inner, inner]
    Dξ3 = ((W*D3 - 6X*D2 - 6D)*S)[inner, inner]
    Dξ4 = ((W*D4 - 8X*D3 - 12D2)*S)[inner, inner]

    yint   = y[inner]
    ξ_y    = Diagonal(  c ./ (yint .+ a).^2)
    ξ_yy   = Diagonal(-2c ./ (yint .+ a).^3)
    ξ_yyy  = Diagonal( 6c ./ (yint .+ a).^4)
    ξ_yyyy = Diagonal(-24c ./ (yint .+ a).^5)
    Dy1 = ξ_y*Dξ1
    Dy2 = ξ_y^2*Dξ2 + ξ_yy*Dξ1
    Dy4 = ξ_y^4*Dξ4 + 6ξ_y^2*ξ_yy*Dξ3 + (3ξ_yy^2 + 4ξ_y*ξ_yyy)*Dξ2 + ξ_yyyy*Dξ1

    # Dirichlet-only second derivative (f = 0 at both ends, f' free) for second-order equations
    # such as Squire, using the standard boundary row/column removal from Trefethen Chapter 7
    Dy2_dirichlet = ξ_y^2*D2[inner, inner] + ξ_yy*D[inner, inner]


    return y, Dy1, Dy2, Dy4, Dy2_dirichlet
end



# inverse of the map y = a(1 + ξ)/(b − ξ) in build_chebyshev_operators, keep the two in step
function y_to_ξ(y,y_half,y_max)
    a = y_half*y_max/(y_max - 2y_half)
    b = 1 + 2a/y_max
    return (b .* y .- a) ./ (y .+ a)
end