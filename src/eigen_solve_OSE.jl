function eigen_solve_OSE(U,Upp,Dy1,Dy2,Dy4,Re_δ,ω)
    n = size(Dy1,1)
    Id = I(n)
    Z  = zeros(ComplexF64, n, n)
    C3 = im*Re_δ*Diagonal(U)
    C2 = -(im*ω*Re_δ*Id + 2Dy2)
    C1 = im*Re_δ*Diagonal(Upp) - im*Re_δ*Diagonal(U)*Dy2
    C0 = Dy4 + im*Re_δ*ω*Dy2
    L_companion = [-C3 -C2 -C1 -C0;  Id Z Z Z;  Z Id Z Z;  Z Z Id Z]
    αs, V = eigen(L_companion)

    idx  = findall(α -> 0 < real(α) < 2 && abs(imag(α)) < 0.5 && ω/real(α) < 0.9, αs)
    k_ts = idx[argmin(imag.(αs[idx]))]

    α    = αs[k_ts]
    println("α = $α")

    peak(v) = v[argmax(abs.(v))]       # complex value where |v| is largest

    φint = V[3n+1:4n, k_ts]
    φ    = [0; φint; 0]
    φ    = φ / peak(φ)
    û    = [0; Dy1*φint; 0]            # u = ∂ψ/∂y
    v̂    = -im*α*[0; φint; 0]          # v = -∂ψ/∂x
    s    = peak(û)                     # same complex factor for both
    û  ./= s
    v̂  ./= s
    return α, û, v̂
end