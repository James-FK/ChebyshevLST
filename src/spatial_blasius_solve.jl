using LinearAlgebra, Plots
include("cheb.jl")
include("base_flow.jl")
N    = 120
Re   = 998     
ω    = 0.1122
yi   = 5.0          
ymax = 60.0

D, ξ = cheb(N)
a = yi*ymax/(ymax - 2yi);  b = 1 + 2a/ymax
y = a*(1 .+ ξ) ./ (b .- ξ)

S  = diagm([0; 1 ./ (1 .- ξ[2:N].^2); 0])
W  = diagm(1 .- ξ.^2);  X = diagm(ξ)
Dk = [I(N+1), D, D^2, D^3, D^4]                                  # Dk[k+1] = D^k
E(k) = ((W*Dk[k+1] - 2k*X*Dk[k] - k*(k-1)*Dk[max(k-1,1)])*S)[2:N,2:N]
E1, E2, E3, E4 = E(1), E(2), E(3), E(4)                          # E4 is old D4

yint = y[2:N];  c = a*(b+1)
M1 = diagm(  c ./ (yint .+ a).^2)
M2 = diagm(-2c ./ (yint .+ a).^3)
M3 = diagm( 6c ./ (yint .+ a).^4)
M4 = diagm(-24c ./ (yint .+ a).^5)
Dy2 = M1^2*E2 + M2*E1
Dy4 = M1^4*E4 + 6M1^2*M2*E3 + (3M2^2 + 4M1*M3)*E2 + M4*E1

# U, Upp = build_FS_profile(.05,1,.1,1.5e-5,yint)
U, Upp = build_FS_profile(yint)

Id = I(N-1);  z = zeros(ComplexF64, N-1, N-1)
C3 = im*Re*diagm(U)
C2 = -(im*ω*Re*Id + 2Dy2)
C1 = im*Re*diagm(Upp) - im*Re*diagm(U)*Dy2
C0 = Dy4 + im*Re*ω*Dy2
A  = [-C3 -C2 -C1 -C0;  Id z z z;  z Id z z;  z z Id z]
λ, V = eigen(A)

idx = findall(α -> 0 < real(α) < 2 && abs(imag(α)) < 0.5 && ω/real(α) < 0.9, λ)
i   = idx[argmin(imag.(λ[idx]))]

α   = λ[i]
println("α = $α")

φ = [0; V[3(N-1)+1:end, i]; 0]
φ = φ / maximum(abs.(φ))
plot(abs.(φ),y, ylims=(0, 10), ylabel="y/δ*", xlabel="|φ|", label="Re=$Re, ω=$ω")



Dy1  = M1*E1
φint = V[3(N-1)+1:end, i]
û    = [0; Dy1*φint; 0]            # u = ∂ψ/∂y
v̂    = -im*α*[0; φint; 0]          # v = -∂ψ/∂x
imx  = argmax(abs.(û))
s    = û[imx]
û ./= s;  v̂ ./= s                  # same complex factor for both

plot(abs.(û),y,ylims = (0,10),ylabel = "y / δ₁  ",xlabel = "|u'|")


# θu = angle.(û);  θv = angle.(v̂)   # phase profiles (plot these vs y)
# plot(θu,y,ylims = (0,5))
