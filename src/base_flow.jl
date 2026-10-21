# Solve the Falkner-Skan ODE for the base flow via a shooting method 



using DifferentialEquations
using Roots


function falkner_skan!(du,u,β,t)
    du[1] = u[2]
    du[2] = u[3]
    du[3] = -u[1] * u[3] - β * (1 - u[2]^2)
    return nothing
end

function miss(c; β, ηmax)
    prob = ODEProblem(falkner_skan!, [0.0, 0.0, c], (0.0, ηmax), β)
    sol  = solve(prob; reltol=1e-8, abstol=1e-10)
    sol(ηmax; idxs=2) - 1.0
end

function solve_fs(β; ηmax=15.0, guess=1.4)
    root = find_zero(c -> miss(c; β=β, ηmax=ηmax), guess)
    prob  = ODEProblem(falkner_skan!, [0.0, 0.0, root], (0.0, ηmax), β)
    solution   = solve(prob; reltol=1e-10, abstol=1e-10)
    return solution
end

function build_FS_profile(y::AbstractArray; β=0.0, ηmax=20.0)
    sol  = solve_fs(β; ηmax=ηmax, guess=0.47)   # f''(0)=0.4696 for Blasius in this scaling
    Iδ   = ηmax - sol(ηmax; idxs=1)             # δ* in η units (≈1.21678 for Blasius)
    η    = min.(y .* Iδ, ηmax)                  # y is in δ* units
    f    = [sol(ηi; idxs=1) for ηi in η]
    fp   = [sol(ηi; idxs=2) for ηi in η]
    fpp  = [sol(ηi; idxs=3) for ηi in η]
    fppp = -f .* fpp .- β .* (1 .- fp.^2)
    return fp, fppp .* Iδ^2                     # U/U∞ and d²U/dy² in δ* units
end