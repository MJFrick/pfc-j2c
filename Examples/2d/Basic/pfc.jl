mutable struct ReArrays
    n
    NL
end

mutable struct CmxArrays
    nk
    NL_k
end

mutable struct halfReArrays
    k2
    C2
    P
    Q
end

mutable struct dims
    Nx
    Ny
end

mutable struct IntParams
    step
end

mutable struct doubleParams
    dt
    dx
    dy
    Bx
    t
    v
    M
    tm
end

mutable struct transforms
    forward
    back
end

mutable struct state
    ReArrays
    CmxArrays
    halfReArrays
    IntParams
    doubleParams
    transforms
end

function Calc_halfReals(s::state)
    Lx = s.dims.Nx * s.doubleParams.dx
    Ly = s.dims.Ny * s.doubleParams.dy
    for i in 0:((s.dims.Nx >> 1))
        kx2 = i < ((s.dims.Ny>>1) + 1) ? (2π * i / Lx)^2 : (2π * (i-s.dims.Nx) / Lx)^2
        for j in 0:(s.dims.Ny-1)
            ky2 = j < ((s.dims.Ny>>1) + 1) ? (2π * j / Ly)^2 : (2π * (j-s.dims.Ny) / Ly)^2
            s.halfReArrays.k2[i+1,j+1] = (kx2 + ky2)
        end
    end

    s.halfReArrays.C2 .= s.doubleParams.Bx.*(2.0.*s.halfReArrays.k2 .- s.halfReArrays.k2.*s.halfReArrays.k2)

    for i in 0:((s.dims.Nx >> 1))
        for j in 0:(s.dims.Ny-1)
            k2 = s.halfReArrays.k2[i+1,j+1]
            M = s.doubleParams.M
            dt = s.doubleParams.dt
            Lambda = 1.0 - s.halfReArrays.C2[i+1,j+1]
            s.halfReArrays.P[i+1,j+1] = -M * dt * k2 * Lambda / (1.0 + M * dt* k2 * Lambda)
            s.halfReArrays.Q[i+1,j+1] = -M * dt * k2 / (1.0 + M * dt * k2 * Lambda)
        end
    end

    return nothing
end #function

function Calc_NL(s::state)
    norm = 1.0 / (s.dims.Nx * s.dims.Ny)
    s.ReArrays.NL .= -0.5 .* s.doubleParams.t .* s.ReArrays.n .* s.ReArrays.n .+ third .* s.doubleParams.v .* s.ReArrays.n .* s.ReArrays.n .* s.ReArrays.n
    s.CmxArrays.NL_k = s.transforms.forward * s.ReArrays.NL
    s.CmxArrays.NL_k .*= norm
    return nothing
end #function

function step(s::state)
    norm = 1.0 / (s.dims.Nx * s.dims.Ny)
    s.CmxArrays.nk = s.transforms.forward * s.ReArrays.n
    s.CmxArrays.nk .*= norm
    Calc_NL(s)
    s.CmxArrays.nk .+= s.halfReArrays.P .* s.CmxArrays.nk .+ s.halfReArrays.Q .* s.CmxArrays.NL_k
    s.ReArrays.n = s.transforms.back * s.CmxArrays.nk

    s.IntParams.step += 1
    s.doubleParams.tm += s.doubleParams.dt
    return nothing
end #function

function seed(s::state, q::Float64, A::Float64)
    for i in 0:((s.dims.Nx - 1))
        for j in 0:((s.dims.Ny - 1))
            x = i * s.doubleParams.dx
            y = j * s.doubleParams.dy
            s.ReArrays.n[i+1,j+1] += exp(-((i - (s.dims.Nx >>1))^2+(j - (s.dims.Ny >> 1))^2) / 300) * 2.0 * A * (cos(q * x) + 2.0 * cos(0.5*(q * x)) * cos(0.5 * sqrt(3) * q * y))
        end
    end

    return nothing
end #function