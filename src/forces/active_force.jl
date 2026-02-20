export ActiveTangentForce

mutable struct ActiveTangentForce
    system_type::Symbol
    n_total::Int64
    n_monomers_1::Int64
    n_monomers_2::Int64
    k::Float64
    f_active::Float64
    forces::Vector{SVector{3,Float64}}
    tang_vecs::Vector{SVector{3,Float64}}
end

function ActiveTangentForce(; system_type::Symbol, k::Float64, f_active::Float64,
                            n_monomers_1::Int64, n_monomers_2::Int64=0)

    if system_type == :single
        n_total = n_monomers_1
    else  # :double
        n_total = n_monomers_1 + n_monomers_2
    end

    forces = [SVector{3,Float64}(0.0, 0.0, 0.0) for _ in 1:n_total]
    tang_vecs = [SVector{3,Float64}(0.0, 0.0, 0.0) for _ in 1:n_total]

    ActiveTangentForce(system_type, n_total, n_monomers_1, n_monomers_2,
                      k, f_active, forces, tang_vecs)
end

function force_active_tangent(coords_i::SVector{3,Float64},
                             coords_j::SVector{3,Float64},
                             coords_k::SVector{3,Float64},
                             k::Float64, boundary)
    ba = vector(coords_j, coords_i, boundary)
    bc = vector(coords_j, coords_k, boundary)
    ba_norm = norm(ba)
    bc_norm = norm(bc)
    
    # set tang vec for b
    tang_vecj = bc / bc_norm - ba / ba_norm
    tang_vecj = normalize(tang_vecj)
    
    cross_ba_bc = ba × bc
    if iszero(cross_ba_bc)
        zf = zero(k .* ba)
        return (tang_vecj, zf, zf, zf)
    end
    
    pa = ba × cross_ba_bc
    pa = pa / norm(pa)
    pc = -bc × cross_ba_bc
    pc = pc / norm(pc)
    
    cosθ = (ba ⋅ bc) / (ba_norm * bc_norm)
    cosθ = clamp(cosθ, -1.0, 1.0)  # Prevent domain errors
    θ = acos(cosθ)
    
    angle_term = k * sin(θ)
    fa = (angle_term / ba_norm) * pa
    fc = (angle_term / bc_norm) * pc
    fb = -fa - fc
    
    return (tang_vecj, fa, fb, fc)
end

function AtomsCalculators.forces!(fs, sys, inter::ActiveTangentForce;
                                   neighbors=nothing, step_n=0, n_threads=Threads.nthreads(),
                                   buffers=nothing, needs_vir=false)
    fill!(inter.forces, SVector{3,Float64}(0.0, 0.0, 0.0))

    if inter.system_type == :single
        compute_forces_single_ring!(inter, sys)
    else
        compute_forces_double_ring!(inter, sys)
    end

    # Apply active forces to active particles
    for id in 1:inter.n_total
        if sys.atoms[id].active == true
            inter.forces[id] += inter.tang_vecs[id] .* inter.f_active
        end
    end

    # Add to force accumulator
    for i in 1:length(fs)
        fs[i] += inter.forces[i]
    end

    return nothing
end

function compute_forces_single_ring!(inter::ActiveTangentForce, sys)
    n = inter.n_monomers_1

    # Ring-periodic triplets: every monomer is a center
    js = collect(1:n)
    is = [mod(j - 2, n) + 1 for j in 1:n]   # previous: [n, 1, 2, ..., n-1]
    ks = [mod(j, n) + 1 for j in 1:n]        # next:     [2, 3, 4, ..., 1]

    tv_forces = force_active_tangent.(sys.coords[is], sys.coords[js], sys.coords[ks],
                                     (inter.k,), (sys.boundary,))

    for idx in 1:n
        inter.tang_vecs[js[idx]] = tv_forces[idx][1]
        inter.forces[is[idx]] += tv_forces[idx][2]
        inter.forces[js[idx]] += tv_forces[idx][3]
        inter.forces[ks[idx]] += tv_forces[idx][4]
    end
end

function compute_forces_double_ring!(inter::ActiveTangentForce, sys)
    n1 = inter.n_monomers_1
    n2 = inter.n_monomers_2

    # First ring: monomers 1..n1, offset=0
    compute_ring_forces!(inter, sys, 0, n1)

    # Second ring: monomers (n1+1)..(n1+n2), offset=n1
    compute_ring_forces!(inter, sys, n1, n2)
end

function compute_ring_forces!(inter::ActiveTangentForce, sys, offset::Int, n::Int)
    if n < 3
        return
    end

    # Ring-periodic triplets: every monomer is a center
    js = [offset + j for j in 1:n]
    is = [offset + mod(j - 2, n) + 1 for j in 1:n]   # previous
    ks = [offset + mod(j, n) + 1 for j in 1:n]        # next

    tv_forces = force_active_tangent.(sys.coords[is], sys.coords[js], sys.coords[ks],
                                     (inter.k,), (sys.boundary,))

    for idx in 1:n
        inter.tang_vecs[js[idx]] = tv_forces[idx][1]
        inter.forces[is[idx]] += tv_forces[idx][2]
        inter.forces[js[idx]] += tv_forces[idx][3]
        inter.forces[ks[idx]] += tv_forces[idx][4]
    end
end