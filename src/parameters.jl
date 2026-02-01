using Random

export AbstractParameters, Parameters

abstract type AbstractParameters end

mutable struct Parameters <: AbstractParameters
    # System identification
    system_type::Symbol  # :single or :double
    
    # Common parameters for all systems
    KT::Float64
    mass::Float64
    kbond::Float64
    kangle::Float64
    factive::Float64
    dt::Float64
    dt_thermal::Float64
    n_steps::Int64
    thermal_steps::Int64
    logger_steps::Int64
    nthreads::Int64
    γ::Float64
    rcut_nf::Float64
    L::Float64
    
    # Single ring parameters (used when system_type == :single)
    n_monomers::Int64
    n_active::Int64
    
    # Double ring parameters (used when system_type == :double)
    n_monomers_1::Int64
    n_monomers_2::Int64
    n_active_1::Int64
    n_active_2::Int64
end

function Parameters(; system_type::Symbol=:single,
                   KT=1.0, mass=1.0, kbond=30.0, kangle=0.0, factive=1.0,
                   dt=0.01, dt_thermal=0.0, n_steps=100_000, thermal_steps=100_000, logger_steps=500,
                   L=0.0, nthreads=0, γ=2.0, rcut_nf=2.0,
                   # Single ring
                   n_monomers=100, n_active=0,
                   # Double ring
                   n_monomers_1=100, n_monomers_2=100,
                   n_active_1=0, n_active_2=0)
    
    # Auto-calculate box size if not provided
    if L <= 0.0
        if system_type == :single
            L = 3.0 * n_monomers^0.5887
        else  # :double
            total_monomers = n_monomers_1 + n_monomers_2
            L = 3.0 * total_monomers^0.5887
        end
    end

    # Default dt_thermal to dt if not specified
    if dt_thermal <= 0.0
        dt_thermal = dt
    end

    # Auto-detect thread count
    if nthreads <= 0
        nthreads = Threads.nthreads()
    end
    
    # Validate system_type
    if system_type ∉ (:single, :double)
        error("system_type must be :single or :double, got $system_type")
    end
    
    Parameters(
        system_type,
        KT, mass, kbond, kangle, factive, dt, dt_thermal, n_steps, thermal_steps, logger_steps,
        nthreads, γ, rcut_nf, L,
        n_monomers, n_active,
        n_monomers_1, n_monomers_2, n_active_1, n_active_2
    )
end

# Helper functions to get number of particles based on system type
function get_n_particles(params::Parameters)
    if params.system_type == :single
        return params.n_monomers
    else
        return params.n_monomers_1 + params.n_monomers_2
    end
end

function get_n_active(params::Parameters)
    if params.system_type == :single
        return params.n_active
    else
        return params.n_active_1 + params.n_active_2
    end
end

# Get activity vector for all particles
function get_activity_vector(params::Parameters)
    if params.system_type == :single
        n_total = params.n_monomers
        activity = vcat(repeat([true], params.n_active),
                       repeat([false], n_total - params.n_active))
        return shuffle(activity)
    else
        # For double rings, create activity vectors for each ring and combine
        total_1 = params.n_monomers_1
        total_2 = params.n_monomers_2
        
        activity_1 = vcat(repeat([true], params.n_active_1),
                         repeat([false], total_1 - params.n_active_1))
        activity_2 = vcat(repeat([true], params.n_active_2),
                         repeat([false], total_2 - params.n_active_2))
        
        return vcat(shuffle(activity_1), shuffle(activity_2))
    end
end