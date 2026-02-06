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
    traj_interval::Int64
    metric_mode::Symbol        # :fixed or :logspaced
    metric_interval::Int64     # fixed interval for metrics (0 = same as traj_interval)
    metric_npoints::Int64      # number of sample points for logspaced
    msd_com::Bool              # compute center-of-mass MSD
    msd_time_averaged::Bool    # compute time-averaged MSD in analysis (requires storing coords)
    export_xyz::Bool           # export XYZ trajectory files
    metrics_format::Symbol     # :jld2 or :csv
    nthreads::Int64
    γ::Float64
    rcut_nf::Float64
    L::Float64

    # Activity distribution
    activity_pattern::Symbol   # :random or :block

    # Initialization method
    init_method::Symbol        # :circle or :fourier
    init_kmax::Int64           # Number of Fourier modes for :fourier init method
    init_adaptive_kmax::Bool   # Auto-scale k_max with system size
    init_thermal_scale::Bool   # Use thermal prefactor sqrt(kT/kbond) in mode amplitudes
    init_rg_calibrate::Bool    # Calibrate to target equilibrium Rg

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
                   dt=0.01, dt_thermal=0.0, n_steps=100_000, thermal_steps=100_000,
                   traj_interval=500, metric_mode=:fixed, metric_interval=0, metric_npoints=1000,
                   msd_com=false, msd_time_averaged=false, export_xyz=false, metrics_format=:jld2,
                   L=0.0, nthreads=0, γ=2.0, rcut_nf=2.0,
                   # Activity distribution
                   activity_pattern=:random,
                   # Initialization method
                   init_method=:fourier,
                   init_kmax=10,
                   init_adaptive_kmax=true,
                   init_thermal_scale=true,
                   init_rg_calibrate=true,
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

    # Default metric_interval to traj_interval if not specified
    if metric_interval <= 0
        metric_interval = traj_interval
    end

    # Auto-detect thread count
    if nthreads <= 0
        nthreads = Threads.nthreads()
    end
    
    # Validate system_type
    if system_type ∉ (:single, :double)
        error("system_type must be :single or :double, got $system_type")
    end
    
    # Validate init_method
    if init_method ∉ (:circle, :fourier)
        error("init_method must be :circle or :fourier, got $init_method")
    end

    Parameters(
        system_type,
        KT, mass, kbond, kangle, factive, dt, dt_thermal, n_steps, thermal_steps,
        traj_interval, metric_mode, metric_interval, metric_npoints,
        msd_com, msd_time_averaged, export_xyz, metrics_format,
        nthreads, γ, rcut_nf, L,
        activity_pattern,
        init_method, init_kmax, init_adaptive_kmax, init_thermal_scale, init_rg_calibrate,
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
        # Shuffle for random pattern, keep as-is for block pattern
        return params.activity_pattern == :random ? shuffle(activity) : activity
    else
        # For double rings, create activity vectors for each ring and combine
        total_1 = params.n_monomers_1
        total_2 = params.n_monomers_2

        activity_1 = vcat(repeat([true], params.n_active_1),
                         repeat([false], total_1 - params.n_active_1))
        activity_2 = vcat(repeat([true], params.n_active_2),
                         repeat([false], total_2 - params.n_active_2))

        if params.activity_pattern == :random
            return vcat(shuffle(activity_1), shuffle(activity_2))
        else
            # Block pattern: active monomers at start of each ring
            return vcat(activity_1, activity_2)
        end
    end
end