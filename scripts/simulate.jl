#!/usr/bin/env julia


using ActiveRings
using JLD2
using ArgParse
using Molly
using ProgressMeter
using CSV
using DataFrames
using Printf

function parse_commandline()
    s = ArgParseSettings(description="Run active polymer ring simulations")
    
    @add_arg_table s begin
        "--system"
            help = "System type: single or double ring"
            arg_type = String
            default = "single"
            range_tester = x -> x in ["single", "double"]
        
        "--n-monomers"
            help = "For single ring: number of monomers"
            arg_type = Int
            default = 100
            dest_name = "n_monomers"

        "--n-active"
            help = "For single ring: number of active monomers"
            arg_type = Int
            default = 0
            dest_name = "n_active"

        "--n-monomers-1"
            help = "For double ring: monomers in ring 1"
            arg_type = Int
            default = 100
            dest_name = "n_monomers_1"

        "--n-monomers-2"
            help = "For double ring: monomers in ring 2"
            arg_type = Int
            default = 100
            dest_name = "n_monomers_2"

        "--n-active-1"
            help = "For double ring: active monomers in ring 1"
            arg_type = Int
            default = 0
            dest_name = "n_active_1"

        "--n-active-2"
            help = "For double ring: active monomers in ring 2"
            arg_type = Int
            default = 0
            dest_name = "n_active_2"
            
        "--fact"
            help = "Magnitude of active force"
            arg_type = Float64
            default = 1.0
            
        "--kangle"
            help = "Angle potential constant"
            arg_type = Float64
            default = 0.0
            
        "--kbond"
            help = "Bond spring constant"
            arg_type = Float64
            default = 30.0
            
        "--dt"
            help = "Time step for active dynamics"
            arg_type = Float64
            default = 0.01

        "--dt-thermal"
            help = "Time step for thermalization (defaults to --dt if not specified)"
            arg_type = Float64
            default = 0.0
            dest_name = "dt_thermal"

        "--n-steps"
            help = "Active simulation steps"
            arg_type = Int
            default = 1_000_000
            dest_name = "n_steps"

        "--thermal-steps"
            help = "Thermalization steps"
            arg_type = Int
            default = 200_000
            dest_name = "thermal_steps"

        "--traj-interval"
            help = "Trajectory logger interval (coords, tangents)"
            arg_type = Int
            default = 500
            dest_name = "traj_interval"

        "--metric-mode"
            help = "Logging mode for metric loggers (MSD, Rg): fixed or logspaced"
            default = "fixed"
            dest_name = "metric_mode"

        "--metric-interval"
            help = "Fixed interval for metric loggers; 0 = same as --traj-interval"
            arg_type = Int
            default = 0
            dest_name = "metric_interval"

        "--metric-npoints"
            help = "Number of sampling points for logspaced metric logging"
            arg_type = Int
            default = 1000
            dest_name = "metric_npoints"

        "--msd-com"
            help = "Enable center-of-mass MSD computation (disabled by default)"
            action = :store_true
            dest_name = "msd_com"

        "--msd-time-averaged"
            help = "Enable time-averaged MSD computation in analysis (stores coords in JLD2)"
            action = :store_true
            dest_name = "msd_time_averaged"

        "--export-xyz"
            help = "Export XYZ trajectory files (disabled by default)"
            action = :store_true
            dest_name = "export_xyz"

        "--metrics-format"
            help = "Output format for metrics (MSD, Rg): jld2 or csv"
            arg_type = String
            default = "jld2"
            range_tester = x -> x in ["jld2", "csv"]
            dest_name = "metrics_format"

        "--L"
            help = "Box size (0 to auto-calculate)"
            arg_type = Float64
            default = 0.0
            
        "--gamma"
            help = "Damping coefficient"
            arg_type = Float64
            default = 2.0
            dest_name = "γ"
            
        "--KT"
            help = "Temperature"
            arg_type = Float64
            default = 1.0
            
        "--rcut-neighbor-finder"
            help = "Neighbor finder cutoff"
            arg_type = Float64
            default = 2.0
            dest_name = "rcut_nf"
            
        "--nthreads"
            help = "Number of threads"
            arg_type = Int
            default = 1
            
        "--simid"
            help = "Simulation identifier"
            arg_type = String
            default = ""
            
        "--no-minimize"
            help = "Skip energy minimization"
            action = :store_true
            dest_name = "no_minimize"
    end
    
    return parse_args(s; as_symbols=true)
end

function init_energy_minimization(params::Parameters, sim_bodies::SimBodies; rcut_nf::Float64=2.0, minimize::Bool=true)
    if !minimize
        return sim_bodies, nothing
    end
    
    loggers = Dict("coords" => Molly.CoordinatesLogger(Float64, params.traj_interval; dims=3))
    sim_bodies.neighbor_finder = create_neighbor_finder(sim_bodies.coords, rcut_nf)
    
    # Use soft sphere system for minimization
    sys = create_soft_system(params, sim_bodies, loggers)
    
    minimizer = Molly.SteepestDescentMinimizer(step_size=0.0001, tol=1e-7)
    Molly.simulate!(sys, minimizer)
    
    updated_bodies = SimBodies(
        atoms=sys.atoms,
        coords=sys.coords,
        velocities=sys.velocities,
        boundary=sys.boundary,
        activity=sim_bodies.activity,
        neighbor_finder=sim_bodies.neighbor_finder
    )
    
    return updated_bodies, sys
end

function create_soft_system(params::Parameters, sim_bodies::SimBodies, loggers)
    n = get_n_particles(params)

    # Create bonds based on system type
    if params.system_type == :single
        # Single ring: 1-2, 2-3, ..., (n-1)-n, n-1
        bonds = Molly.InteractionList2Atoms(
            append!(collect(1:n-1), n),
            append!(collect(2:n), 1),
            [Molly.HarmonicBond(r0=1.05, k=100) for _ in 1:n]
        )
    else  # :double
        # Two separate rings
        n1 = params.n_monomers_1
        n2 = params.n_monomers_2
        n_total = n1 + n2

        bonds = Molly.InteractionList2Atoms(
            collect(1:n_total),
            append!(push!(collect(2:n1), 1),
                    push!(collect(n1+2:n_total), n1+1)),
            [Molly.HarmonicBond(r0=1.05, k=100) for _ in 1:n_total]
        )
    end

    sys = Molly.System(
        atoms=sim_bodies.atoms,
        pairwise_inters=(Molly.SoftSphere(use_neighbors=true),),
        specific_inter_lists=(bonds,),
        coords=sim_bodies.coords,
        velocities=sim_bodies.velocities,
        boundary=sim_bodies.boundary,
        neighbor_finder=sim_bodies.neighbor_finder,
        loggers=loggers,
        force_units=Molly.NoUnits,
        energy_units=Molly.NoUnits
    )

    return sys
end

# Custom logger for progress bar
mutable struct ProgressLogger
    n_steps::Int
    progress::Progress
end

function ProgressLogger(n_steps::Int, total_steps::Int, description::String)
    p = Progress(total_steps; desc=description, barlen=50)
    ProgressLogger(n_steps, p)
end

function Molly.log_property!(logger::ProgressLogger, sys, buffers, neighbors=nothing, step_n::Integer=0; n_threads=1, kwargs...)
    if step_n % logger.n_steps == 0
        update!(logger.progress, step_n)
    end
end

function run_thermalization(params::Parameters, sim_bodies::SimBodies; rcut_nf::Float64=2.0)
    n_particles = get_n_particles(params)

    # Create progress bar
    progress_logger = ProgressLogger(params.traj_interval, params.thermal_steps, "Thermalization")

    # Compute metric logging schedule for this phase
    metric_schedule = make_logging_schedule(
        params.metric_mode, params.thermal_steps;
        fixed_interval=params.metric_interval, n_points=params.metric_npoints
    )

    loggers = Dict(
        "coords" => Molly.CoordinatesLogger(Float64, params.traj_interval; dims=3),
        "rg" => RgLogger(metric_schedule, params.system_type,
                          params.system_type == :single ? params.n_monomers : params.n_monomers_1,
                          params.system_type == :single ? 0 : params.n_monomers_2),
        "tangents" => TangentLogger(params.traj_interval, params),
        "msd" => MSDLogger(metric_schedule, params.system_type,
                           params.system_type == :single ? params.n_monomers : params.n_monomers_1,
                           params.system_type == :single ? 0 : params.n_monomers_2,
                           sim_bodies.coords, sim_bodies.boundary;
                           compute_com=params.msd_com),
        "progress" => progress_logger
    )

    sim_bodies.neighbor_finder = create_neighbor_finder(sim_bodies.coords, rcut_nf)
    sys = create_initial_system(params, sim_bodies; loggers=loggers)

    vverlet = Molly.VelocityVerlet(dt=params.dt_thermal, remove_CM_motion=true)
    Molly.simulate!(sys, vverlet, params.thermal_steps, n_threads=params.nthreads)

    finish!(progress_logger.progress)
    println()  # Ensure progress bar output is flushed
    
    updated_bodies = SimBodies(
        atoms=sys.atoms,
        coords=sys.coords,
        velocities=sys.velocities,
        boundary=sys.boundary,
        activity=sim_bodies.activity,
        neighbor_finder=sim_bodies.neighbor_finder
    )
    
    return sys, updated_bodies
end

function run_active_dynamics(params::Parameters, sim_bodies::SimBodies; rcut_nf::Float64=2.0)
    # Create progress bar
    progress_logger = ProgressLogger(params.traj_interval, params.n_steps, "Active Dynamics")

    # Compute metric logging schedule for this phase
    metric_schedule = make_logging_schedule(
        params.metric_mode, params.n_steps;
        fixed_interval=params.metric_interval, n_points=params.metric_npoints
    )

    loggers = Dict(
        "coords" => Molly.CoordinatesLogger(Float64, params.traj_interval; dims=3),
        "rg" => RgLogger(metric_schedule, params.system_type,
                          params.system_type == :single ? params.n_monomers : params.n_monomers_1,
                          params.system_type == :single ? 0 : params.n_monomers_2),
        "tangents" => TangentLogger(params.traj_interval, params),
        "msd" => MSDLogger(metric_schedule, params.system_type,
                           params.system_type == :single ? params.n_monomers : params.n_monomers_1,
                           params.system_type == :single ? 0 : params.n_monomers_2,
                           sim_bodies.coords, sim_bodies.boundary;
                           compute_com=params.msd_com),
        "progress" => progress_logger
    )

    sim_bodies.neighbor_finder = create_neighbor_finder(sim_bodies.coords, rcut_nf)
    sys = create_initial_system(params, sim_bodies; loggers=loggers)

    vverlet = VelocityVerlet(dt=params.dt, remove_CM_motion=true)
    simulate!(sys, vverlet, params.n_steps, n_threads=params.nthreads)

    finish!(progress_logger.progress)
    println()  # Ensure progress bar output is flushed

    return sys
end

function save_results(params::Parameters, thermal_sys, active_sys, sim_bodies; simid::String="")
    system_type = params.system_type

    # Create output directories
    mkpath("_data")
    mkpath("_data/sims")
    mkpath("_data/jld2")
    mkpath("_data/csv")

    # Generate filename with descriptive pattern
    if system_type == :single
        fname = "single_N$(params.n_monomers)_Nact$(params.n_active)_kang$(params.kangle)_fact$(params.factive)"
    else
        fname = "double_N1_$(params.n_monomers_1)_N2_$(params.n_monomers_2)_Nact1_$(params.n_active_1)_Nact2_$(params.n_active_2)_kang$(params.kangle)_fact$(params.factive)"
    end

    if !isempty(simid)
        fname *= "_$(simid)"
    end

    if params.metrics_format == :jld2
        # Save metrics in JLD2 format
        excluded_loggers = ["progress"]
        if !params.msd_time_averaged
            push!(excluded_loggers, "coords")  # coords only needed for time-averaged MSD
        end

        jldopen("_data/jld2/$(fname).jld2", "w"; compress=true) do file
            file["params"] = params
            file["thermal/loggers"] = Dict(k => v for (k, v) in thermal_sys.loggers if k ∉ excluded_loggers)
            file["active/loggers"] = Dict(k => v for (k, v) in active_sys.loggers if k ∉ excluded_loggers)
        end

        println("\n✓ Results saved:")
        println("  JLD2: _data/jld2/$(fname).jld2")
        if params.msd_time_averaged
            println("        (includes coords for time-averaged MSD)")
        end
    else
        # Save metrics in CSV format
        save_metrics_csv(params, active_sys, "_data/csv", fname, :active)
        save_metrics_csv(params, thermal_sys, "_data/csv", fname, :thermal)

        # Save params to JLD2 (for reference)
        jldopen("_data/jld2/$(fname).jld2", "w"; compress=true) do file
            file["params"] = params
        end

        println("\n✓ Results saved:")
        println("  CSV:  _data/csv/$(fname)_*.csv")
        println("  JLD2: _data/jld2/$(fname).jld2 (params only)")
    end

    # Save XYZ files (only if enabled)
    if params.export_xyz
        save_xyz(thermal_sys, sim_bodies, "_data/sims/$(fname)_thermal.xyz")
        save_xyz(active_sys, sim_bodies, "_data/sims/$(fname)_active.xyz")
        println("  XYZ:  _data/sims/$(fname)_{thermal,active}.xyz")
    end
end

"""
Save metrics (MSD, Rg) to CSV files.
"""
function save_metrics_csv(params::Parameters, sys, output_dir::String, fname::String, phase::Symbol)
    phase_str = string(phase)
    dt = phase == :active ? params.dt : params.dt_thermal

    # Format function for values
    format_val(x) = @sprintf("%.6f", x)

    # Save MSD data
    msd_logger = sys.loggers["msd"]
    if hasproperty(msd_logger, :step_indices) && !isempty(msd_logger.step_indices)
        lag_times = msd_logger.step_indices .* dt
    else
        traj_int = params.traj_interval
        lag_times = collect(1:length(msd_logger.msd_monomer)) .* (dt * traj_int)
    end

    # Monomer MSD
    df_msd = DataFrame(
        lag_time = lag_times,
        msd_monomer = [format_val(x) for x in msd_logger.msd_monomer]
    )
    if params.msd_com && !isempty(msd_logger.msd_com)
        df_msd.msd_com = [format_val(x) for x in msd_logger.msd_com]
    end
    CSV.write(joinpath(output_dir, "$(fname)_$(phase_str)_msd.csv"), df_msd)

    # Save Rg data
    rg_logger = sys.loggers["rg"]
    if hasproperty(rg_logger, :step_indices) && !isempty(rg_logger.step_indices)
        rg_times = rg_logger.step_indices .* dt
    else
        traj_int = params.traj_interval
        rg_times = collect(1:length(rg_logger.Rg_total)) .* (dt * traj_int)
    end

    df_rg = DataFrame(
        time = rg_times,
        Rg = [format_val(x) for x in rg_logger.Rg_total],
        Rg1 = [format_val(x) for x in rg_logger.Rg1],
        Rg2 = [format_val(x) for x in rg_logger.Rg2],
        Rg3 = [format_val(x) for x in rg_logger.Rg3]
    )
    CSV.write(joinpath(output_dir, "$(fname)_$(phase_str)_rg.csv"), df_rg)
end

function save_xyz(sys, sim_bodies, filename)
    n = length(sim_bodies.atoms)
    pol_bonds = sys.specific_inter_lists[1]
    L = sim_bodies.boundary.side_lengths[1]  # Cubic box

    open(filename, "w") do io
        for (frame_idx, coords) in enumerate(sys.loggers["coords"].history)
            # Get tangent vectors for this frame
            tangents = sys.loggers["tangents"].tang_vecs[frame_idx]

            # Wrap coordinates into the box (Molly saves unwrapped coordinates)
            wrapped_coords = [Molly.wrap_coords(c, sim_bodies.boundary) for c in coords]

            # Write atom count
            println(io, n)

            # Write comment line with Properties specification and frame metadata
            # Properties: species(1) + pos(3) + tangent(3) + active(1) = 8 columns total
            bond_str = join(["$(pol_bonds.is[i])-$(pol_bonds.js[i])" for i in 1:length(pol_bonds.is)], ",")
            println(io, "Properties=species:S:1:pos:R:3:tangent:R:3:active:I:1 ",
                    "Lattice=\"$L 0 0 0 $L 0 0 0 $L\" ",
                    "bonds=\"$bond_str\" ",
                    "step=$frame_idx")

            # Write atom lines with all properties as columns
            for i in 1:n
                pos = wrapped_coords[i]
                tang = tangents[i]
                is_active = Int(sim_bodies.atoms[i].active)

                # Format: species x y z tang_x tang_y tang_z active
                println(io, "Monomer $(pos[1]) $(pos[2]) $(pos[3]) ",
                        "$(tang[1]) $(tang[2]) $(tang[3]) $is_active")
            end
        end
    end
end

function main(args=ARGS)
    parsed = parse_commandline()
    
    # Determine system type
    system_type = Symbol(parsed[:system])
    
    # Create parameters
    if system_type == :single
        params = Parameters(
            system_type=:single,
            n_monomers=parsed[:n_monomers],
            n_active=parsed[:n_active],
            factive=parsed[:fact],
            kangle=parsed[:kangle],
            kbond=parsed[:kbond],
            dt=parsed[:dt],
            dt_thermal=parsed[:dt_thermal],
            n_steps=parsed[:n_steps],
            thermal_steps=parsed[:thermal_steps],
            traj_interval=parsed[:traj_interval],
            metric_mode=Symbol(parsed[:metric_mode]),
            metric_interval=parsed[:metric_interval],
            metric_npoints=parsed[:metric_npoints],
            msd_com=parsed[:msd_com],
            msd_time_averaged=parsed[:msd_time_averaged],
            export_xyz=parsed[:export_xyz],
            metrics_format=Symbol(parsed[:metrics_format]),
            L=parsed[:L],
            γ=parsed[:γ],
            KT=parsed[:KT],
            rcut_nf=parsed[:rcut_nf],
            nthreads=parsed[:nthreads]
        )
    else
        params = Parameters(
            system_type=:double,
            n_monomers_1=parsed[:n_monomers_1],
            n_monomers_2=parsed[:n_monomers_2],
            n_active_1=parsed[:n_active_1],
            n_active_2=parsed[:n_active_2],
            factive=parsed[:fact],
            kangle=parsed[:kangle],
            kbond=parsed[:kbond],
            dt=parsed[:dt],
            dt_thermal=parsed[:dt_thermal],
            n_steps=parsed[:n_steps],
            thermal_steps=parsed[:thermal_steps],
            traj_interval=parsed[:traj_interval],
            metric_mode=Symbol(parsed[:metric_mode]),
            metric_interval=parsed[:metric_interval],
            metric_npoints=parsed[:metric_npoints],
            msd_com=parsed[:msd_com],
            msd_time_averaged=parsed[:msd_time_averaged],
            export_xyz=parsed[:export_xyz],
            metrics_format=Symbol(parsed[:metrics_format]),
            L=parsed[:L],
            γ=parsed[:γ],
            KT=parsed[:KT],
            rcut_nf=parsed[:rcut_nf],
            nthreads=parsed[:nthreads]
        )
    end
    
    println("\n🚀 Starting $(system_type) ring simulation")
    println("   System: $(system_type == :single ? "1 ring" : "2 rings")")
    println("   Total monomers: $(get_n_particles(params))")
    println("   Active monomers: $(get_n_active(params))")
    println("   Active force: $(params.factive)")
    println("   Temperature: $(params.KT)")
    println("   Angle constant: $(params.kangle)")
    println()
    
    # Create bodies
    sim_bodies = create_bodies(params)
    
    # Energy minimization
    println("⚡ Step 1: Energy minimization")
    sim_bodies, _ = init_energy_minimization(params, sim_bodies; 
                                            rcut_nf=parsed[:rcut_nf], 
                                            minimize=!parsed[:no_minimize])
    println("   ✓ Minimization complete\n")
    
    # Thermalization
    println("🔥 Step 2: Thermalization ($(params.thermal_steps) steps)")
    thermal_sys, sim_bodies = run_thermalization(params, sim_bodies; rcut_nf=parsed[:rcut_nf])
    println("   ✓ Thermalization complete\n")
    println("   Final Rg: $(thermal_sys.loggers["rg"].Rg_total[end])")
    println()
    
    # Active dynamics
    println("⚛️  Step 3: Active dynamics ($(params.n_steps) steps)")
    t₀ = time()
    active_sys = run_active_dynamics(params, sim_bodies; rcut_nf=parsed[:rcut_nf])
    t₁ = time()
    println("   ✓ Active dynamics complete")
    println("   Time: $(round(t₁ - t₀, digits=1)) seconds")
    println("   Final Rg: $(active_sys.loggers["rg"].Rg_total[end])")
    println()
    
    # Save results
    println("💾 Step 4: Saving results")
    save_results(params, thermal_sys, active_sys, sim_bodies; simid=parsed[:simid])
    
    println("\n✅ Simulation complete!")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
