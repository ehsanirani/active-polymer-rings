#!/usr/bin/env julia

using ActiveRings
using JLD2
using ArgParse
using Statistics
using Printf
using CSV
using DataFrames

function parse_commandline()
    s = ArgParseSettings(description="Check equilibration of polymer simulation using autocorrelation analysis")

    @add_arg_table s begin
        "input"
            help = "JLD2 file to analyze (or use --config to run simulation first)"
            arg_type = String
            default = ""

        "--config"
            help = "TOML config file to run simulation (if no input file provided)"
            arg_type = String
            default = ""

        "--phase"
            help = "Simulation phase to analyze: active or thermal"
            arg_type = String
            default = "active"

        "--max-lag"
            help = "Maximum lag for ACF (0 = auto, n/2)"
            arg_type = Int
            default = 0
            dest_name = "max_lag"

        "--save-acf"
            help = "Save ACF data to CSV file"
            action = :store_true
            dest_name = "save_acf"

        "--output"
            help = "Output CSV filename (default: acf_<input_basename>.csv)"
            arg_type = String
            default = ""

        "--quiet"
            help = "Suppress detailed output"
            action = :store_true

        # Simulation parameters (passed through to simulate.jl)
        "--n-monomers"
            help = "Number of monomers"
            arg_type = Int
            default = nothing
            dest_name = "n_monomers"

        "--n-active"
            help = "Number of active monomers"
            arg_type = Int
            default = nothing
            dest_name = "n_active"

        "--fact"
            help = "Active force magnitude"
            arg_type = Float64
            default = nothing

        "--kangle"
            help = "Angle potential constant"
            arg_type = Float64
            default = nothing

        "--kbond"
            help = "Bond spring constant"
            arg_type = Float64
            default = nothing

        "--dt"
            help = "Integration timestep"
            arg_type = Float64
            default = nothing

        "--KT"
            help = "Temperature"
            arg_type = Float64
            default = nothing

        "--gamma"
            help = "Damping coefficient"
            arg_type = Float64
            default = nothing

        "--n-steps"
            help = "Number of active dynamics steps (0 = skip active phase)"
            arg_type = Int
            default = nothing
            dest_name = "n_steps"

        "--thermal-steps"
            help = "Number of thermalization steps (0 = skip thermal phase)"
            arg_type = Int
            default = nothing
            dest_name = "thermal_steps"

        "--traj-interval"
            help = "Trajectory logging interval"
            arg_type = Int
            default = nothing
            dest_name = "traj_interval"

        "--metric-interval"
            help = "Metric logging interval"
            arg_type = Int
            default = nothing
            dest_name = "metric_interval"

        "--activity-pattern"
            help = "Activity pattern: random or block"
            arg_type = String
            default = nothing
            dest_name = "activity_pattern"

        "--save-state"
            help = "Save final simulation state to this file"
            arg_type = String
            default = ""
            dest_name = "save_state"

        "--load-state"
            help = "Load initial state from this file (skips thermalization)"
            arg_type = String
            default = ""
            dest_name = "load_state"

        "--simid"
            help = "Simulation identifier for output files"
            arg_type = String
            default = nothing

        "--data-dir"
            help = "Data directory for simulation output"
            arg_type = String
            default = "_data"
            dest_name = "data_dir"
    end

    return parse_args(s; as_symbols=true)
end

"""
    analyze_equilibration(rg_values, dt_sim, metric_interval; max_lag=0, verbose=true)

Analyze equilibration using autocorrelation of Rg time series.

Arguments:
- rg_values: Rg time series
- dt_sim: simulation timestep (params.dt)
- metric_interval: steps between Rg samples

Returns a NamedTuple with ACF, correlation times, and equilibration assessment.
All times are reported in units of dt (same as n_steps).
"""
function analyze_equilibration(rg_values::Vector{Float64}, dt_sim::Float64, metric_interval::Int;
                               max_lag::Int=0, verbose::Bool=true)
    n = length(rg_values)

    # Total steps covered by this data
    total_steps = n * metric_interval

    if verbose
        println("\nAnalyzing Rg time series...")
        println("  Rg samples: $n")
        println("  Metric interval: $metric_interval steps")
        println("  Total steps: $total_steps (= $n × $metric_interval)")
        println("  dt: $dt_sim")
        println("  Mean Rg: $(@sprintf("%.4f", mean(rg_values)))")
        println("  Std Rg: $(@sprintf("%.4f", std(rg_values)))")
    end

    # Compute ACF
    acf = compute_autocorrelation(rg_values, max_lag)
    actual_max_lag = length(acf) - 1

    if verbose
        println("\nComputing autocorrelation (max_lag = $actual_max_lag samples)...")
    end

    # Estimate correlation times with different methods (in samples)
    τ_first_zero = estimate_correlation_time(acf; method=:first_zero)
    τ_e_folding = estimate_correlation_time(acf; method=:e_folding)
    τ_integrated = estimate_correlation_time(acf; method=:integrated)

    # Convert to steps (same unit as n_steps)
    τ_first_zero_steps = τ_first_zero * metric_interval
    τ_e_folding_steps = τ_e_folding * metric_interval
    τ_integrated_steps = τ_integrated * metric_interval

    # Calculate number of correlation times in simulation
    n_corr_times = total_steps / τ_integrated_steps

    # Effective number of independent samples
    n_effective = n / (2 * τ_integrated)

    if verbose
        println("\nCorrelation Time Estimates (τ_corr):")
        println("  " * "-"^55)
        println("  Method          | Lag [samples] |    τ [steps]")
        println("  " * "-"^55)
        @printf("  First zero      | %13.1f | %13.1f\n", τ_first_zero, τ_first_zero_steps)
        @printf("  e-folding (1/e) | %13.1f | %13.1f\n", τ_e_folding, τ_e_folding_steps)
        @printf("  Integrated      | %13.1f | %13.1f\n", τ_integrated, τ_integrated_steps)
        println("  " * "-"^55)
        println("  (τ [steps] = Lag [samples] × metric_interval)")

        println("\nEquilibration Assessment:")
        println("  Total steps: $total_steps")
        println("  Correlation time (τ_integrated): $(@sprintf("%.0f", τ_integrated_steps)) steps")
        println("  Number of correlation times: $(@sprintf("%.1f", n_corr_times)) (= total_steps / τ)")
        println("  Effective independent samples: $(@sprintf("%.1f", n_effective))")

        if n_corr_times >= 20
            println("  Status: WELL EQUILIBRATED (>= 20 τ_corr)")
        elseif n_corr_times >= 10
            println("  Status: LIKELY EQUILIBRATED (>= 10 τ_corr)")
        elseif n_corr_times >= 5
            println("  Status: MARGINAL (5-10 τ_corr) - consider longer simulation")
        else
            println("  Status: NOT EQUILIBRATED (< 5 τ_corr) - need longer simulation")
        end
    end

    return (
        acf = acf,
        τ_first_zero_samples = τ_first_zero,
        τ_e_folding_samples = τ_e_folding,
        τ_integrated_samples = τ_integrated,
        τ_first_zero_steps = τ_first_zero_steps,
        τ_e_folding_steps = τ_e_folding_steps,
        τ_integrated_steps = τ_integrated_steps,
        n_correlation_times = n_corr_times,
        n_effective_samples = n_effective,
        dt_sim = dt_sim,
        metric_interval = metric_interval,
        total_steps = total_steps,
        n_samples = n,
        mean_rg = mean(rg_values),
        std_rg = std(rg_values)
    )
end

"""
    save_acf_csv(result, output_file; rg_values=nothing)

Save ACF data and summary statistics to CSV files.
Creates two files:
- <output>_acf.csv: lag samples, lag steps, and ACF values
- <output>_summary.csv: correlation times and equilibration metrics
"""
function save_acf_csv(result, output_base::String; rg_values::Union{Vector{Float64},Nothing}=nothing)
    acf = result.acf
    metric_interval = result.metric_interval
    lag_samples = collect(0:length(acf)-1)
    lag_steps = lag_samples .* metric_interval

    # Save ACF data
    acf_file = "$(output_base)_acf.csv"
    df_acf = DataFrame(
        lag_samples = lag_samples,
        lag_steps = lag_steps,
        acf = acf
    )
    CSV.write(acf_file, df_acf)
    println("\nACF data saved to: $acf_file")

    # Save summary statistics
    summary_file = "$(output_base)_summary.csv"
    df_summary = DataFrame(
        metric = [
            "n_samples",
            "metric_interval",
            "total_steps",
            "dt_sim",
            "mean_rg",
            "std_rg",
            "tau_first_zero_samples",
            "tau_first_zero_steps",
            "tau_e_folding_samples",
            "tau_e_folding_steps",
            "tau_integrated_samples",
            "tau_integrated_steps",
            "n_correlation_times",
            "n_effective_samples"
        ],
        value = [
            result.n_samples,
            result.metric_interval,
            result.total_steps,
            result.dt_sim,
            result.mean_rg,
            result.std_rg,
            result.τ_first_zero_samples,
            result.τ_first_zero_steps,
            result.τ_e_folding_samples,
            result.τ_e_folding_steps,
            result.τ_integrated_samples,
            result.τ_integrated_steps,
            result.n_correlation_times,
            result.n_effective_samples
        ],
        unit = [
            "",
            "steps",
            "steps",
            "",
            "length",
            "length",
            "samples",
            "steps",
            "samples",
            "steps",
            "samples",
            "steps",
            "",
            ""
        ]
    )
    CSV.write(summary_file, df_summary)
    println("Summary saved to: $summary_file")

    # Optionally save Rg time series
    if rg_values !== nothing
        rg_file = "$(output_base)_rg.csv"
        steps = collect(1:length(rg_values)) .* metric_interval
        df_rg = DataFrame(
            sample = collect(1:length(rg_values)),
            step = steps,
            rg = rg_values
        )
        CSV.write(rg_file, df_rg)
        println("Rg time series saved to: $rg_file")
    end
end

"""
    load_rg_from_jld2(jld2_file, phase)

Load Rg time series, dt_sim, and metric_interval from a JLD2 file.

Returns: (rg_values, dt_sim, metric_interval, params)
"""
function load_rg_from_jld2(jld2_file::String, phase::Symbol)
    data = jldopen(jld2_file, "r")

    phase_str = string(phase)
    if !haskey(data, "$(phase_str)/loggers")
        close(data)
        error("Phase '$phase' not found in JLD2 file. Available: $(keys(data))")
    end

    loggers = data["$(phase_str)/loggers"]
    rg_logger = loggers["rg"]
    params = data["params"]
    close(data)

    rg_values = rg_logger.Rg_total

    # Determine dt based on phase
    dt_sim = phase == :active ? params.dt : params.dt_thermal
    # Handle dt_thermal = 0 case (means use dt)
    if dt_sim == 0.0
        dt_sim = params.dt
    end

    # Get logging interval in steps
    if hasproperty(rg_logger, :step_indices) && !isempty(rg_logger.step_indices) && length(rg_logger.step_indices) > 1
        # Variable spacing - use average interval
        step_diffs = diff(rg_logger.step_indices)
        metric_interval = round(Int, mean(step_diffs))
    else
        # Fixed interval
        metric_interval = hasproperty(params, :metric_interval) && params.metric_interval > 0 ?
                          params.metric_interval : params.traj_interval
    end

    return rg_values, dt_sim, metric_interval, params
end

"""
    run_simulation_and_analyze(config_path, args)

Run a simulation and then analyze its equilibration.
"""
function run_simulation_and_analyze(config_path::String, args::Dict{Symbol,Any})
    # Build command line arguments for simulate.jl
    sim_args = ["--config", config_path]

    # Pass through all simulation parameters
    param_mappings = [
        (:n_monomers, "--n-monomers"),
        (:n_active, "--n-active"),
        (:fact, "--fact"),
        (:kangle, "--kangle"),
        (:kbond, "--kbond"),
        (:dt, "--dt"),
        (:KT, "--KT"),
        (:gamma, "--gamma"),
        (:n_steps, "--n-steps"),
        (:thermal_steps, "--thermal-steps"),
        (:traj_interval, "--traj-interval"),
        (:metric_interval, "--metric-interval"),
        (:activity_pattern, "--activity-pattern"),
        (:simid, "--simid"),
    ]

    for (key, flag) in param_mappings
        if haskey(args, key) && args[key] !== nothing
            push!(sim_args, flag, string(args[key]))
        end
    end

    # String parameters with empty string as default
    if !isempty(args[:save_state])
        push!(sim_args, "--save-state", args[:save_state])
    end
    if !isempty(args[:load_state])
        push!(sim_args, "--load-state", args[:load_state])
    end

    push!(sim_args, "--data-dir", args[:data_dir])

    # Force JLD2 format so we can access Rg data for analysis
    push!(sim_args, "--metrics-format", "jld2")

    println("="^70)
    println("RUNNING SIMULATION")
    println("="^70)
    println("Config: $config_path")

    # Record time before simulation to find the new file
    start_time = time()

    # Run simulation as subprocess
    script_dir = dirname(@__FILE__)
    simulate_script = joinpath(script_dir, "simulate.jl")

    # Build the full command
    cmd = `$(Base.julia_cmd()) --project=. $simulate_script $sim_args`
    println("Running: julia --project=. scripts/simulate.jl $(join(sim_args, " "))")
    println()

    # Run and check for errors
    run(cmd)

    # Find the output JLD2 file (created after start_time)
    jld2_dir = joinpath(args[:data_dir], "jld2")
    mkpath(jld2_dir)  # Ensure directory exists

    jld2_files = filter(f -> endswith(f, ".jld2"), readdir(jld2_dir, join=true))

    if isempty(jld2_files)
        error("No JLD2 files found in $jld2_dir after simulation")
    end

    # Find files created after we started the simulation
    new_files = filter(f -> mtime(f) >= start_time, jld2_files)

    if isempty(new_files)
        # Fallback to most recent file if timing check fails
        jld2_file = sort(jld2_files, by=mtime, rev=true)[1]
        println("\nWarning: Could not identify new file by timestamp, using most recent")
    else
        # Get the most recently modified new file
        jld2_file = sort(new_files, by=mtime, rev=true)[1]
    end

    println("\nAnalyzing: $jld2_file")

    return jld2_file
end

function main()
    args = parse_commandline()

    verbose = !args[:quiet]

    # Determine phase to analyze
    # Auto-detect based on simulation settings if running a new simulation
    phase_str = args[:phase]
    if !isempty(args[:config])
        # Running a simulation - auto-detect phase
        if args[:n_steps] !== nothing && args[:n_steps] == 0
            # No active phase, analyze thermal
            phase_str = "thermal"
            if verbose
                println("Auto-detected: analyzing thermal phase (--n-steps 0)")
            end
        elseif args[:thermal_steps] !== nothing && args[:thermal_steps] == 0
            # No thermal phase, analyze active
            phase_str = "active"
        elseif !isempty(args[:load_state])
            # Loading from state typically means active phase
            phase_str = "active"
        end
    end
    phase = Symbol(phase_str)

    if verbose
        println("\n" * "="^70)
        println("EQUILIBRATION CHECK (Autocorrelation Analysis)")
        println("="^70)
    end

    # Determine input file
    jld2_file = args[:input]

    if isempty(jld2_file)
        # No input file - need to run simulation
        if isempty(args[:config])
            error("Either provide an input JLD2 file or --config to run a simulation")
        end
        jld2_file = run_simulation_and_analyze(args[:config], args)
    end

    if verbose
        println("File: $jld2_file")
        println("Phase: $phase")
    end

    # Load Rg data
    rg_values, dt_sim, metric_interval, params = load_rg_from_jld2(jld2_file, phase)

    # Analyze equilibration
    result = analyze_equilibration(rg_values, dt_sim, metric_interval;
                                   max_lag=args[:max_lag],
                                   verbose=verbose)

    # Save ACF to CSV if requested
    if args[:save_acf]
        output_base = args[:output]
        if isempty(output_base)
            base_name = splitext(basename(jld2_file))[1]
            output_base = "acf_$(base_name)"
        end
        save_acf_csv(result, output_base; rg_values=rg_values)
    end

    if verbose
        println("\n" * "="^70)
    end

    return result
end

# Run if called directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
