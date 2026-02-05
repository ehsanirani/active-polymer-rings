#!/usr/bin/env julia

using JLD2
using ActiveRings
using CSV
using DataFrames
using Printf

"""
Export MSD data from JLD2 file to CSV files.

Creates separate CSV files for:
- MSD_monomer_noavg.csv (non-time-averaged monomer MSD)
- MSD_com_noavg.csv (non-time-averaged COM MSD)
- MSD_monomer_avg.csv (time-averaged monomer MSD)
- MSD_com_avg.csv (time-averaged COM MSD)
"""
function export_msd_to_csv(jld2_file::String; phase::Symbol=:active, output_dir::String="_data/csv")
    println("\n" * "="^70)
    println("MSD DATA EXPORT TO CSV")
    println("="^70)
    println("File: $jld2_file")
    println("Phase: $phase")
    println("Output directory: $output_dir\n")

    # Create output directory
    mkpath(output_dir)

    # Load data
    data = jldopen(jld2_file, "r")
    phase_str = string(phase)

    # Check if MSD logger exists
    if !haskey(data["$(phase_str)/loggers"], "msd")
        close(data)
        error("MSD logger not found in file. Was this simulation run with MSD tracking enabled?")
    end

    coords_history = data["$(phase_str)/loggers"]["coords"].history
    msd_logger = data["$(phase_str)/loggers"]["msd"]
    params = data["params"]
    close(data)

    # Get simulation metadata
    dt = phase == :active ? params.dt : params.dt_thermal
    traj_int = hasproperty(params, :traj_interval) ? params.traj_interval : params.logger_steps
    time_interval = dt * traj_int

    # Extract run name from filename
    run_name = splitext(basename(jld2_file))[1]

    println("Loading MSD data...")

    # Get non-time-averaged MSD from logger
    msd_monomer_noavg = msd_logger.msd_monomer
    msd_com_noavg = msd_logger.msd_com

    # Compute time-averaged MSD
    println("Computing time-averaged MSD...")
    msd_monomer_avg = compute_msd(coords_history)
    msd_com_avg = compute_msd_com_timeaveraged(coords_history)

    # Format function for values
    format_values(arr) = [@sprintf("%.6f", x) for x in arr]

    # Save non-time-averaged monomer MSD
    println("Saving non-time-averaged monomer MSD...")
    filename = joinpath(output_dir, "$(run_name)_MSD_monomer_noavg.csv")
    # Use step_indices for non-averaged lag times if available
    if hasproperty(msd_logger, :step_indices) && !isempty(msd_logger.step_indices)
        lag_times = msd_logger.step_indices .* dt
    else
        lag_times = collect(1:length(msd_monomer_noavg)) .* time_interval
    end

    df = DataFrame(lag_time = lag_times, msd = format_values(msd_monomer_noavg))
    CSV.write(filename, df)
    println("  ✓ Saved to: $filename")

    # Save non-time-averaged COM MSD
    println("Saving non-time-averaged COM MSD...")
    filename = joinpath(output_dir, "$(run_name)_MSD_com_noavg.csv")

    df = DataFrame(lag_time = lag_times, msd = format_values(msd_com_noavg))
    CSV.write(filename, df)
    println("  ✓ Saved to: $filename")

    # Save time-averaged monomer MSD
    println("Saving time-averaged monomer MSD...")
    filename = joinpath(output_dir, "$(run_name)_MSD_monomer_avg.csv")
    lag_times_avg = collect(1:length(msd_monomer_avg)) .* time_interval

    df = DataFrame(lag_time = lag_times_avg, msd = format_values(msd_monomer_avg))
    CSV.write(filename, df)
    println("  ✓ Saved to: $filename")

    # Save time-averaged COM MSD
    println("Saving time-averaged COM MSD...")
    filename = joinpath(output_dir, "$(run_name)_MSD_com_avg.csv")

    df = DataFrame(lag_time = lag_times_avg, msd = format_values(msd_com_avg))
    CSV.write(filename, df)
    println("  ✓ Saved to: $filename")

    println("\n" * "="^70)
    println("EXPORT COMPLETE")
    println("="^70)
    println("Summary:")
    println("  Non-averaged frames: $(length(msd_monomer_noavg))")
    println("  Time-averaged frames: $(length(msd_monomer_avg))")
    println("  Time interval: $time_interval")
    println("  Max lag time: $(lag_times[end])")
    println()

    return (msd_monomer_noavg=msd_monomer_noavg,
            msd_com_noavg=msd_com_noavg,
            msd_monomer_avg=msd_monomer_avg,
            msd_com_avg=msd_com_avg)
end

# Command line interface
if length(ARGS) > 0
    jld2_file = ARGS[1]
    phase = length(ARGS) > 1 ? Symbol(ARGS[2]) : :active
    output_dir = length(ARGS) > 2 ? ARGS[3] : "_data/csv"

    export_msd_to_csv(jld2_file; phase=phase, output_dir=output_dir)
else
    println("Usage: julia export_msd_to_csv.jl <jld2_file> [phase] [output_dir]")
    println()
    println("Arguments:")
    println("  jld2_file   : Path to JLD2 simulation file")
    println("  phase       : :active (default) or :thermal")
    println("  output_dir  : Output directory (default: _data/csv)")
    println()
    println("Example:")
    println("  julia export_msd_to_csv.jl _data/jld2/single_N100_Nact50_kang0.0_fact5.0.jld2")
    println("  julia export_msd_to_csv.jl _data/jld2/single_N100_Nact50_kang0.0_fact5.0.jld2 active _data/my_csv")
end
