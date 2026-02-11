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
- MSD_com_frame_noavg.csv (non-time-averaged COM-frame MSD)
- MSD_monomer_avg.csv (time-averaged monomer MSD)
- MSD_com_avg.csv (time-averaged COM MSD)
- MSD_com_frame_avg.csv (time-averaged COM-frame MSD)
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

    loggers = data["$(phase_str)/loggers"]
    msd_logger = loggers["msd"]
    params = data["params"]
    # Coords may not be present if msd_time_averaged was disabled
    has_coords = haskey(loggers, "coords")
    coords_history = has_coords ? loggers["coords"].history : nothing
    close(data)

    # Get simulation metadata
    dt = phase == :active ? params.dt : params.dt_thermal
    traj_int = hasproperty(params, :traj_interval) ? params.traj_interval : params.logger_steps
    time_interval = dt * traj_int

    # Check flags (with backward-compatible defaults)
    do_com = hasproperty(params, :msd_com) ? params.msd_com : true
    # Time-averaged MSD requires coords - disable if not available
    do_timeavg = hasproperty(params, :msd_time_averaged) ? params.msd_time_averaged : true
    if do_timeavg && !has_coords
        println("Note: Coords not stored in JLD2, skipping time-averaged MSD export")
        do_timeavg = false
    end

    # Extract run name from filename
    run_name = splitext(basename(jld2_file))[1]

    println("Loading MSD data...")

    # Get non-time-averaged MSD from logger
    msd_monomer_noavg = msd_logger.msd_monomer
    msd_com_noavg = msd_logger.msd_com
    msd_com_frame_noavg = hasproperty(msd_logger, :msd_com_frame) ? msd_logger.msd_com_frame : Float64[]

    # Check COM-frame flag (with backward-compatible default)
    do_com_frame = hasproperty(params, :msd_com_frame) ? params.msd_com_frame : false

    # Format function for values
    format_values(arr) = [@sprintf("%.6f", x) for x in arr]

    # Non-averaged lag times
    if hasproperty(msd_logger, :step_indices) && !isempty(msd_logger.step_indices)
        lag_times = msd_logger.step_indices .* dt
    else
        lag_times = collect(1:length(msd_monomer_noavg)) .* time_interval
    end

    # Save non-time-averaged monomer MSD
    println("Saving non-time-averaged monomer MSD...")
    filename = joinpath(output_dir, "$(run_name)_MSD_monomer_noavg.csv")
    df = DataFrame(lag_time = lag_times, msd = format_values(msd_monomer_noavg))
    CSV.write(filename, df)
    println("  ✓ Saved to: $filename")

    # Save non-time-averaged COM MSD (if enabled)
    if do_com && !isempty(msd_com_noavg)
        println("Saving non-time-averaged COM MSD...")
        filename = joinpath(output_dir, "$(run_name)_MSD_com_noavg.csv")
        df = DataFrame(lag_time = lag_times, msd = format_values(msd_com_noavg))
        CSV.write(filename, df)
        println("  ✓ Saved to: $filename")
    end

    # Save non-time-averaged COM-frame MSD (if enabled)
    if do_com_frame && !isempty(msd_com_frame_noavg)
        println("Saving non-time-averaged COM-frame MSD...")
        filename = joinpath(output_dir, "$(run_name)_MSD_com_frame_noavg.csv")
        df = DataFrame(lag_time = lag_times, msd = format_values(msd_com_frame_noavg))
        CSV.write(filename, df)
        println("  ✓ Saved to: $filename")
    end

    # Compute and save time-averaged MSD (if enabled)
    msd_monomer_avg = Float64[]
    msd_com_avg = Float64[]
    msd_com_frame_avg = Float64[]
    if do_timeavg
        println("Computing time-averaged MSD...")
        msd_monomer_avg = compute_msd(coords_history)
        lag_times_avg = collect(1:length(msd_monomer_avg)) .* time_interval

        println("Saving time-averaged monomer MSD...")
        filename = joinpath(output_dir, "$(run_name)_MSD_monomer_avg.csv")
        df = DataFrame(lag_time = lag_times_avg, msd = format_values(msd_monomer_avg))
        CSV.write(filename, df)
        println("  ✓ Saved to: $filename")

        if do_com
            msd_com_avg = compute_msd_com_timeaveraged(coords_history)
            println("Saving time-averaged COM MSD...")
            filename = joinpath(output_dir, "$(run_name)_MSD_com_avg.csv")
            df = DataFrame(lag_time = lag_times_avg, msd = format_values(msd_com_avg))
            CSV.write(filename, df)
            println("  ✓ Saved to: $filename")
        end

        if do_com_frame
            msd_com_frame_avg = compute_msd_com_frame(coords_history)
            println("Saving time-averaged COM-frame MSD...")
            filename = joinpath(output_dir, "$(run_name)_MSD_com_frame_avg.csv")
            df = DataFrame(lag_time = lag_times_avg, msd = format_values(msd_com_frame_avg))
            CSV.write(filename, df)
            println("  ✓ Saved to: $filename")
        end
    end

    println("\n" * "="^70)
    println("EXPORT COMPLETE")
    println("="^70)
    println("Summary:")
    println("  Non-averaged frames: $(length(msd_monomer_noavg))")
    if do_timeavg
        println("  Time-averaged frames: $(length(msd_monomer_avg))")
    end
    println("  Time interval: $time_interval")
    println("  Max lag time: $(lag_times[end])")
    println("  COM MSD: $(do_com ? "enabled" : "disabled")")
    println("  COM-frame MSD: $(do_com_frame ? "enabled" : "disabled")")
    println("  Time-averaged MSD: $(do_timeavg ? "enabled" : "disabled")")
    println()

    return (msd_monomer_noavg=msd_monomer_noavg,
            msd_com_noavg=msd_com_noavg,
            msd_com_frame_noavg=msd_com_frame_noavg,
            msd_monomer_avg=msd_monomer_avg,
            msd_com_avg=msd_com_avg,
            msd_com_frame_avg=msd_com_frame_avg)
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
