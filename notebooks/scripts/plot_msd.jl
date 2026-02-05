#!/usr/bin/env julia

"""
Plot MSD data in log-log scale from JLD2 or CSV files.

Usage:
  julia --project=notebooks plot_msd.jl <input_file> [options]

Examples:
  # From JLD2 file
  julia --project=notebooks scripts/plot_msd.jl ../_data/jld2/single_N100_Nact50_kang0.0_fact5.0.jld2

  # From CSV file
  julia --project=notebooks scripts/plot_msd.jl ../_data/csv/MSD_monomer_noavg.csv --csv
"""

using CairoMakie
using JLD2
using CSV
using DataFrames
using Statistics

# Load the plotting utilities
include("../src/PlotUtils.jl")
using .PlotUtils

function load_msd_from_jld2(filepath::String; phase::Symbol=:active)
    """Load MSD data from JLD2 file."""
    println("Loading MSD data from JLD2: $filepath")

    data = jldopen(filepath, "r")
    phase_str = string(phase)

    # Check if MSD logger exists
    if !haskey(data["$(phase_str)/loggers"], "msd")
        close(data)
        error("MSD logger not found. Was this simulation run with MSD tracking?")
    end

    # Load non-time-averaged MSD
    msd_logger = data["$(phase_str)/loggers"]["msd"]
    msd_monomer_noavg = msd_logger.msd_monomer
    msd_com_noavg = msd_logger.msd_com

    # Load coordinates for time-averaged computation
    coords_history = data["$(phase_str)/loggers"]["coords"].history
    params = data["params"]
    close(data)

    # Compute time-averaged MSD
    println("Computing time-averaged MSD...")
    msd_monomer_avg = compute_msd(coords_history)
    msd_com_avg = compute_msd_com_timeaveraged(coords_history)

    # Calculate lag times
    dt = phase == :active ? params.dt : params.dt_thermal
    time_interval = dt * params.logger_steps
    lag_times_noavg = collect(1:length(msd_monomer_noavg)) .* time_interval
    lag_times_avg = collect(1:length(msd_monomer_avg)) .* time_interval

    return (
        lag_times_noavg = lag_times_noavg,
        msd_monomer_noavg = msd_monomer_noavg,
        msd_com_noavg = msd_com_noavg,
        lag_times_avg = lag_times_avg,
        msd_monomer_avg = msd_monomer_avg,
        msd_com_avg = msd_com_avg
    )
end

function load_msd_from_csv(filepath::String)
    """Load MSD data from CSV file."""
    println("Loading MSD data from CSV: $filepath")

    df = CSV.read(filepath, DataFrame)

    # Assume first column is lag_time, rest are data columns
    lag_times = df[:, 1]

    # Get data from second column (first data column)
    msd_data = df[:, 2]

    return lag_times, msd_data, names(df)[2]
end

function main(args=ARGS)
    if length(args) == 0
        println("Usage: julia --project=notebooks scripts/plot_msd.jl <input_file> [options]")
        println()
        println("Options:")
        println("  --csv            Input is CSV file (default: JLD2)")
        println("  --phase <name>   Phase to analyze: active (default) or thermal")
        println("  --output <file>  Output filename (default: msd_plot.png)")
        println("  --title <text>   Custom plot title")
        println()
        return
    end

    input_file = args[1]
    is_csv = "--csv" in args
    phase = :active
    output_file = "msd_plot.png"
    custom_title = nothing

    # Parse optional arguments
    for i in 2:length(args)
        if args[i] == "--phase" && i < length(args)
            phase = Symbol(args[i+1])
        elseif args[i] == "--output" && i < length(args)
            output_file = args[i+1]
        elseif args[i] == "--title" && i < length(args)
            custom_title = args[i+1]
        end
    end

    # Set publication theme
    set_publication_theme!()

    # Load data and create plots
    if is_csv
        lag_times, msd_data, data_label = load_msd_from_csv(input_file)
        title = custom_title !== nothing ? custom_title : "MSD - $data_label"

        println("Creating log-log plot...")
        fig, ax = plot_msd_loglog(lag_times, msd_data;
                                   label=data_label,
                                   title=title,
                                   add_diffusion_line=true)
    else
        # Load from JLD2
        data = load_msd_from_jld2(input_file; phase=phase)

        title = custom_title !== nothing ? custom_title : "Mean-Squared Displacement"

        println("Creating comprehensive MSD plots...")
        fig = plot_msd_all_types(
            data.lag_times_noavg,
            data.msd_monomer_noavg,
            data.msd_com_noavg,
            data.lag_times_avg,
            data.msd_monomer_avg,
            data.msd_com_avg;
            title=title
        )
    end

    # Save figure
    output_path = joinpath(@__DIR__, "../output", output_file)
    println("Saving figure to: $output_path")
    save(output_path, fig, px_per_unit=2)

    println("✓ Plot saved successfully!")

    return fig
end

# Run if called as script
if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
