#!/usr/bin/env julia

using JLD2
using ActiveRings
using Plots

function analyze_msd(jld2_file::String; phase::Symbol=:active)
    println("\n" * "="^70)
    println("MSD ANALYSIS")
    println("="^70)
    println("File: $jld2_file")
    println("Phase: $phase\n")

    # Load data
    data = jldopen(jld2_file, "r")
    phase_str = string(phase)
    coords_history = data["$(phase_str)/loggers"]["coords"].history
    msd_logger = data["$(phase_str)/loggers"]["msd"]
    params = data["params"]
    close(data)

    n_frames = length(coords_history)
    dt = phase == :active ? params.dt : params.dt_thermal
    traj_int = hasproperty(params, :traj_interval) ? params.traj_interval : params.logger_steps
    time_interval = dt * traj_int

    println("Number of frames: $n_frames")
    println("Trajectory time interval: $time_interval")
    println()

    # Check flags (with backward-compatible defaults)
    do_com = hasproperty(params, :msd_com) ? params.msd_com : true
    do_timeavg = hasproperty(params, :msd_time_averaged) ? params.msd_time_averaged : true

    # Get non-time-averaged MSD from logger (already computed)
    println("Loading non-time-averaged MSD from simulation...")
    msd_monomer = msd_logger.msd_monomer
    msd_com = msd_logger.msd_com

    # Compute time-averaged MSD (post-processing, if enabled)
    msd_monomer_avg = Float64[]
    msd_com_avg = Float64[]
    if do_timeavg
        println("Computing time-averaged MSD...")
        msd_monomer_avg = compute_msd(coords_history)
        if do_com
            msd_com_avg = compute_msd_com_timeaveraged(coords_history)
        end
    end

    # Time arrays — use step_indices for non-averaged if available
    if hasproperty(msd_logger, :step_indices) && !isempty(msd_logger.step_indices)
        lag_times = msd_logger.step_indices .* dt
    else
        lag_times = collect(1:length(msd_monomer)) .* time_interval
    end
    lag_times_avg = do_timeavg ? collect(1:length(msd_monomer_avg)) .* time_interval : Float64[]

    # Display results
    println("\nResults Summary:")
    println("-"^70)
    println("Non-averaged MSD (monomers) at max lag: $(msd_monomer[end])")
    if do_com && !isempty(msd_com)
        println("Non-averaged MSD (COM) at max lag: $(msd_com[end])")
    end
    if do_timeavg && !isempty(msd_monomer_avg)
        println("Time-averaged MSD (monomers) at max lag: $(msd_monomer_avg[end])")
        if do_com && !isempty(msd_com_avg)
            println("Time-averaged MSD (COM) at max lag: $(msd_com_avg[end])")
        end
    end

    # Plot
    p1 = plot(lag_times, msd_monomer,
              label="Monomers (single t₀)",
              xlabel="Lag time", ylabel="MSD",
              title="Mean-Squared Displacement",
              linewidth=2, legend=:topleft)
    if do_com && !isempty(msd_com)
        plot!(p1, lag_times, msd_com, label="COM (single t₀)", linewidth=2)
    end
    if do_timeavg && !isempty(msd_monomer_avg)
        plot!(p1, lag_times_avg, msd_monomer_avg,
              label="Monomers (time-avg)", linewidth=2, linestyle=:dash)
        if do_com && !isempty(msd_com_avg)
            plot!(p1, lag_times_avg, msd_com_avg,
                  label="COM (time-avg)", linewidth=2, linestyle=:dash)
        end
    end

    # Log-log plot
    p2 = plot(lag_times, msd_monomer,
              label="Monomers (single t₀)",
              xlabel="Lag time", ylabel="MSD",
              title="MSD (log-log)",
              linewidth=2, legend=:topleft,
              xscale=:log10, yscale=:log10)
    if do_com && !isempty(msd_com)
        plot!(p2, lag_times, msd_com, label="COM (single t₀)", linewidth=2)
    end
    if do_timeavg && !isempty(msd_monomer_avg)
        plot!(p2, lag_times_avg, msd_monomer_avg,
              label="Monomers (time-avg)", linewidth=2, linestyle=:dash)
        if do_com && !isempty(msd_com_avg)
            plot!(p2, lag_times_avg, msd_com_avg,
                  label="COM (time-avg)", linewidth=2, linestyle=:dash)
        end
    end

    plot(p1, p2, layout=(1,2), size=(1200, 500))

    # Generate output filename from input filename
    base_name = splitext(basename(jld2_file))[1]
    output_file = "msd_$(base_name).png"

    savefig(output_file)
    println("\nPlot saved to: $output_file")
    println("="^70)

    return (msd_monomer=msd_monomer,
            msd_com=msd_com,
            msd_monomer_avg=msd_monomer_avg,
            msd_com_avg=msd_com_avg,
            lag_times=lag_times,
            lag_times_avg=lag_times_avg)
end

# Command line interface
if length(ARGS) > 0
    phase = length(ARGS) > 1 ? Symbol(ARGS[2]) : :active
    analyze_msd(ARGS[1]; phase=phase)
else
    println("Usage: julia analyze_msd.jl <jld2_file> [phase]")
    println("  phase: :active (default) or :thermal")
end
