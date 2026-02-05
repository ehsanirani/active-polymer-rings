module PlotUtils

using CairoMakie
using Statistics
using LinearAlgebra
using Random

export plot_msd_loglog, plot_msd_linear, plot_msd_comparison
export set_publication_theme!
export compute_msd, compute_msd_com_timeaveraged

"""
    set_publication_theme!()

Set CairoMakie theme for publication-quality figures.
"""
function set_publication_theme!()
    set_theme!(
        fontsize = 16,
        linewidth = 2,
        markersize = 8,
        Axis = (
            xlabelsize = 18,
            ylabelsize = 18,
            titlesize = 20,
            xticklabelsize = 14,
            yticklabelsize = 14,
            xgridvisible = true,
            ygridvisible = true,
            xminorticksvisible = true,
            yminorticksvisible = true,
        ),
        Legend = (
            labelsize = 14,
            framevisible = true,
        )
    )
end

"""
    plot_msd_loglog(lag_times, msd_data;
                    label="MSD",
                    title="Mean-Squared Displacement",
                    add_diffusion_line=false,
                    diffusion_coeff=nothing)

Create log-log plot of MSD data.

# Arguments
- `lag_times`: Array of lag times
- `msd_data`: Array of MSD values
- `label`: Legend label
- `title`: Plot title
- `add_diffusion_line`: If true, adds reference line with slope 1 (diffusive behavior)
- `diffusion_coeff`: If provided, plots MSD = 6D*t line

# Returns
- `fig, ax`: Figure and axis objects
"""
function plot_msd_loglog(lag_times, msd_data;
                         label="MSD",
                         title="Mean-Squared Displacement (log-log)",
                         add_diffusion_line=false,
                         diffusion_coeff=nothing,
                         color=:blue)
    fig = Figure(resolution=(800, 600))
    ax = Axis(fig[1, 1],
              xlabel="Lag time",
              ylabel="MSD",
              title=title,
              xscale=log10,
              yscale=log10)

    # Plot MSD data
    lines!(ax, lag_times, msd_data, label=label, color=color, linewidth=2)

    # Add reference line for diffusive behavior (slope = 1)
    if add_diffusion_line
        # Find a nice range for the reference line
        t_ref = [minimum(lag_times), maximum(lag_times)]
        # Normalize to pass through the data at mid-point
        mid_idx = length(lag_times) ÷ 2
        msd_mid = msd_data[mid_idx]
        t_mid = lag_times[mid_idx]
        msd_ref = msd_mid .* (t_ref ./ t_mid)
        lines!(ax, t_ref, msd_ref, label="slope = 1", color=:black,
               linestyle=:dash, linewidth=1.5)
    end

    # Add theoretical diffusion line if coefficient provided
    if diffusion_coeff !== nothing
        t_theory = range(minimum(lag_times), maximum(lag_times), length=100)
        msd_theory = 6 * diffusion_coeff .* t_theory
        lines!(ax, t_theory, msd_theory, label="6D⋅t (D=$diffusion_coeff)",
               color=:red, linestyle=:dot, linewidth=2)
    end

    axislegend(ax, position=:lt)

    return fig, ax
end

"""
    plot_msd_linear(lag_times, msd_data;
                    label="MSD",
                    title="Mean-Squared Displacement")

Create linear plot of MSD data.

# Returns
- `fig, ax`: Figure and axis objects
"""
function plot_msd_linear(lag_times, msd_data;
                        label="MSD",
                        title="Mean-Squared Displacement",
                        color=:blue)
    fig = Figure(resolution=(800, 600))
    ax = Axis(fig[1, 1],
              xlabel="Lag time",
              ylabel="MSD",
              title=title)

    lines!(ax, lag_times, msd_data, label=label, color=color, linewidth=2)
    axislegend(ax, position=:lt)

    return fig, ax
end

"""
    plot_msd_comparison(lag_times_list, msd_data_list, labels;
                       title="MSD Comparison",
                       loglog=true)

Plot multiple MSD curves for comparison.

# Arguments
- `lag_times_list`: Vector of lag time arrays
- `msd_data_list`: Vector of MSD data arrays
- `labels`: Vector of labels for each curve
- `title`: Plot title
- `loglog`: If true, use log-log scale

# Returns
- `fig, ax`: Figure and axis objects
"""
function plot_msd_comparison(lag_times_list, msd_data_list, labels;
                            title="MSD Comparison",
                            loglog=true)
    fig = Figure(resolution=(1000, 700))

    if loglog
        ax = Axis(fig[1, 1],
                  xlabel="Lag time",
                  ylabel="MSD",
                  title=title,
                  xscale=log10,
                  yscale=log10)
    else
        ax = Axis(fig[1, 1],
                  xlabel="Lag time",
                  ylabel="MSD",
                  title=title)
    end

    colors = [:blue, :red, :green, :orange, :purple, :cyan, :magenta]

    for (i, (lag_times, msd_data, label)) in enumerate(zip(lag_times_list, msd_data_list, labels))
        color = colors[mod1(i, length(colors))]
        lines!(ax, lag_times, msd_data, label=label, color=color, linewidth=2)
    end

    axislegend(ax, position=:lt)

    return fig, ax
end

"""
    plot_msd_all_types(lag_times_noavg, msd_monomer_noavg, msd_com_noavg,
                      lag_times_avg, msd_monomer_avg, msd_com_avg;
                      title="MSD - All Types")

Create comprehensive plot showing all 4 MSD types.

# Returns
- `fig`: Figure with 2 subplots (linear and log-log)
"""
function plot_msd_all_types(lag_times_noavg, msd_monomer_noavg, msd_com_noavg,
                           lag_times_avg, msd_monomer_avg, msd_com_avg;
                           title="Mean-Squared Displacement")
    fig = Figure(resolution=(1400, 600))

    # Linear plot
    ax1 = Axis(fig[1, 1],
               xlabel="Lag time",
               ylabel="MSD",
               title="$title (Linear)")

    lines!(ax1, lag_times_noavg, msd_monomer_noavg,
           label="Monomers (single t₀)", color=:blue, linewidth=2)
    lines!(ax1, lag_times_noavg, msd_com_noavg,
           label="COM (single t₀)", color=:red, linewidth=2)
    lines!(ax1, lag_times_avg, msd_monomer_avg,
           label="Monomers (time-avg)", color=:blue, linewidth=2, linestyle=:dash)
    lines!(ax1, lag_times_avg, msd_com_avg,
           label="COM (time-avg)", color=:red, linewidth=2, linestyle=:dash)

    axislegend(ax1, position=:lt)

    # Log-log plot
    ax2 = Axis(fig[1, 2],
               xlabel="Lag time",
               ylabel="MSD",
               title="$title (Log-Log)",
               xscale=log10,
               yscale=log10)

    lines!(ax2, lag_times_noavg, msd_monomer_noavg,
           label="Monomers (single t₀)", color=:blue, linewidth=2)
    lines!(ax2, lag_times_noavg, msd_com_noavg,
           label="COM (single t₀)", color=:red, linewidth=2)
    lines!(ax2, lag_times_avg, msd_monomer_avg,
           label="Monomers (time-avg)", color=:blue, linewidth=2, linestyle=:dash)
    lines!(ax2, lag_times_avg, msd_com_avg,
           label="COM (time-avg)", color=:red, linewidth=2, linestyle=:dash)

    # Add slope = 1 reference line
    t_ref = [minimum(lag_times_noavg), maximum(lag_times_noavg)]
    mid_idx = length(msd_monomer_noavg) ÷ 2
    msd_mid = msd_monomer_noavg[mid_idx]
    t_mid = lag_times_noavg[mid_idx]
    msd_ref = msd_mid .* (t_ref ./ t_mid)
    lines!(ax2, t_ref, msd_ref, label="slope = 1", color=:black,
           linestyle=:dot, linewidth=1.5)

    axislegend(ax2, position=:lt)

    return fig
end

"""
    compute_msd(coords_history, max_samples=1000)

Calculate time-averaged mean squared displacement (MSD) as a function of time lag.

Uses random sampling of frame pairs to improve performance for long trajectories.
"""
function compute_msd(coords_history, max_samples::Int=1000)
    n_frames = length(coords_history)
    N = length(coords_history[1])

    msd_array = Vector{Float64}(undef, n_frames-1)

    for t in 1:n_frames-1
        total_msd = 0.0

        # Sample frame pairs to avoid O(n_frames^2) complexity
        m = min(max_samples, n_frames - t)
        idxs = randperm(n_frames - t)[1:m]

        for j in idxs
            diff = coords_history[j+t] .- coords_history[j]
            total_msd += sum(d -> dot(d, d), diff)
        end

        msd_array[t] = total_msd / (N * m)
    end

    return msd_array
end

"""
    compute_msd_com_timeaveraged(coords_history, max_samples=1000)

Compute time-averaged MSD for the center of mass.

Uses multiple time origins to improve statistics.
"""
function compute_msd_com_timeaveraged(coords_history, max_samples::Int=1000)
    n_frames = length(coords_history)

    # Compute center of mass for each frame
    com_history = [sum(coords) / length(coords) for coords in coords_history]

    msd_array = Vector{Float64}(undef, n_frames-1)

    for t in 1:n_frames-1
        total_msd = 0.0

        # Sample time origins
        m = min(max_samples, n_frames - t)
        idxs = randperm(n_frames - t)[1:m]

        for j in idxs
            diff = com_history[j+t] - com_history[j]
            total_msd += dot(diff, diff)
        end

        msd_array[t] = total_msd / m
    end

    return msd_array
end

end # module
