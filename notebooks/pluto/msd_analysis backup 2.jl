### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ Cell order:
# ╟─a1b2c3d4-1234-5678-90ab-cdef12345678
# ╠═b2c3d4e5-2345-6789-01bc-def123456789
# ╠═c3d4e5f6-3456-7890-12cd-ef1234567890
# ╟─d4e5f6g7-4567-8901-23de-f12345678901
# ╠═e5f6g7h8-5678-9012-34ef-123456789012
# ╠═f6g7h8i9-6789-0123-45f1-234567890123
# ╟─g7h8i9j0-7890-1234-56f1-345678901234
# ╠═h8i9j0k1-8901-2345-67f1-456789012345
# ╠═i9j0k1l2-9012-3456-78f1-567890123456

# ╔═╡ a1b2c3d4-1234-5678-90ab-cdef12345678
md"""
# MSD Analysis - Active Polymer Rings

This notebook loads and analyzes Mean-Squared Displacement (MSD) data from simulation JLD2 files.

## Features
- Load MSD data from JLD2 files
- Interactive log-log and linear plots
- Compare different simulations
- Export publication-quality figures
"""

# ╔═╡ b2c3d4e5-2345-6789-01bc-def123456789
begin
    using CairoMakie
    using JLD2
    using Statistics

    # Load the parent project to access ActiveRings
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "../.."))
    using ActiveRings

    # Load plotting utilities
    Pkg.activate(@__DIR__)
    include("../src/PlotUtils.jl")
    using .PlotUtils

    set_publication_theme!()
end

# ╔═╡ c3d4e5f6-3456-7890-12cd-ef1234567890
md"""
## Configuration

Specify the JLD2 file path and analysis phase:
"""

# ╔═╡ d4e5f6g7-4567-8901-23de-f12345678901
begin
    # ===== USER CONFIGURATION =====
    jld2_file = "../../_data/jld2/single_100_50_0.0_5.0.jld2"  # Change this path
    phase = :active  # :active or :thermal
    # ==============================
end

# ╔═╡ e5f6g7h8-5678-9012-34ef-123456789012
md"""
## Load Data
"""

# ╔═╡ f6g7h8i9-6789-0123-45f1-234567890123
begin
    println("Loading data from: $jld2_file")

    data = jldopen(jld2_file, "r")
    phase_str = string(phase)

    # Load MSD logger data
    msd_logger = data["$(phase_str)/loggers"]["msd"]
    msd_monomer_noavg = msd_logger.msd_monomer
    msd_com_noavg = msd_logger.msd_com

    # Load coordinates for time-averaged computation
    coords_history = data["$(phase_str)/loggers"]["coords"].history
    params = data["params"]
    close(data)

    # Compute time-averaged MSD
    msd_monomer_avg = compute_msd(coords_history)
    msd_com_avg = compute_msd_com_timeaveraged(coords_history)

    # Calculate lag times
    dt = phase == :active ? params.dt : params.dt_thermal
    time_interval = dt * params.logger_steps
    lag_times_noavg = collect(1:length(msd_monomer_noavg)) .* time_interval
    lag_times_avg = collect(1:length(msd_monomer_avg)) .* time_interval

    println("✓ Data loaded successfully")
    println("  Frames (non-averaged): $(length(msd_monomer_noavg))")
    println("  Frames (time-averaged): $(length(msd_monomer_avg))")
    println("  Time interval: $time_interval")
end

# ╔═╡ g7h8i9j0-7890-1234-56f1-345678901234
md"""
## MSD Plots

### Comprehensive view (Linear + Log-Log)
"""

# ╔═╡ h8i9j0k1-8901-2345-67f1-456789012345
plot_msd_all_types(
    lag_times_noavg,
    msd_monomer_noavg,
    msd_com_noavg,
    lag_times_avg,
    msd_monomer_avg,
    msd_com_avg;
    title="Mean-Squared Displacement"
)

# ╔═╡ i9j0k1l2-9012-3456-78f1-567890123456
md"""
## Individual Plots

### Monomer MSD (Log-Log)
"""

# ╔═╡ j0k1l2m3-0123-4567-89f1-678901234567
begin
    fig_mono, ax_mono = plot_msd_loglog(
        lag_times_noavg,
        msd_monomer_noavg;
        label="Monomers (single t₀)",
        title="Monomer MSD",
        add_diffusion_line=true,
        color=:blue
    )

    lines!(ax_mono, lag_times_avg, msd_monomer_avg,
           label="Monomers (time-avg)",
           color=:blue,
           linestyle=:dash,
           linewidth=2)

    axislegend(ax_mono, position=:lt)
    fig_mono
end

# ╔═╡ k1l2m3n4-1234-5678-90f1-789012345678
md"""
### Center of Mass MSD (Log-Log)
"""

# ╔═╡ l2m3n4o5-2345-6789-01f1-890123456789
begin
    fig_com, ax_com = plot_msd_loglog(
        lag_times_noavg,
        msd_com_noavg;
        label="COM (single t₀)",
        title="Center of Mass MSD",
        add_diffusion_line=true,
        color=:red
    )

    lines!(ax_com, lag_times_avg, msd_com_avg,
           label="COM (time-avg)",
           color=:red,
           linestyle=:dash,
           linewidth=2)

    axislegend(ax_com, position=:lt)
    fig_com
end

# ╔═╡ m3n4o5p6-3456-7890-12f1-901234567890
md"""
## Export Figures

To save figures, evaluate the following cells:
"""

# ╔═╡ n4o5p6q7-4567-8901-23f1-012345678901
begin
    # Save comprehensive plot
    save("../output/msd_comprehensive.png",
         plot_msd_all_types(lag_times_noavg, msd_monomer_noavg, msd_com_noavg,
                           lag_times_avg, msd_monomer_avg, msd_com_avg),
         px_per_unit=2)

    println("✓ Saved: ../output/msd_comprehensive.png")
end

# ╔═╡ o5p6q7r8-5678-9012-34f1-123456789012
md"""
## Summary Statistics
"""

# ╔═╡ p6q7r8s9-6789-0123-45f1-234567890123
begin
    println("MSD Summary Statistics")
    println("="^60)
    println("Max lag time: $(maximum(lag_times_noavg))")
    println()
    println("Final MSD values:")
    println("  Monomer (non-avg):   $(msd_monomer_noavg[end])")
    println("  Monomer (time-avg):  $(msd_monomer_avg[end])")
    println("  COM (non-avg):       $(msd_com_noavg[end])")
    println("  COM (time-avg):      $(msd_com_avg[end])")
    println()
    println("Ratio MSD_monomer / MSD_COM:")
    println("  Non-averaged:  $(msd_monomer_noavg[end] / msd_com_noavg[end])")
    println("  Time-averaged: $(msd_monomer_avg[end] / msd_com_avg[end])")
end
