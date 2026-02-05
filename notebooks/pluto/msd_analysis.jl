### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 8f3e2a40-a3d1-11ee-1234-0123456789ab
begin
    using CairoMakie
    using JLD2
    using Statistics

    # Load plotting utilities
    include("../src/PlotUtils.jl")
    using .PlotUtils

    set_publication_theme!()

    md"✓ Packages loaded successfully"
end

# ╔═╡ 9f4e3b51-a3d1-11ee-2345-123456789abc
md"""
# MSD Analysis - Active Polymer Rings

This notebook loads and analyzes Mean-Squared Displacement (MSD) data from simulation JLD2 files.

## Features
- Load MSD data from JLD2 files
- Interactive log-log and linear plots
- Compare different simulations
- Export publication-quality figures
"""

# ╔═╡ af5e4c62-a3d1-11ee-3456-23456789abcd
md"""
## Configuration

Specify the JLD2 file path and analysis phase:
"""

# ╔═╡ bf6e5d73-a3d1-11ee-4567-3456789abcde
begin
    # ===== USER CONFIGURATION =====
    jld2_file = "../../_data/jld2/single_N100_Nact50_kang0.0_fact5.0.jld2"  # Change this path
    phase = :active  # :active or :thermal
    # ==============================
end

# ╔═╡ cf7e6e84-a3d1-11ee-5678-456789abcdef
md"""
## Load Data
"""

# ╔═╡ df8e7f95-a3d1-11ee-6789-56789abcdef0
begin
    if isfile(jld2_file)
        data = jldopen(jld2_file, "r")
        phase_str = string(phase)

        # Load MSD logger data (non-time-averaged, already computed)
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

        md"""
        ✓ **Data loaded successfully**
        - Frames (non-averaged): $(length(msd_monomer_noavg))
        - Frames (time-averaged): $(length(msd_monomer_avg))
        - Time interval: $time_interval
        """
    else
        md"""
        ⚠️ **File not found**: `$jld2_file`

        Please update the `jld2_file` variable above with the correct path.
        """
    end
end

# ╔═╡ ef9e8fa6-a3d1-11ee-789a-6789abcdef01
md"""
## MSD Plots

### Comprehensive view (Linear + Log-Log)
"""

# ╔═╡ ffae9fb7-a3d1-11ee-89ab-789abcdef012
if @isdefined(msd_monomer_noavg)
    plot_msd_all_types(
        lag_times_noavg,
        msd_monomer_noavg,
        msd_com_noavg,
        lag_times_avg,
        msd_monomer_avg,
        msd_com_avg;
        title="Mean-Squared Displacement"
    )
end

# ╔═╡ 0fbeafc8-a3d2-11ee-9abc-89abcdef0123
md"""
## Individual Plots

### Monomer MSD (Log-Log)
"""

# ╔═╡ 1fcebfd9-a3d2-11ee-abcd-9abcdef01234
if @isdefined(msd_monomer_noavg)
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
end

# ╔═╡ 2fdecfea-a3d2-11ee-bcde-abcdef012345
md"""
### Center of Mass MSD (Log-Log)
"""

# ╔═╡ 3feedffb-a3d2-11ee-cdef-bcdef0123456
if @isdefined(msd_com_noavg)
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
end

# ╔═╡ 4ffef00c-a3d2-11ee-def0-cdef01234567
md"""
## Summary Statistics
"""

# ╔═╡ 5f0e011d-a3d2-11ee-ef01-def012345678
if @isdefined(msd_monomer_noavg)
    md"""
    ### MSD Summary

    **Max lag time:** $(maximum(lag_times_noavg))

    **Final MSD values:**
    - Monomer (non-avg):   $(round(msd_monomer_noavg[end], digits=4))
    - Monomer (time-avg):  $(round(msd_monomer_avg[end], digits=4))
    - COM (non-avg):       $(round(msd_com_noavg[end], digits=4))
    - COM (time-avg):      $(round(msd_com_avg[end], digits=4))

    **Ratio MSD_monomer / MSD_COM:**
    - Non-averaged:  $(round(msd_monomer_noavg[end] / msd_com_noavg[end], digits=3))
    - Time-averaged: $(round(msd_monomer_avg[end] / msd_com_avg[end], digits=3))
    """
end

# ╔═╡ 6f1e122e-a3d2-11ee-f012-ef0123456789
md"""
## Export Figures

Run this cell to save the comprehensive plot:
"""

# ╔═╡ 7f2e233f-a3d2-11ee-0123-f01234567890
if @isdefined(msd_monomer_noavg)
    let
        output_path = "../output/msd_comprehensive.png"
        mkpath(dirname(output_path))

        save(output_path,
             plot_msd_all_types(lag_times_noavg, msd_monomer_noavg, msd_com_noavg,
                               lag_times_avg, msd_monomer_avg, msd_com_avg),
             px_per_unit=2)

        md"✓ **Saved:** `$output_path`"
    end
end

# ╔═╡ Cell order:
# ╟─9f4e3b51-a3d1-11ee-2345-123456789abc
# ╠═8f3e2a40-a3d1-11ee-1234-0123456789ab
# ╟─af5e4c62-a3d1-11ee-3456-23456789abcd
# ╠═bf6e5d73-a3d1-11ee-4567-3456789abcde
# ╟─cf7e6e84-a3d1-11ee-5678-456789abcdef
# ╠═df8e7f95-a3d1-11ee-6789-56789abcdef0
# ╟─ef9e8fa6-a3d1-11ee-789a-6789abcdef01
# ╠═ffae9fb7-a3d1-11ee-89ab-789abcdef012
# ╟─0fbeafc8-a3d2-11ee-9abc-89abcdef0123
# ╠═1fcebfd9-a3d2-11ee-abcd-9abcdef01234
# ╟─2fdecfea-a3d2-11ee-bcde-abcdef012345
# ╠═3feedffb-a3d2-11ee-cdef-bcdef0123456
# ╟─4ffef00c-a3d2-11ee-def0-cdef01234567
# ╠═5f0e011d-a3d2-11ee-ef01-def012345678
# ╟─6f1e122e-a3d2-11ee-f012-ef0123456789
# ╠═7f2e233f-a3d2-11ee-0123-f01234567890
