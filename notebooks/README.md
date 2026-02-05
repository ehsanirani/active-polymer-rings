# Notebooks - Visualization and Analysis

This subproject contains notebooks, scripts, and utilities for creating publication-quality figures from simulation data.

## Structure

```
notebooks/
├── Project.toml           # Independent dependencies (CairoMakie, etc.)
├── src/
│   └── PlotUtils.jl      # Plotting utilities module
├── scripts/
│   └── plot_msd.jl       # Batch plotting script
├── jupyter/               # Jupyter notebooks
├── pluto/                # Pluto notebooks
│   └── msd_analysis.jl   # Interactive MSD analysis
└── output/               # Generated figures
```

## Setup

### First Time Setup

1. Navigate to the notebooks directory:
   ```bash
   cd notebooks
   ```

2. Activate and instantiate the project:
   ```bash
   julia --project=. -e 'using Pkg; Pkg.instantiate()'
   ```

This will install CairoMakie and other visualization dependencies **without affecting the main simulation project**.

## Usage

### Option 1: Batch Plotting Script

Plot MSD data from JLD2 files:

```bash
# From notebooks directory
julia --project=. scripts/plot_msd.jl ../_data/jld2/single_N100_Nact50_kang0.0_fact5.0.jld2

# With custom output name
julia --project=. scripts/plot_msd.jl ../_data/jld2/single_N100_Nact50_kang0.0_fact5.0.jld2 --output my_msd.png

# From CSV file
julia --project=. scripts/plot_msd.jl ../_data/csv/MSD_monomer_noavg.csv --csv --output msd_csv.png
```

### Option 2: Pluto Notebook (Interactive)

1. Start Pluto:
   ```bash
   cd notebooks
   julia --project=. -e 'using Pluto; Pluto.run()'
   ```

2. Open `pluto/msd_analysis.jl` in the browser

3. Edit the configuration cell to point to your JLD2 file:
   ```julia
   jld2_file = "../../_data/jld2/single_N100_Nact50_kang0.0_fact5.0.jld2"
   phase = :active
   ```

4. The notebook will automatically update all plots

### Option 3: Jupyter Notebook

1. Start Jupyter:
   ```bash
   cd notebooks
   julia --project=. -e 'using IJulia; notebook(dir=pwd())'
   ```

2. Create a new notebook in the `jupyter/` folder

3. Load the plotting utilities:
   ```julia
   include("../src/PlotUtils.jl")
   using .PlotUtils
   ```

## Plotting Functions

The `PlotUtils` module is completely self-contained and provides:
- MSD computation functions (`compute_msd`, `compute_msd_com_timeaveraged`)
- Publication-quality plotting functions
- No dependency on the main ActiveRings project

Available functions:

### Basic MSD Plots

```julia
using CairoMakie
include("src/PlotUtils.jl")
using .PlotUtils

# Set publication theme
set_publication_theme!()

# Log-log plot
fig, ax = plot_msd_loglog(lag_times, msd_data;
                           label="My MSD",
                           add_diffusion_line=true)

# Linear plot
fig, ax = plot_msd_linear(lag_times, msd_data; label="My MSD")

# Save
save("output/my_figure.png", fig, px_per_unit=2)
```

### Comprehensive MSD Plot

Show all 4 MSD types (monomer/COM × non-averaged/time-averaged):

```julia
fig = plot_msd_all_types(
    lag_times_noavg,
    msd_monomer_noavg,
    msd_com_noavg,
    lag_times_avg,
    msd_monomer_avg,
    msd_com_avg
)

save("output/msd_comprehensive.png", fig, px_per_unit=2)
```

### Compare Multiple Simulations

```julia
lag_times_list = [lag_times1, lag_times2, lag_times3]
msd_data_list = [msd1, msd2, msd3]
labels = ["Simulation 1", "Simulation 2", "Simulation 3"]

fig, ax = plot_msd_comparison(lag_times_list, msd_data_list, labels;
                               title="MSD Comparison",
                               loglog=true)
```

## Loading Data

### From JLD2 Files

All necessary functions are in `PlotUtils.jl` - no need to load the main project:

```julia
using JLD2
include("src/PlotUtils.jl")
using .PlotUtils

# Load MSD from JLD2 file
data = jldopen("../_data/jld2/single_N100_Nact50_kang0.0_fact5.0.jld2", "r")
msd_logger = data["active/loggers"]["msd"]
coords_history = data["active/loggers"]["coords"].history
params = data["params"]
close(data)

# Non-time-averaged MSD (already computed during simulation)
msd_monomer_noavg = msd_logger.msd_monomer
msd_com_noavg = msd_logger.msd_com  # may be empty unless --msd-com was used

# Check flags (backward-compatible with old data files)
do_com = hasproperty(params, :msd_com) ? params.msd_com : true
do_timeavg = hasproperty(params, :msd_time_averaged) ? params.msd_time_averaged : true

# Coords may not be present if --msd-time-averaged was not used
has_coords = haskey(loggers, "coords")
if do_timeavg && !has_coords
    do_timeavg = false
end
coords_history = has_coords ? loggers["coords"].history : nothing

# Compute time-averaged MSD (if enabled and coords available)
msd_monomer_avg = do_timeavg ? compute_msd(coords_history) : Float64[]
msd_com_avg = (do_timeavg && do_com) ? compute_msd_com_timeaveraged(coords_history) : Float64[]

# Calculate lag times — use step_indices for non-averaged (supports logspaced)
dt = params.dt
if hasproperty(msd_logger, :step_indices) && !isempty(msd_logger.step_indices)
    lag_times = msd_logger.step_indices .* dt
else
    traj_int = hasproperty(params, :traj_interval) ? params.traj_interval : params.logger_steps
    lag_times = collect(1:length(msd_monomer_noavg)) .* (dt * traj_int)
end
```

### From CSV Files

```julia
using CSV, DataFrames

df = CSV.read("../_data/csv/MSD_monomer_noavg.csv", DataFrame)
lag_times = df[:, 1]
msd_data = df[:, 2]  # First data column
```

## Dependencies

This subproject has its own `Project.toml` with:
- **CairoMakie**: Publication-quality plotting
- **JLD2**: Load simulation data
- **CSV, DataFrames**: Load CSV data
- **Statistics, LinearAlgebra, Random**: MSD computation and analysis

These dependencies are **completely independent** from the main simulation project and won't affect simulation performance or precompilation time.

The `PlotUtils` module includes all necessary MSD computation functions, so you don't need to load the main ActiveRings project.

## Tips

1. **High-DPI figures**: Use `px_per_unit=2` when saving for high resolution
   ```julia
   save("output/figure.png", fig, px_per_unit=2)
   ```

2. **PDF for publications**: Save as PDF for vector graphics
   ```julia
   save("output/figure.pdf", fig)
   ```

3. **Customize theme**: Modify `set_publication_theme!()` in `src/PlotUtils.jl` for consistent styling

4. **Multiple panels**: CairoMakie makes it easy to create multi-panel figures
   ```julia
   fig = Figure(resolution=(1400, 600))
   ax1 = Axis(fig[1, 1], title="Panel A")
   ax2 = Axis(fig[1, 2], title="Panel B")
   ```
