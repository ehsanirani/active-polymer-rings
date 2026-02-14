# ActiveRings

Molecular dynamics simulations for active polymer ring systems. Supports both single rings and double (catenated) ring configurations using Julia and Molly.jl.

## Quick Start

### Installation

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### Run Your First Simulation

```bash
# Single ring with 100 monomers, 10 active
julia --project=. scripts/simulate.jl --system single \
  --n-monomers 100 --n-active 10 --fact 1.0

# Double ring system
julia --project=. scripts/simulate.jl --system double \
  --n-monomers-1 100 --n-monomers-2 100 \
  --n-active-1 10 --n-active-2 10 --fact 1.0
```

### What Happens
- **Energy minimization** → **Thermalization** → **Active dynamics**
- Output saved to `_data/jld2/` and `_data/sims/`

## Configuration

### Using Config Files

Instead of passing many CLI arguments, you can use a TOML config file:

```bash
# Use a config file for single ring
julia --project=. scripts/simulate.jl --config config/single_ring.toml

# Use a config file for double ring
julia --project=. scripts/simulate.jl --config config/double_ring.toml

# Override specific values from config
julia --project=. scripts/simulate.jl --config config/single_ring.toml --n-active 30

# CLI arguments always override config file values
```

Config files are stored in the `config/` folder:
- `config/single_ring.toml` — Complete config for single ring simulations
- `config/double_ring.toml` — Complete config for double (catenated) ring simulations
- `config/default.toml` — Reference with all available options

### Priority Order
1. CLI arguments (highest priority)
2. Config file values
3. Default values (lowest priority)

## Parameters

### System Configuration

**`--system`** (default: `single`)
Choose simulation type: `single` for one polymer ring or `double` for two catenated (linked) rings.

**`--n-monomers`** (default: 100) — *Single ring only*
Number of monomers in the polymer ring. Typical range: 50-500.
Larger values increase computational cost but allow studying longer-wavelength fluctuations and more realistic polymer behavior.

**`--n-active`** (default: 0) — *Single ring only*
Number of active (self-propelled) monomers. Range: 0 to `--n-monomers`.
- `0`: Purely passive polymer
- `10-50`: Partially active system with mixed dynamics
- Equal to `--n-monomers`: Fully active ring

**`--activity-pattern`** (default: `random`)
Distribution of active monomers along the polymer chain:
- `random`: Active monomers randomly distributed (shuffled)
- `block`: Active monomers form a contiguous block starting from monomer 1

**`--n-monomers-1`** (default: 100) — *Double ring only*
Number of monomers in the first polymer ring. Typical range: 50-500.

**`--n-monomers-2`** (default: 100) — *Double ring only*
Number of monomers in the second polymer ring. Typical range: 50-500.
Can differ from `--n-monomers-1` to study asymmetric catenanes.

**`--n-active-1`** (default: 0) — *Double ring only*
Number of active monomers in the first ring. Range: 0 to `--n-monomers-1`.

**`--n-active-2`** (default: 0) — *Double ring only*
Number of active monomers in the second ring. Range: 0 to `--n-monomers-2`.
Can differ from `--n-active-1` to create heterogeneous activity patterns.

### Physics Parameters

**`--fact`** (default: 1.0) — *Dimensionless*
Active force magnitude applied to active monomers along their local tangent direction.
The ratio `fact/KT` determines whether activity or thermal fluctuations dominate dynamics.

**`--kangle`** (default: 0.0) — *Energy units*
Angle bending rigidity constant. Controls polymer stiffness via a cosine angle potential.

**`--kbond`** (default: 30.0) — *Energy units*
Bond spring constant for FENE (finitely extensible nonlinear elastic) bonds connecting monomers.

**`--KT`** (default: 1.0) — *Energy units*
Thermal energy (temperature). All other energies are measured relative to this value.
Increasing `KT` enhances thermal fluctuations.
Key ratio: `fact/KT` determines the balance between active forcing and thermal noise.

**`--gamma`** (default: 2.0) — *Inverse time units*
Damping coefficient (friction) for Langevin dynamics.
Controls momentum relaxation timescales and affects dynamics at short times.

### Simulation Control

**`--dt`** (default: 0.01) — *Time units*
Integration timestep for the active dynamics phase using the Velocity Verlet algorithm.
Must be small enough for numerical stability (typically `dt < 0.1`). If the simulation crashes or shows instabilities, decrease this value.

**`--dt-thermal`** (default: same as `--dt`) — *Time units*
Integration timestep for the thermalization (equilibration) phase.
Defaults to the same value as `--dt` if not specified. You may want to use a larger timestep during thermalization for faster equilibration since there are no active forces that might require smaller timesteps for stability.
Example: `--dt 0.01 --dt-thermal 0.02` uses a larger timestep during equilibration.

**`--n-steps`** (default: 1,000,000)
Number of timesteps for the active dynamics phase (production run).
Typical range: 10⁵ - 10⁷ steps. This is your main data collection period after equilibration.

**`--thermal-steps`** (default: 200,000)
Number of timesteps for thermalization (equilibration) before turning on active forces.
Allows the system to relax to a typical equilibrium configuration. Typical: 10⁴ - 10⁶ steps.

**`--traj-interval`** (default: 500)
Trajectory logging interval. Coordinates and tangent vectors are saved every N steps.
Lower values = more frequent logging = larger output files. Higher values = less temporal resolution.

### Metric Logging

Metric loggers (MSD, Rg) can run on a separate schedule from trajectory loggers (coordinates, tangents). This is useful for capturing early-time dynamics with logarithmic spacing while keeping trajectory files manageable.

**`--metric-mode`** (default: `fixed`)
Logging mode for metric loggers (MSD, Rg): `fixed` or `logspaced`.
- `fixed`: Metrics logged at regular intervals (same behavior as trajectory loggers)
- `logspaced`: Metrics logged at logarithmically spaced steps, providing dense sampling at early times and sparse sampling at late times

**`--metric-interval`** (default: 0)
Fixed interval for metric loggers. Only used when `--metric-mode fixed`.
- `0`: Use the same interval as `--traj-interval`
- Any positive integer: Use this interval instead

**`--metric-npoints`** (default: 1000)
Number of sampling points for logarithmically spaced metric logging. Only used when `--metric-mode logspaced`.
Higher values give denser temporal sampling.

### Per-Metric Logging

MSD and Rg can use different logging schedules. This is useful when you want logspaced MSD (to capture early-time ballistic regime) but fixed-interval Rg (to enable autocorrelation analysis with uniform sampling).

**MSD-specific options** (override global metric settings for MSD only):

**`--msd-mode`** (default: use `--metric-mode`)
Logging mode for MSD: `fixed` or `logspaced`.

**`--msd-interval`** (default: use `--metric-interval`)
Fixed interval for MSD logging.

**`--msd-npoints`** (default: use `--metric-npoints`)
Number of sampling points for logspaced MSD logging.

**Rg-specific options** (override global metric settings for Rg only):

**`--rg-mode`** (default: use `--metric-mode`)
Logging mode for Rg: `fixed` or `logspaced`.

**`--rg-interval`** (default: use `--metric-interval`)
Fixed interval for Rg logging.

**`--rg-npoints`** (default: use `--metric-npoints`)
Number of sampling points for logspaced Rg logging.

**Example: Logspaced MSD with fixed Rg**
```bash
julia --project=. scripts/simulate.jl --config config/single_ring.toml \
  --msd-mode logspaced --msd-npoints 500 \
  --rg-mode fixed --rg-interval 100
```

These options can also be set in config files under `[msd]` and `[rg]` sections:
```toml
[msd]
mode = "logspaced"
npoints = 500

[rg]
mode = "fixed"
interval = 100
```

### MSD Options

**`--msd-com`**
Enable center-of-mass MSD computation during simulation and analysis.
By default, only monomer MSD is computed. Use this flag to also compute COM MSD.

**`--msd-time-averaged`**
Enable time-averaged MSD computation using multiple time origins for better statistics.
By default, only the non-time-averaged (single t₀) MSD from the simulation logger is used.

With `--metrics-format jld2`: Coordinate trajectories are stored in the JLD2 file for post-processing.
With `--metrics-format csv`: Time-averaged MSD is computed and exported to a separate `*_msd_timeaveraged.csv` file.

### Output Options

**`--metrics-format`** (default: `jld2`)
Output format for metrics (MSD, Rg): `jld2` or `csv`.
- `jld2`: Store metrics in binary JLD2 format (compact, preserves full logger structure)
- `csv`: Export metrics to separate CSV files (human-readable, easy to import in other tools)

When using `csv` format, the following files are created:
- `{params}_active_msd.csv` — MSD data (see columns below)
- `{params}_active_msd_timeaveraged.csv` — Time-averaged MSD (if `--msd-time-averaged` enabled)
- `{params}_active_rg.csv` — Rg data (time, Rg, Rg1, Rg2, Rg3)
- `{params}_thermal_msd.csv` — same for thermal phase
- `{params}_thermal_rg.csv` — same for thermal phase

**MSD CSV columns:**
- `lag_time` — time lag
- `msd_monomer` — monomer MSD (always included)
- `msd_com` — center-of-mass MSD (if `--msd-com` enabled)
- `msd_com_frame` — MSD in COM reference frame (if `--msd-com-frame` enabled)

For double ring systems, per-ring columns are automatically included:
- `msd_monomer_1`, `msd_monomer_2` — monomer MSD per ring
- `msd_com_1`, `msd_com_2` — COM MSD per ring (if `--msd-com` enabled)
- `msd_com_frame_1`, `msd_com_frame_2` — COM-frame MSD per ring (if `--msd-com-frame` enabled)

**`--export-xyz`**
Export XYZ trajectory files for visualization (disabled by default).
Creates `{params}_thermal.xyz` and `{params}_active.xyz` files in `_data/sims/`.
Use this flag when you need trajectories for visualization in OVITO, VMD, or similar tools.

### Advanced Options

**`--init-method`** (default: `fourier`)
Ring initialization method: `circle` or `fourier`.
- `circle`: Start with a perfect circular ring configuration
- `fourier`: Generate ring using low-frequency Fourier modes. Produces smooth, globally curved configurations with O(N) complexity. Guaranteed unknotted for `k_max ≤ 15`.

**`--init-kmax`** (default: 10)
Number of Fourier modes for `fourier` initialization method.
Higher values produce more random/rougher configurations. Keep `k_max ≤ 15` to guarantee unknotted rings.

**`--L`** (default: 0.0) — *Length units*
Periodic boundary box side length (cubic box).
- `0.0`: Auto-calculate based on polymer size using scaling law `L = 3.0 × n^0.5887`
- Manual value: Ensure `L > 2 × Rg` (twice the radius of gyration) to avoid self-interactions

**`--rcut-neighbor-finder`** (default: 2.0) — *Length units*
Cutoff distance for neighbor list construction (optimization for pairwise interactions).
Usually no need to change. Must be larger than the interaction range of soft-sphere potentials.

**`--no-minimize`**
Flag to skip initial energy minimization step.
By default, the system performs steepest descent minimization to remove bad contacts. Only skip if continuing from a pre-equilibrated configuration.

**`--simid`** (default: `""`)
Simulation identifier appended to output filenames.
Useful for running parameter sweeps or multiple replicates: e.g., `--simid run01` creates files like `single_100_10_0.0_1.0_run01.jld2`.

**`--nthreads`** (default: 1)
Number of CPU threads for parallel force calculations.
Limited speedup for small systems. Consider using for `n_monomers > 500`.

### Checkpointing and State Management

ActiveRings supports saving and loading simulation states, enabling:
- Periodic checkpoints during long simulations
- Resuming crashed simulations from the last checkpoint
- Reusing equilibrated states as starting points for multiple runs

**`--checkpoint-interval`** (default: 0)
Save checkpoint every N steps during simulation. Set to 0 to disable checkpointing.
Example: `--checkpoint-interval 100000` saves state every 100k steps.

**`--checkpoint-dir`** (default: `_data/checkpoints`)
Directory for checkpoint files. Checkpoints are saved with descriptive filenames:
- `checkpoint_thermal_step100000.jld2`
- `checkpoint_active_step200000.jld2`

**`--checkpoint-keep`** (default: 2)
Number of checkpoints to keep per phase. Older checkpoints are automatically deleted to save disk space (rolling checkpoints).

**`--save-state`**
Save final simulation state to specified file. Useful for creating reusable equilibrated configurations.
Example: `--save-state "_data/states/equilibrated_N200.jld2"`

**`--load-state`**
Load initial state from a previously saved file. Skips initialization and uses coordinates, velocities, and optionally activity pattern from the saved state.
Example: `--load-state "_data/states/equilibrated_N200.jld2"`

**`--load-activity`** (default: `new`)
How to handle activity pattern when loading state:
- `new`: Generate fresh activity pattern from current params (`--n-active`, `--activity-pattern`)
- `keep`: Use the activity pattern from the saved state

**`--resume`**
Resume from the latest checkpoint in `--checkpoint-dir`. Automatically finds the most recent checkpoint and continues the simulation. Implies `--load-activity keep`.

#### Checkpointing Examples

**Run with periodic checkpoints**
```bash
julia --project=. scripts/simulate.jl \
  --n-monomers 200 --n-active 50 --fact 5.0 \
  --n-steps 3000000 \
  --checkpoint-interval 100000 \
  --checkpoint-dir "_data/checkpoints/run1"
```

**Resume a crashed simulation**
```bash
julia --project=. scripts/simulate.jl \
  --n-monomers 200 --n-active 50 --fact 5.0 \
  --n-steps 3000000 \
  --resume \
  --checkpoint-dir "_data/checkpoints/run1"
```

**Save equilibrated passive state for reuse**
```bash
julia --project=. scripts/simulate.jl \
  --n-monomers 200 --n-active 0 \
  --thermal-steps 2000000 --n-steps 0 \
  --save-state "_data/states/passive_N200.jld2"
```

**Load saved state with different activity**
```bash
julia --project=. scripts/simulate.jl \
  --load-state "_data/states/passive_N200.jld2" \
  --n-monomers 200 --n-active 100 --fact 5.0 \
  --load-activity new \
  --n-steps 3000000
```

### Common Parameter Combinations

**Equilibrium polymer (thermal fluctuations only)**
```bash
julia --project=. scripts/simulate.jl --system single \
  --n-monomers 100 --n-active 0
```

**Weakly active flexible ring**
```bash
julia --project=. scripts/simulate.jl --system single \
  --n-monomers 100 --n-active 20 --fact 1.0 --kangle 0.0
```

**Strongly active semi-flexible ring**
```bash
julia --project=. scripts/simulate.jl --system single \
  --n-monomers 100 --n-active 50 --fact 5.0 --kangle 5.0
```

**Double ring with asymmetric activity**
```bash
julia --project=. scripts/simulate.jl --system double \
  --n-monomers-1 100 --n-monomers-2 150 \
  --n-active-1 30 --n-active-2 10 --fact 3.0
```

**Quick test run (small system, short time)**
```bash
julia --project=. scripts/simulate.jl --system single \
  --n-monomers 50 --n-active 10 --fact 2.0 \
  --n-steps 10000 --thermal-steps 5000 --traj-interval 100
```

**Log-spaced metrics for long simulation**
```bash
julia --project=. scripts/simulate.jl --system single \
  --n-monomers 100 --n-active 50 --fact 5.0 \
  --n-steps 10000000 --traj-interval 5000 \
  --metric-mode logspaced --metric-npoints 2000
```

## Batch Simulations

The `tools/sweep.sh` script provides a unified interface for parameter sweeps with base state support.

### Single Parameter Sweep

Sweep over one parameter while keeping others fixed:

```bash
# Sweep active force
./tools/sweep.sh param --config config/single_ring.toml \
    --sweep --fact "1.0 2.0 3.0 4.0 5.0"

# Override some params and sweep another
./tools/sweep.sh param --config config/single_ring.toml \
    --n-monomers 200 --kangle 5.0 \
    --sweep --fact "1.0 2.0 3.0"

# Run 10 replicates with same parameters (for statistics)
./tools/sweep.sh param --config config/single_ring.toml --runs 10

# Sweep parameter with 5 replicates each (3 values × 5 runs = 15 sims)
./tools/sweep.sh param --config config/single_ring.toml \
    --sweep --fact "1.0 2.0 3.0" --runs 5

# Run 4 simulations in parallel
./tools/sweep.sh param --config config/single_ring.toml \
    --sweep --n-active "10 20 30 40" --parallel 4

# Use range notation (expands to 0 25 50 75 100)
./tools/sweep.sh param --config config/single_ring.toml \
    --sweep --n-active "0:25:100"

# Dry run to preview commands
./tools/sweep.sh param --config config/single_ring.toml \
    --sweep --fact "1.0 2.0 3.0" --dry-run
```

### Grid Parameter Sweep

Sweep over multiple parameters (all combinations):

```bash
# Sweep n-active and fact (3x3 = 9 simulations)
./tools/sweep.sh grid --config config/single_ring.toml \
    --sweep --n-active "10 30 50" --sweep --fact "1.0 3.0 5.0"

# Override some params and sweep others
./tools/sweep.sh grid --config config/single_ring.toml \
    --n-monomers 200 --kangle 5.0 \
    --sweep --n-active "10 30 50" --sweep --fact "1.0 3.0"

# Grid sweep with 5 replicates each (2x3x5 = 30 simulations)
./tools/sweep.sh grid --config config/single_ring.toml \
    --sweep --n-active "10 30" --sweep --fact "1.0 2.0 3.0" --runs 5

# With parallel execution
./tools/sweep.sh grid --config config/single_ring.toml \
    --sweep --n-active "10 20 30" --sweep --fact "1.0 2.0" --parallel 4
```

### Base State Support

Base states allow you to pre-equilibrate a polymer configuration and reuse it across multiple simulations. This is useful when you want all runs to start from the same equilibrated state, saving thermalization time.

```bash
# Create and use a passive (n-active=0) base state
# All simulations will load from the same equilibrated configuration
./tools/sweep.sh param --config config/single_ring.toml \
    --base-state passive --n-monomers 200 \
    --sweep --n-active "0 50 100 150 200" --runs 5

# Use an active base state (fully active equilibration)
./tools/sweep.sh param --config config/single_ring.toml \
    --base-state active --n-monomers 200 \
    --sweep --n-active "50 100 150" --runs 3
```

Base states are automatically created if they don't exist and cached in `_data/base_states/`:
- `passive_N200_thermal2M.jld2` — passive equilibration with 2M thermal steps
- `active_N200_thermal2M.jld2` — active equilibration with 2M thermal steps

### Per-Run Base States

For statistical independence, you may want each replica to start from a different equilibrated configuration. The `--base-state-per-run` option creates separate base states for each run index:

```bash
# Each replica gets its own independent equilibrated state
./tools/sweep.sh param --config config/single_ring.toml \
    --base-state passive --base-state-per-run --n-monomers 200 \
    --sweep --n-active "0 50 100" --runs 3
```

This creates:
- `passive_N200_thermal2M_run0.jld2` — used by all run0 simulations (n-active=0, 50, 100)
- `passive_N200_thermal2M_run1.jld2` — used by all run1 simulations
- `passive_N200_thermal2M_run2.jld2` — used by all run2 simulations

**When to use per-run base states:**
- When you need statistically independent replicas
- When studying ensemble properties that shouldn't share initial conditions
- When running multiple independent trajectories for error estimation

**When to use shared base state (default):**
- When comparing different parameter values from the same starting point
- When initial configuration shouldn't vary between runs
- When you want faster setup (only one equilibration needed)

### Sweep Script Options

| Option | Description |
|--------|-------------|
| `--config FILE` | Config file (optional) |
| `--data-dir DIR` | Base data directory (default: _data) |
| `--sweep --param "values"` | Parameter to sweep (space-separated or range notation) |
| `--param value` | Fixed override (applied to all runs) |
| `--runs N` | Number of replicates per parameter combination (default: 1) |
| `--run-start N` | Starting run index (default: 0) |
| `--parallel N` | Run N simulations in parallel |
| `--dry-run` | Preview commands without executing |
| `--prefix STR` | Prefix for simulation IDs (default: run) |

Any `simulate.jl` option can be passed through to individual simulations (e.g., `--metrics-format csv`, `--kangle 5.0`).

**Tip**: Use `--metrics-format csv` when you plan to aggregate results later (see [Aggregating Results](#aggregating-results)).

**Base state options:**

| Option | Description |
|--------|-------------|
| `--base-state TYPE` | Base state type: `passive`, `active`, or file path |
| `--base-state-per-run` | Create separate base state for each replica |
| `--n-monomers N` | Number of monomers (required with passive/active) |
| `--thermal-steps N` | Thermal steps for base state creation (default: 2000000) |
| `--fact F` | Active force for base state creation (default: 5.0) |

**Range notation:**
- `"0:25:100"` expands to `"0 25 50 75 100"`
- `"1.0:0.5:3.0"` expands to `"1.0 1.5 2.0 2.5 3.0"`

### Aggregating Results

After running parameter sweeps with multiple replicas, use `aggregate_sweep.jl` to compute statistics across runs.

**Important**: Aggregation works on CSV metric files. Make sure to use `--metrics-format csv` when running sweeps:

```bash
# Run sweep with CSV output (required for aggregation)
./tools/sweep.sh param --config config/single_ring.toml \
    --metrics-format csv \
    --sweep --n-active "0 50 100" --runs 5
```

This creates CSV files in `_data/csv/` with naming pattern:
- `single_N{n_monomers}_Nact{n_active}_kangle{kangle}_fact{fact}_{simid}_{phase}_{metric}.csv`
- Example: `single_N200_Nact50_kangle0.0_fact5.0_run0_active_rg.csv`

#### Time-Series Aggregation

Combine replicas at each time point to get mean, std, and sem:

```bash
# Aggregate all Rg files from a sweep
julia --project=. scripts/aggregate_sweep.jl \
    --pattern "_data/csv/single_N200_Nact*_run*_active_rg.csv" \
    --mode timeseries \
    --output-dir _data/aggregated

# Group by parameter value (creates one output file per n_active)
julia --project=. scripts/aggregate_sweep.jl \
    --pattern "_data/csv/single_N200_Nact*_run*_active_rg.csv" \
    --mode timeseries \
    --group-by n_active \
    --output-dir _data/aggregated
```

**Time-series output** (`_aggregated.csv`):
```csv
time,Rg_mean,Rg_std,Rg_sem,n_runs
0.5,8.234,0.312,0.099,10
1.0,8.156,0.298,0.094,10
```

#### Steady-State Summary

Extract tail statistics for each file and aggregate across replicas:

```bash
# Summary table with one row per parameter combination
julia --project=. scripts/aggregate_sweep.jl \
    --pattern "_data/csv/single_N200_Nact*_run*_active_rg.csv" \
    --mode summary \
    --group-by n_active \
    --tail-percent 10 \
    --output _data/aggregated/rg_vs_nactive_summary.csv
```

**Summary output** (`_summary.csv`):
```csv
n_active,n_runs,Rg_mean,Rg_std,Rg_sem
0,10,8.52,0.31,0.10
50,10,8.21,0.28,0.09
100,10,7.84,0.35,0.11
```

#### Aggregation Options

| Option | Description |
|--------|-------------|
| `--pattern` | Glob pattern to match files (required, quote to prevent shell expansion) |
| `--mode` | `timeseries`, `summary`, or `both` (default: `timeseries`) |
| `--group-by` | Parameter(s) to group by (comma-separated, e.g., `n_active,kangle`) |
| `--tail-percent` | Percentage of data for tail statistics (default: 20) |
| `--output` | Output file path (for summary mode without per-group files) |
| `--output-dir` | Output directory (default: `_data/aggregated`) |
| `--interpolate` | Interpolate to common time points if times don't match |
| `--verbose` | Print detailed progress information |

**Note**: Quote the pattern to prevent shell expansion of `*` wildcards:
```bash
--pattern "_data/csv/single_N200_Nact*_run*_active_rg.csv"  # correct
--pattern _data/csv/single_N200_Nact*_run*_active_rg.csv    # may fail
```

#### Complete Workflow Example

Here's a full workflow from sweep to aggregated results:

```bash
# 1. Run a parameter sweep with CSV output
./tools/sweep.sh param --config config/single_ring.toml \
    --n-monomers 200 --metrics-format csv \
    --sweep --n-active "0 50 100 150 200" --runs 5 --parallel 4

# 2. Verify output files were created
ls _data/csv/single_N200_Nact*_run*_active_rg.csv

# 3. Aggregate time-series data (one file per n_active value)
julia --project=. scripts/aggregate_sweep.jl \
    --pattern "_data/csv/single_N200_Nact*_run*_active_rg.csv" \
    --mode timeseries \
    --group-by n_active \
    --output-dir _data/aggregated

# 4. Generate steady-state summary table
julia --project=. scripts/aggregate_sweep.jl \
    --pattern "_data/csv/single_N200_Nact*_run*_active_rg.csv" \
    --mode summary \
    --group-by n_active \
    --tail-percent 20 \
    --output _data/aggregated/rg_vs_nactive.csv
```

Output files:
- `_data/aggregated/single_N200_Nact0_*_aggregated.csv` — time-series for n_active=0
- `_data/aggregated/single_N200_Nact50_*_aggregated.csv` — time-series for n_active=50
- ...
- `_data/aggregated/rg_vs_nactive.csv` — summary table with one row per n_active

#### Legacy Aggregation Script

The older `aggregate.jl` script is still available for simple aggregation:

```bash
julia --project=. scripts/aggregate.jl \
    --pattern "_data/csv/single_100_50_0.0_5.0_*_active_msd.csv" \
    --output "_data/csv/aggregated_msd.csv"
```

## Output Files

Simulations create organized output in `_data/`:

### JLD2 Files (`_data/jld2/`)

Binary data files. Contents depend on `--metrics-format`:
- With `--metrics-format jld2` (default): Parameters + all metric loggers (MSD, Rg, tangents)
- With `--metrics-format csv`: Parameters only (metrics are in CSV files)

Additional data stored in JLD2:
- Coordinate trajectories (only if `--msd-time-averaged` is enabled)

**Filename format**:
- Single: `single_{n_monomers}_{n_active}_{kangle}_{factive}[_simid].jld2`
- Double: `double_{n1}_{n2}_{nact1}_{nact2}_{kangle}_{factive}[_simid].jld2`

### CSV Files (`_data/csv/`) — With `--metrics-format csv`

Human-readable metric files:
- `{params}_active_msd.csv` — MSD data (see [MSD CSV columns](#output-options) for details)
- `{params}_active_rg.csv` — Rg data with columns: time, Rg, Rg1, Rg2, Rg3
- `{params}_thermal_msd.csv`, `{params}_thermal_rg.csv` — same for thermal phase

For double ring systems, MSD files include per-ring columns (`msd_monomer_1`, `msd_monomer_2`, etc.).

### XYZ Trajectories (`_data/sims/`) — Optional

Extended XYZ format files for visualization (only created if `--export-xyz` is enabled):
- `{params}_thermal.xyz` - Thermalization phase
- `{params}_active.xyz` - Active dynamics phase

**Features**:
- Periodic boundary conditions (lattice vectors)
- Tangent vectors (per-particle property)
- Activity flags (0=passive, 1=active)
- Bond connectivity list
- Compatible with: OVITO, VMD, ASE, MDAnalysis

**Note**: Bonds are enforced during simulation but must be inferred from distance (~1.0-1.6) or the bonds property during visualization.

## Analysis

### Quick Start

```julia
using ActiveRings
analyze_simulation("_data/jld2/your_file.jld2"; phase=:active)
```

### Available Metrics

The analysis module automatically computes:
- **Rg, Rg1, Rg2, Rg3**: Radius of gyration and principal components
- **MSD**: Mean squared displacement
- **Rs**: End-to-end distance vs contour distance
- **beta**: Tangent-tangent correlation function

### Options

```julia
analyze_simulation(
    "path/to/file.jld2";
    phase=:active,              # or :thermal
    run_name="my_analysis",     # custom name
    output_dir="_data/csv"      # output directory
)
```

### Output Format

CSV files organized in `_data/csv/`:
```
_data/csv/
├── Rg_ring1.csv, Rg_ring2.csv
├── Rg1_ring1.csv, Rg2_ring1.csv, Rg3_ring1.csv
├── MSD_ring1.csv, MSD_ring2.csv
├── Rs_ring1.csv, Rs_ring2.csv
└── beta_ring1.csv, beta_ring2.csv
```

### Custom Analysis

```julia
using JLD2

data = load("_data/jld2/your_file.jld2")
coords = data["active/loggers"]["coords"].history
rg = data["active/loggers"]["rg"].Rg_total

# Your custom analysis here
```

### Equilibration Check

The `check_equilibration.jl` script verifies that a simulation has reached steady state using autocorrelation analysis of the radius of gyration (Rg).

**Analyze an existing simulation:**
```bash
julia --project=. scripts/check_equilibration.jl _data/jld2/simulation.jld2 --save-acf
```

**Run a new simulation and analyze:**
```bash
julia --project=. scripts/check_equilibration.jl \
    --config config/single_ring.toml \
    --n-monomers 200 \
    --n-steps 2000000 \
    --save-acf
```

**Run only thermal phase (passive polymer):**
```bash
julia --project=. scripts/check_equilibration.jl \
    --config config/single_ring.toml \
    --n-monomers 200 \
    --n-steps 0 \
    --thermal-steps 2000000 \
    --save-state _data/base_states/passive_N200.jld2 \
    --save-acf
```

**Run active phase from a base state:**
```bash
julia --project=. scripts/check_equilibration.jl \
    --config config/single_ring.toml \
    --load-state _data/base_states/passive_N200.jld2 \
    --n-active 50 --fact 5.0 \
    --n-steps 2000000 \
    --save-acf
```

**Output interpretation:**
```
Correlation Time Estimates (τ_corr):
  -------------------------------------------------------
  Method          | Lag [samples] |    τ [steps]
  -------------------------------------------------------
  First zero      |          45.0 |       45000.0
  e-folding (1/e) |          32.0 |       32000.0
  Integrated      |          38.5 |       38500.0
  -------------------------------------------------------

Equilibration Assessment:
  Total steps: 2000000
  Correlation time (τ_integrated): 38500 steps
  Number of correlation times: 51.9 (= total_steps / τ)
  Effective independent samples: 26.0
  Status: WELL EQUILIBRATED (>= 20 τ_corr)
```

- **τ [steps]**: Correlation time in the same units as `--n-steps`
- **Number of correlation times**: Should be ≥10 for equilibration, ≥20 for well-equilibrated
- **Effective independent samples**: Actual number of statistically independent measurements

**CSV output** (with `--save-acf`):
- `acf_<basename>_acf.csv` — ACF values vs lag
- `acf_<basename>_summary.csv` — Correlation times and metrics
- `acf_<basename>_rg.csv` — Rg time series

**Options:**
| Option | Description |
|--------|-------------|
| `--config FILE` | Config file to run simulation |
| `--phase PHASE` | Phase to analyze: `active` or `thermal` (auto-detected) |
| `--save-acf` | Save ACF data to CSV files |
| `--output BASE` | Output filename base (default: `acf_<input>`) |
| `--max-lag N` | Maximum lag for ACF (default: n/2) |
| `--quiet` | Suppress detailed output |

All simulation parameters from `simulate.jl` are supported (e.g., `--n-monomers`, `--n-active`, `--fact`, `--kangle`, etc.).

## Citation

If you use this code in your research, please cite:

[To be added]
