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

### Advanced Options

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

## Output Files

Simulations create organized output in `_data/`:

### JLD2 Files (`_data/jld2/`)

Binary data files containing:
- Parameters
- Coordinate trajectories (thermal and active phases)
- Radius of gyration time series
- Tangent vectors and activity flags

**Filename format**:
- Single: `single_{n_monomers}_{n_active}_{kangle}_{factive}[_simid].jld2`
- Double: `double_{n1}_{n2}_{nact1}_{nact2}_{kangle}_{factive}[_simid].jld2`

### XYZ Trajectories (`_data/sims/`)

Extended XYZ format files for visualization:
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

## Citation

If you use this code in your research, please cite:

[To be added]
