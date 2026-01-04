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

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--system` | System type: `single` or `double` | `single` |
| `--n-monomers` | Monomers in single ring | 100 |
| `--n-active` | Active monomers in single ring | 0 |
| `--n-monomers-1` | Monomers in ring 1 (double) | 100 |
| `--n-monomers-2` | Monomers in ring 2 (double) | 100 |
| `--n-active-1` | Active monomers in ring 1 | 0 |
| `--n-active-2` | Active monomers in ring 2 | 0 |

### Physics Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--fact` | Active force magnitude | 1.0 |
| `--kangle` | Angle bending constant | 0.0 |
| `--kbond` | Bond spring constant | 30.0 |
| `--KT` | Temperature | 1.0 |
| `--gamma` | Damping coefficient | 2.0 |

### Simulation Control

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--dt` | Integration timestep | 0.01 |
| `--n-steps` | Active dynamics steps | 1,000,000 |
| `--thermal-steps` | Thermalization steps | 200,000 |
| `--logger-steps` | Logging frequency | 500 |

### Advanced Options

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--L` | Box size (0=auto-calculate) | 0.0 |
| `--rcut-neighbor-finder` | Neighbor cutoff distance | 2.0 |
| `--no-minimize` | Skip energy minimization | `false` |
| `--simid` | Simulation identifier | `""` |

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
