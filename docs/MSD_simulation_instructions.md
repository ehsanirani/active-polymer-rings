# MSD Simulation Instructions

## Overview

This document provides instructions for running MSD simulations of random active polymer rings.

### Parameters (all simulations)
- **N = 200** monomers
- **n_active** = 0, 50, 100, 150, 200 (sweep)
- **f_act = 5**
- **k_angle = 1**
- **100 runs** per parameter combination

### MSD Variants (all 3 enabled)
1. **Monomer MSD** — default
2. **COM MSD** — center-of-mass motion (`--msd-com`)
3. **COM-frame MSD** — internal fluctuations (`--msd-com-frame`)

### Time Scales
| Scale | dt | n_steps | Total time |
|-------|-----|---------|------------|
| Short | 1e-5 | 1,000,000 | 10 |
| Long | 0.01 | 100,000 | 1,000 |

### Systems
1. **Single ring**
2. **Double symmetric** — both rings have same n_active
3. **Double asymmetric** — ring1 active, ring2 passive (n_active_2 = 0)

---

## Directory Structure

Simulations are organized by system type and time scale:

```
_data/
├── single_ring/
│   ├── short_time/
│   │   ├── jld2/
│   │   └── csv/
│   └── long_time/
│       ├── jld2/
│       └── csv/
├── double_symmetric/
│   ├── short_time/
│   └── long_time/
├── double_asymmetric/
│   ├── short_time/
│   └── long_time/
└── aggregated/
    ├── single_ring/
    ├── double_symmetric/
    └── double_asymmetric/
```

---

## 1. Single Ring Simulations

### Short-time MSD

```bash
./tools/sweep.sh param --config config/single_ring.toml \
    --data-dir _data/single_ring/short_time \
    --n-monomers 200 --kangle 1.0 --fact 5.0 \
    --dt 1e-5 --n-steps 1000000 \
    --msd-mode logspaced --msd-npoints 2000 \
    --msd-com true --msd-com-frame true \
    --metrics-format csv \
    --sweep --n-active "0 50 100 150 200" \
    --runs 100 --parallel 8
```

### Long-time MSD

```bash
./tools/sweep.sh param --config config/single_ring.toml \
    --data-dir _data/single_ring/long_time \
    --n-monomers 200 --kangle 1.0 --fact 5.0 \
    --dt 0.01 --n-steps 100000 \
    --msd-mode logspaced --msd-npoints 2000 \
    --msd-com true --msd-com-frame true \
    --metrics-format csv \
    --sweep --n-active "0 50 100 150 200" \
    --runs 100 --parallel 8
```

---

## 2. Double Symmetric Ring Simulations

Both rings have the same number of active monomers.

### Short-time MSD

```bash
./tools/sweep.sh param --config config/double_ring.toml \
    --data-dir _data/double_symmetric/short_time \
    --n-monomers-1 200 --n-monomers-2 200 \
    --kangle 1.0 --fact 5.0 \
    --dt 1e-5 --n-steps 1000000 \
    --msd-mode logspaced --msd-npoints 2000 \
    --msd-com true --msd-com-frame true \
    --metrics-format csv \
    --sweep --n-active-1 "0 50 100 150 200" \
    --sweep --n-active-2 "0 50 100 150 200" \
    --runs 100 --parallel 8
```

**Note**: This creates a grid (5×5 = 25 parameter combinations). For symmetric rings where `n_active_1 = n_active_2`, you need to use `grid` mode or filter combinations.

**Alternative — using grid mode with matching values:**

Since we only want symmetric cases (`n_active_1 = n_active_2`), we need a workaround. Run separate sweeps for each value:

```bash
for nact in 0 50 100 150 200; do
    ./tools/sweep.sh param --config config/double_ring.toml \
        --data-dir _data/double_symmetric/short_time \
        --n-monomers-1 200 --n-monomers-2 200 \
        --n-active-1 $nact --n-active-2 $nact \
        --kangle 1.0 --fact 5.0 \
        --dt 1e-5 --n-steps 1000000 \
        --msd-mode logspaced --msd-npoints 2000 \
        --msd-com true --msd-com-frame true \
        --metrics-format csv \
        --runs 100 --parallel 8
done
```

### Long-time MSD

```bash
for nact in 0 50 100 150 200; do
    ./tools/sweep.sh param --config config/double_ring.toml \
        --data-dir _data/double_symmetric/long_time \
        --n-monomers-1 200 --n-monomers-2 200 \
        --n-active-1 $nact --n-active-2 $nact \
        --kangle 1.0 --fact 5.0 \
        --dt 0.01 --n-steps 100000 \
        --msd-mode logspaced --msd-npoints 2000 \
        --msd-com true --msd-com-frame true \
        --metrics-format csv \
        --runs 100 --parallel 8
done
```

---

## 3. Double Asymmetric Ring Simulations

Ring 1 is active (n_active_1 varies), Ring 2 is passive (n_active_2 = 0).

### Short-time MSD

```bash
./tools/sweep.sh param --config config/double_ring.toml \
    --data-dir _data/double_asymmetric/short_time \
    --n-monomers-1 200 --n-monomers-2 200 \
    --n-active-2 0 \
    --kangle 1.0 --fact 5.0 \
    --dt 1e-5 --n-steps 1000000 \
    --msd-mode logspaced --msd-npoints 2000 \
    --msd-com true --msd-com-frame true \
    --metrics-format csv \
    --sweep --n-active-1 "0 50 100 150 200" \
    --runs 100 --parallel 8
```

### Long-time MSD

```bash
./tools/sweep.sh param --config config/double_ring.toml \
    --data-dir _data/double_asymmetric/long_time \
    --n-monomers-1 200 --n-monomers-2 200 \
    --n-active-2 0 \
    --kangle 1.0 --fact 5.0 \
    --dt 0.01 --n-steps 100000 \
    --msd-mode logspaced --msd-npoints 2000 \
    --msd-com true --msd-com-frame true \
    --metrics-format csv \
    --sweep --n-active-1 "0 50 100 150 200" \
    --runs 100 --parallel 8
```

---

## Aggregating Results

After simulations complete, aggregate MSD data across runs.

### Single Ring Aggregation

```bash
# Short-time
julia --project=. scripts/aggregate_sweep.jl \
    --pattern "_data/single_ring/short_time/csv/*_active_msd.csv" \
    --mode timeseries \
    --group-by n_active \
    --output-dir _data/aggregated/single_ring/short_time

# Long-time
julia --project=. scripts/aggregate_sweep.jl \
    --pattern "_data/single_ring/long_time/csv/*_active_msd.csv" \
    --mode timeseries \
    --group-by n_active \
    --output-dir _data/aggregated/single_ring/long_time

# Summary tables
julia --project=. scripts/aggregate_sweep.jl \
    --pattern "_data/single_ring/short_time/csv/*_active_msd.csv" \
    --mode summary --group-by n_active --tail-percent 20 \
    --output _data/aggregated/single_ring/msd_short_summary.csv

julia --project=. scripts/aggregate_sweep.jl \
    --pattern "_data/single_ring/long_time/csv/*_active_msd.csv" \
    --mode summary --group-by n_active --tail-percent 20 \
    --output _data/aggregated/single_ring/msd_long_summary.csv
```

### Double Symmetric Aggregation

```bash
# Short-time
julia --project=. scripts/aggregate_sweep.jl \
    --pattern "_data/double_symmetric/short_time/csv/*_active_msd.csv" \
    --mode timeseries \
    --group-by n_active_1 \
    --output-dir _data/aggregated/double_symmetric/short_time

# Long-time
julia --project=. scripts/aggregate_sweep.jl \
    --pattern "_data/double_symmetric/long_time/csv/*_active_msd.csv" \
    --mode timeseries \
    --group-by n_active_1 \
    --output-dir _data/aggregated/double_symmetric/long_time
```

### Double Asymmetric Aggregation

```bash
# Short-time
julia --project=. scripts/aggregate_sweep.jl \
    --pattern "_data/double_asymmetric/short_time/csv/*_active_msd.csv" \
    --mode timeseries \
    --group-by n_active_1 \
    --output-dir _data/aggregated/double_asymmetric/short_time

# Long-time
julia --project=. scripts/aggregate_sweep.jl \
    --pattern "_data/double_asymmetric/long_time/csv/*_active_msd.csv" \
    --mode timeseries \
    --group-by n_active_1 \
    --output-dir _data/aggregated/double_asymmetric/long_time
```

---

## Output File Formats

### MSD CSV Columns

**Single ring:**
- `lag_time` — time lag (τ)
- `msd_monomer` — average monomer MSD
- `msd_com` — center-of-mass MSD (if `--msd-com` enabled)
- `msd_com_frame` — MSD in COM frame (if `--msd-com-frame` enabled)

**Double rings** (all columns above plus per-ring columns):
- `msd_monomer_1`, `msd_monomer_2` — monomer MSD for each ring (automatic)
- `msd_com_1`, `msd_com_2` — COM MSD for each ring (if `--msd-com` enabled)
- `msd_com_frame_1`, `msd_com_frame_2` — COM-frame MSD for each ring (if `--msd-com-frame` enabled)

**Note**: Per-ring monomer MSD is computed automatically for double ring systems. The COM variants require their respective flags.

### Aggregated Output
- `*_aggregated.csv` — time-series with mean, std, sem across runs
- `*_summary.csv` — steady-state statistics

---

## Tips

1. **Parallelization**: Adjust `--parallel N` based on available CPU cores. For 100 runs, `--parallel 8` or higher is recommended.

2. **Monitoring progress**: Check simulation progress with:
   ```bash
   ls _data/single_ring/short_time/csv/ | wc -l
   ```

3. **Resume interrupted sweeps**: The sweep script skips existing files by default. Just re-run the same command to continue.

4. **Dry run**: Add `--dry-run` to preview commands without executing:
   ```bash
   ./tools/sweep.sh param --config config/single_ring.toml \
       --sweep --n-active "0 50 100 150 200" --runs 100 --dry-run
   ```

5. **Log-spaced MSD**: We use `--msd-mode logspaced` to capture early-time ballistic behavior with dense sampling at short times.

---

## Estimated Run Count

| System | Values | Runs | Short-time | Long-time | Total |
|--------|--------|------|------------|-----------|-------|
| Single | 5 | 100 | 500 | 500 | 1,000 |
| Double Symmetric | 5 | 100 | 500 | 500 | 1,000 |
| Double Asymmetric | 5 | 100 | 500 | 500 | 1,000 |
| **Total** | | | | | **3,000** |
