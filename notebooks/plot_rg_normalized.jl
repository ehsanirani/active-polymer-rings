#!/usr/bin/env julia
# =============================================================================
# Plot normalized mean Rg vs n_active
# =============================================================================
#
# Plots mean Rg / Rg_passive versus n_active, where:
#   - Mean Rg is computed from the tail percentage of the active phase
#   - Rg_passive is the mean Rg at n_active=0 (passive case)
#
# Supports two normalization modes:
#   - Global (default): All runs normalized by mean passive Rg across all runs
#   - Per-run (--per-run-normalization): Each run normalized by its own passive Rg
#
# Per-run normalization is useful when using --base-state-per-run in sweep.sh,
# where each replica has its own independent equilibrated base state.
#
# Usage:
#   julia --project=. scripts/plot_rg_normalized.jl [options]
#
# Examples:
#   julia --project=. scripts/plot_rg_normalized.jl --pattern "block_*"
#   julia --project=. scripts/plot_rg_normalized.jl --tail-percent 20
#   julia --project=. scripts/plot_rg_normalized.jl --pattern "block_*" --runs
#   julia --project=. scripts/plot_rg_normalized.jl --pattern "block_*" --per-run-normalization
#
# =============================================================================

using CSV
using DataFrames
using CairoMakie
using Statistics
using ArgParse
using Printf
using JLD2

function parse_commandline()
    s = ArgParseSettings(description="Plot normalized mean Rg vs n_active",
                         autofix_names=true)

    @add_arg_table s begin
        "--input-dir"
            help = "Directory containing CSV files"
            arg_type = String
            default = "_data/csv"
        "--output-dir"
            help = "Directory for output plots"
            arg_type = String
            default = "_data/plots"
        "--pattern"
            help = "Filename pattern to match (e.g., block_*)"
            arg_type = String
            default = "block_*"
        "--tail-percent"
            help = "Percentage of data points from end to use for mean calculation"
            arg_type = Float64
            default = 10.0
        "--runs"
            help = "Data has multiple runs per n_active (aggregate them)"
            action = :store_true
        "--per-run-normalization"
            help = "Normalize each run by its own passive Rg (for per-run base states)"
            action = :store_true
        "--export-csv"
            help = "Export aggregated data to CSV file"
            action = :store_true
        "--dpi"
            help = "Resolution for PNG output"
            arg_type = Int
            default = 300
        "--base-states-dir"
            help = "Directory containing base state JLD2 files (fallback for passive Rg)"
            arg_type = String
            default = ""
    end

    return parse_args(s; as_symbols=true)
end

function load_rg_from_base_state(jld2_file::String)
    """Load Rg timeseries from a base state JLD2 file (thermal phase)."""
    try
        data = jldopen(jld2_file, "r") do file
            if haskey(file, "thermal/loggers")
                loggers = file["thermal/loggers"]
                if haskey(loggers, "rg")
                    rg_logger = loggers["rg"]
                    return DataFrame(Rg = rg_logger.Rg_total)
                end
            end
            return nothing
        end
        return data
    catch e
        @warn "Failed to load Rg from $jld2_file: $e"
        return nothing
    end
end

function find_base_state_file(base_states_dir::String, run_id::Int, n_monomers::Int)
    """Find the base state file for a given run."""
    if isempty(base_states_dir) || !isdir(base_states_dir)
        return nothing
    end

    # Try common naming patterns (as regex patterns)
    patterns = [
        Regex("passive_N$(n_monomers)_thermal.*_run$(run_id)\\.jld2"),
        Regex("passive_N$(n_monomers)_.*_run$(run_id)\\.jld2"),
        Regex(".*_run$(run_id)\\.jld2")
    ]

    all_files = readdir(base_states_dir)
    for pattern in patterns
        files = filter(f -> occursin(pattern, f), all_files)
        if !isempty(files)
            return joinpath(base_states_dir, first(files))
        end
    end
    return nothing
end

function extract_nactive(filename::String)
    # Match patterns like "Nact50" or "block_50" or "block_50_run0"
    m = match(r"Nact(\d+)", filename)
    if m !== nothing
        return parse(Int, m.captures[1])
    end
    m = match(r"block_(\d+)(?:_run\d+)?_", filename)
    if m !== nothing
        return parse(Int, m.captures[1])
    end
    return nothing
end

function extract_nmonomers(filename::String)
    m = match(r"_N(\d+)_", filename)
    return m !== nothing ? parse(Int, m.captures[1]) : 200
end

function extract_run_id(filename::String)
    m = match(r"_run(\d+)_", filename)
    return m !== nothing ? parse(Int, m.captures[1]) : 0
end

function match_pattern(filename::String, pattern::String)
    regex_str = replace(pattern, "*" => ".*", "?" => ".")
    return occursin(Regex(regex_str), filename)
end

function load_rg_data(input_dir::String, pattern::String)
    all_files = readdir(input_dir, join=true)

    # Filter for active phase Rg files matching pattern
    files = filter(all_files) do f
        bn = basename(f)
        endswith(bn, "_active_rg.csv") && match_pattern(bn, pattern)
    end

    if isempty(files)
        error("No active phase Rg CSV files found in $input_dir matching pattern '$pattern'")
    end

    # Group by n_active (may have multiple runs per n_active)
    data = Dict{Int, Vector{DataFrame}}()
    n_monomers = 200

    for file in files
        bn = basename(file)
        n_active = extract_nactive(bn)
        if n_active === nothing
            continue
        end
        n_monomers = extract_nmonomers(bn)

        df = CSV.read(file, DataFrame)

        # For n_active=0, if active phase data is empty, try thermal phase
        if n_active == 0 && (nrow(df) == 0 || all(ismissing, df.Rg))
            thermal_file = replace(file, "_active_rg.csv" => "_thermal_rg.csv")
            if isfile(thermal_file)
                df = CSV.read(thermal_file, DataFrame)
                if nrow(df) > 0
                    println("  Loaded n_active=$n_active (thermal): $(basename(thermal_file))")
                    if !haskey(data, n_active)
                        data[n_active] = DataFrame[]
                    end
                    push!(data[n_active], df)
                    continue
                end
            end
            # No thermal file or empty - skip
            println("  Skipped n_active=$n_active: empty active data, no thermal fallback")
            continue
        end

        if !haskey(data, n_active)
            data[n_active] = DataFrame[]
        end
        push!(data[n_active], df)
        println("  Loaded n_active=$n_active: $bn")
    end

    # Fallback: if no n_active=0 data found, look for base state thermal CSV
    if !haskey(data, 0)
        # Look for base_state_creation thermal files (passive reference)
        base_thermal_pattern = Regex(".*_base_state_creation_thermal_rg\\.csv")
        base_thermal_files = filter(f -> occursin(base_thermal_pattern, basename(f)), all_files)

        if !isempty(base_thermal_files)
            for base_file in base_thermal_files
                df = CSV.read(base_file, DataFrame)
                if nrow(df) > 0
                    # Try to extract n_monomers from base state file
                    bn = basename(base_file)
                    m = match(r"_N(\d+)_", bn)
                    if m !== nothing
                        n_monomers = parse(Int, m.captures[1])
                    end

                    println("  Loaded n_active=0 (base state thermal): $bn")
                    if !haskey(data, 0)
                        data[0] = DataFrame[]
                    end
                    push!(data[0], df)
                end
            end
        end
    end

    return data, n_monomers
end

function load_rg_data_by_run(input_dir::String, pattern::String; base_states_dir::String="")
    all_files = readdir(input_dir, join=true)

    # Filter for active phase Rg files matching pattern
    files = filter(all_files) do f
        bn = basename(f)
        endswith(bn, "_active_rg.csv") && match_pattern(bn, pattern)
    end

    if isempty(files)
        error("No active phase Rg CSV files found in $input_dir matching pattern '$pattern'")
    end

    # Group by run_id, then by n_active: Dict{run_id => Dict{n_active => DataFrame}}
    data_by_run = Dict{Int, Dict{Int, DataFrame}}()
    n_monomers = 200

    for file in files
        bn = basename(file)
        n_active = extract_nactive(bn)
        if n_active === nothing
            continue
        end
        run_id = extract_run_id(bn)
        n_monomers = extract_nmonomers(bn)

        df = CSV.read(file, DataFrame)

        # For n_active=0, if active phase data is empty, try fallbacks
        if n_active == 0 && (nrow(df) == 0 || all(ismissing, df.Rg))
            # Try 1: thermal phase CSV (same filename pattern)
            thermal_file = replace(file, "_active_rg.csv" => "_thermal_rg.csv")
            if isfile(thermal_file)
                df = CSV.read(thermal_file, DataFrame)
                if nrow(df) > 0
                    println("  Loaded run$run_id, n_active=$n_active (thermal CSV): $(basename(thermal_file))")
                    if !haskey(data_by_run, run_id)
                        data_by_run[run_id] = Dict{Int, DataFrame}()
                    end
                    data_by_run[run_id][n_active] = df
                    continue
                end
            end

            # Try 2: base_state_creation thermal CSV (per-run base state)
            # Look for files like *_base_state_creation_run{run_id}_thermal_rg.csv
            base_state_thermal_pattern = Regex(".*_base_state_creation_run$(run_id)_thermal_rg\\.csv")
            base_thermal_files = filter(f -> occursin(base_state_thermal_pattern, basename(f)), all_files)
            if !isempty(base_thermal_files)
                df = CSV.read(first(base_thermal_files), DataFrame)
                if nrow(df) > 0
                    println("  Loaded run$run_id, n_active=$n_active (base state thermal): $(basename(first(base_thermal_files)))")
                    if !haskey(data_by_run, run_id)
                        data_by_run[run_id] = Dict{Int, DataFrame}()
                    end
                    data_by_run[run_id][n_active] = df
                    continue
                end
            end

            # Try 3: base state JLD2 file (if it contains thermal loggers - unlikely)
            base_state_file = find_base_state_file(base_states_dir, run_id, n_monomers)
            if base_state_file !== nothing
                df = load_rg_from_base_state(base_state_file)
                if df !== nothing && nrow(df) > 0
                    println("  Loaded run$run_id, n_active=$n_active (base state JLD2): $(basename(base_state_file))")
                    if !haskey(data_by_run, run_id)
                        data_by_run[run_id] = Dict{Int, DataFrame}()
                    end
                    data_by_run[run_id][n_active] = df
                    continue
                end
            end

            # No fallback available
            println("  Skipped run$run_id, n_active=$n_active: empty active data, no fallback found")
            continue
        end

        if !haskey(data_by_run, run_id)
            data_by_run[run_id] = Dict{Int, DataFrame}()
        end
        data_by_run[run_id][n_active] = df
        println("  Loaded run$run_id, n_active=$n_active: $bn")
    end

    return data_by_run, n_monomers
end

function compute_tail_mean(rg_values::Vector{Float64}, tail_percent::Float64)
    n_points = length(rg_values)
    start_idx = max(1, n_points - round(Int, n_points * tail_percent / 100) + 1)
    return mean(rg_values[start_idx:end])
end

function compute_tail_stats(rg_values::Vector{Float64}, tail_percent::Float64)
    n_points = length(rg_values)
    start_idx = max(1, n_points - round(Int, n_points * tail_percent / 100) + 1)
    tail_values = rg_values[start_idx:end]
    return mean(tail_values), std(tail_values)
end

function process_global_normalization(args, tail_percent)
    println("Loading Rg data...")
    data, n_monomers = load_rg_data(args[:input_dir], args[:pattern])
    println("Loaded data for $(length(data)) n_active values")
    println()

    n_active_values = sort(collect(keys(data)))

    # Check that we have passive (n_active=0) data
    if !haskey(data, 0)
        error("No data for n_active=0 (passive case) found. Cannot normalize.")
    end

    # Compute passive Rg (mean over all runs if multiple)
    passive_rg_values = Float64[]
    for df in data[0]
        rg_values = Vector{Float64}(df.Rg)
        push!(passive_rg_values, compute_tail_mean(rg_values, tail_percent))
    end
    passive_rg = mean(passive_rg_values)
    println("Passive Rg (n_active=0, $(length(data[0])) runs): $(round(passive_rg, digits=4))")
    println()

    # Compute normalized mean Rg for each n_active
    n_active_list = Int[]
    mean_rg_norm = Float64[]
    std_rg_norm = Float64[]
    n_runs_list = Int[]

    println("Normalized mean Rg values (using last $(tail_percent)% of data):")
    for n_active in n_active_values
        dfs = data[n_active]
        n_runs = length(dfs)

        # Compute mean Rg for each run
        run_means = Float64[]
        for df in dfs
            rg_values = Vector{Float64}(df.Rg)
            push!(run_means, compute_tail_mean(rg_values, tail_percent))
        end

        # Aggregate across runs
        rg_mean = mean(run_means)
        rg_std = n_runs > 1 ? std(run_means) : 0.0

        # Normalize
        rg_mean_norm = rg_mean / passive_rg
        rg_std_norm = rg_std / passive_rg

        push!(n_active_list, n_active)
        push!(mean_rg_norm, rg_mean_norm)
        push!(std_rg_norm, rg_std_norm)
        push!(n_runs_list, n_runs)

        if n_runs > 1
            @printf("  n_active=%3d (%d runs): Rg/Rg_passive = %.4f ± %.4f\n",
                    n_active, n_runs, rg_mean_norm, rg_std_norm)
        else
            @printf("  n_active=%3d: Rg/Rg_passive = %.4f\n", n_active, rg_mean_norm)
        end
    end

    return n_active_list, mean_rg_norm, std_rg_norm, n_monomers
end

function process_per_run_normalization(args, tail_percent)
    println("Loading Rg data (grouped by run)...")
    base_states_dir = get(args, :base_states_dir, "")
    data_by_run, n_monomers = load_rg_data_by_run(args[:input_dir], args[:pattern];
                                                   base_states_dir=base_states_dir)
    run_ids = sort(collect(keys(data_by_run)))
    println("Loaded data for $(length(run_ids)) runs")
    println()

    # Filter to runs that have passive (n_active=0) data, skip others with warning
    valid_run_ids = Int[]
    skipped_runs = Int[]
    for run_id in run_ids
        if haskey(data_by_run[run_id], 0)
            push!(valid_run_ids, run_id)
        else
            push!(skipped_runs, run_id)
        end
    end

    if !isempty(skipped_runs)
        println("Warning: Skipping runs without passive (n_active=0) data: $(join(["run$r" for r in skipped_runs], ", "))")
        println()
    end

    if isempty(valid_run_ids)
        error("No runs have n_active=0 (passive) data. Cannot normalize.")
    end

    # Get all n_active values across valid runs only
    all_n_active = Set{Int}()
    for run_id in valid_run_ids
        union!(all_n_active, keys(data_by_run[run_id]))
    end
    n_active_values = sort(collect(all_n_active))

    # Compute per-run passive Rg values (only for valid runs)
    println("Per-run passive Rg values:")
    passive_rg_per_run = Dict{Int, Float64}()
    invalid_passive_runs = Int[]
    for run_id in valid_run_ids
        passive_df = data_by_run[run_id][0]
        rg_values = Vector{Float64}(passive_df.Rg)
        # Filter out NaN values
        rg_values = filter(!isnan, rg_values)
        if isempty(rg_values)
            @printf("  run%d: Rg_passive = NaN (no valid data, skipping)\n", run_id)
            push!(invalid_passive_runs, run_id)
            continue
        end
        passive_rg = compute_tail_mean(rg_values, tail_percent)
        if isnan(passive_rg)
            @printf("  run%d: Rg_passive = NaN (skipping)\n", run_id)
            push!(invalid_passive_runs, run_id)
            continue
        end
        passive_rg_per_run[run_id] = passive_rg
        @printf("  run%d: Rg_passive = %.4f\n", run_id, passive_rg)
    end

    # Remove runs with invalid passive Rg
    valid_run_ids = filter(r -> r ∉ invalid_passive_runs, valid_run_ids)

    if isempty(valid_run_ids)
        error("No runs have valid passive Rg data. Cannot normalize.")
    end

    if !isempty(invalid_passive_runs)
        println("  Skipped $(length(invalid_passive_runs)) runs with invalid passive data")
    end
    println()

    # Compute normalized Rg for each n_active
    # Each run is normalized by its own passive Rg
    n_active_list = Int[]
    mean_rg_norm = Float64[]
    std_rg_norm = Float64[]
    n_runs_list = Int[]

    println("Normalized mean Rg values (per-run normalization, last $(tail_percent)% of data):")
    for n_active in n_active_values
        # Collect normalized Rg from each valid run that has this n_active
        norm_values = Float64[]
        for run_id in valid_run_ids
            if haskey(data_by_run[run_id], n_active)
                df = data_by_run[run_id][n_active]
                rg_values = Vector{Float64}(df.Rg)
                # Filter out NaN values
                rg_values = filter(!isnan, rg_values)
                if isempty(rg_values)
                    continue
                end
                rg_mean = compute_tail_mean(rg_values, tail_percent)
                if isnan(rg_mean)
                    continue
                end
                rg_norm = rg_mean / passive_rg_per_run[run_id]
                push!(norm_values, rg_norm)
            end
        end

        if isempty(norm_values)
            continue
        end

        n_runs = length(norm_values)
        rg_mean_norm = mean(norm_values)
        rg_std_norm = n_runs > 1 ? std(norm_values) : 0.0

        push!(n_active_list, n_active)
        push!(mean_rg_norm, rg_mean_norm)
        push!(std_rg_norm, rg_std_norm)
        push!(n_runs_list, n_runs)

        if n_runs > 1
            @printf("  n_active=%3d (%d runs): Rg/Rg_passive = %.4f ± %.4f\n",
                    n_active, n_runs, rg_mean_norm, rg_std_norm)
        else
            @printf("  n_active=%3d: Rg/Rg_passive = %.4f\n", n_active, rg_mean_norm)
        end
    end

    return n_active_list, mean_rg_norm, std_rg_norm, n_monomers
end

function main()
    args = parse_commandline()

    mkpath(args[:output_dir])

    tail_percent = args[:tail_percent]
    per_run_norm = args[:per_run_normalization]

    if per_run_norm
        println("Using per-run normalization (each run normalized by its own passive Rg)")
        println()
        n_active_list, mean_rg_norm, std_rg_norm, n_monomers = process_per_run_normalization(args, tail_percent)
    else
        n_active_list, mean_rg_norm, std_rg_norm, n_monomers = process_global_normalization(args, tail_percent)
    end
    println()

    # Export CSV if requested
    if args[:export_csv]
        tail_str = @sprintf("tail%g", tail_percent)
        norm_str = per_run_norm ? "_perrun" : ""
        csv_base = "rg_normalized_vs_nactive_$(tail_str)$(norm_str)"
        csv_file = joinpath(args[:output_dir], "$(csv_base).csv")

        # Compute SEM (standard error of mean)
        # Note: std_rg_norm is already std across runs, SEM would need n_runs
        # For now, export std as the error measure
        df_out = DataFrame(
            n_active = n_active_list,
            rg_norm_mean = mean_rg_norm,
            rg_norm_std = std_rg_norm
        )
        CSV.write(csv_file, df_out)
        println("Data exported to: $csv_file")
        println()
    end

    # Create plot
    fig = Figure(size=(700, 500), fontsize=14)
    ax = Axis(fig[1, 1],
              xlabel="Number of active monomers (n_active)",
              ylabel="⟨Rg⟩ / Rg_passive",
              title="Normalized Radius of Gyration (N=$n_monomers, f=5.0)")

    scatter!(ax, n_active_list, mean_rg_norm, color=:steelblue, markersize=12)
    if any(std_rg_norm .> 0)
        errorbars!(ax, n_active_list, mean_rg_norm, std_rg_norm,
                   color=:steelblue, linewidth=1.5, whiskerwidth=8)
    end
    lines!(ax, n_active_list, mean_rg_norm, color=:steelblue, linewidth=2)

    # Reference line at 1.0
    hlines!(ax, [1.0], color=:gray, linestyle=:dash, linewidth=1, label="passive")

    xlims!(ax, -5, n_monomers + 5)

    # Save
    tail_str = @sprintf("tail%g", tail_percent)
    norm_str = per_run_norm ? "_perrun" : ""
    base = "rg_normalized_vs_nactive_$(tail_str)$(norm_str)"

    png_file = joinpath(args[:output_dir], "$(base).png")
    pdf_file = joinpath(args[:output_dir], "$(base).pdf")

    save(png_file, fig; px_per_unit=args[:dpi]/72)
    save(pdf_file, fig)

    println("Plots saved:")
    println("  $png_file")
    println("  $pdf_file")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
