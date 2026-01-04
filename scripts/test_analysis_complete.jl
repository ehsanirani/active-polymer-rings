#!/usr/bin/env julia

using ActiveRings
using JLD2
using StaticArrays
using Random
using Statistics
using CSV
using DataFrames
using Printf

# Load legacy functions
include("../2Rings-Codes/Three_base_function.jl")
include("../2Rings-Codes/Read_Simres.jl")

println("="^60)
println("COMPREHENSIVE ANALYSIS MODULE TEST")
println("="^60)

# Test file
filepath = "_data/jld2/double_10_8_3_2_0.0_1.0_test_double.jld2"

println("\nLoading simulation data from: $filepath")

# Load parameters for metadata
params = jldopen(filepath, "r") do f
    f["params"]
end

println("  System: $(params.system_type)")
println("  Ring 1: $(params.n_monomers_1) monomers, $(params.n_active_1) active")
println("  Ring 2: $(params.n_monomers_2) monomers, $(params.n_active_2) active")
println("  dt = $(params.dt), logger_steps = $(params.logger_steps)")

# Create output directories
mkpath("_data/csv_legacy")
mkpath("_data/csv_new")

println("\n" * "="^60)
println("RUNNING LEGACY ANALYSIS")
println("="^60)

# Load data with legacy method
coords_active = jldopen(filepath, "r") do f
    f["active/loggers"]["coords"].history
end

tangents_active = jldopen(filepath, "r") do f
    tang_raw = f["active/loggers"]["tangents"].tang_vecs
    [[Vector{Float64}(t) for t in frame] for frame in tang_raw]
end

# Split for double ring
n_monomers = length(coords_active[1])
n1 = params.n_monomers_1
n2 = params.n_monomers_2

coords_ring1_legacy = [frame[1:n1] for frame in coords_active]
coords_ring2_legacy = [frame[(n1+1):(n1+n2)] for frame in coords_active]

tangents_ring1_legacy = [frame[1:n1] for frame in tangents_active]
tangents_ring2_legacy = [frame[(n1+1):(n1+n2)] for frame in tangents_active]

# Helper function to save legacy CSV (adapted from Save_CSV_Functions.jl)
function save_legacy_rgs(run_name::String, ring_number::Int, dt::Float64,
    logger_steps::Int, Rg::Vector{Float64}, Rg1::Vector{Float64},
    Rg2::Vector{Float64}, Rg3::Vector{Float64})

    base_path = "_data/csv_legacy/"
    mkpath(base_path)

    files = Dict(
        :Rg => base_path * "Rg_ring$(ring_number).csv",
        :Rg1 => base_path * "Rg1_ring$(ring_number).csv",
        :Rg2 => base_path * "Rg2_ring$(ring_number).csv",
        :Rg3 => base_path * "Rg3_ring$(ring_number).csv",
    )

    format4(arr) = [@sprintf("%.4f", x) for x in arr]

    for (key, data) in zip([:Rg, :Rg1, :Rg2, :Rg3], [Rg, Rg1, Rg2, Rg3])
        df = DataFrame(time=dt * collect(1:length(data)) * logger_steps)
        df[!, run_name] = format4(data)
        CSV.write(files[key], df)
    end
end

function save_legacy_msd(run_name::String, ring_number::Int, dt::Float64,
    logger_steps::Int, msd_array::Vector{Float64})
    base_path = "_data/csv_legacy/"
    filename = base_path * "MSD_ring$(ring_number).csv"
    format4(arr) = [@sprintf("%.4f", x) for x in arr]

    n = length(msd_array)
    df = DataFrame(time=dt * collect(1:n) * logger_steps)
    df[!, run_name] = format4(msd_array)
    CSV.write(filename, df)
end

function save_legacy_rs(run_name::String, ring_number::Int, Rs_array::Vector{Float64})
    base_path = "_data/csv_legacy/"
    filename = base_path * "Rs_ring$(ring_number).csv"
    format4(arr) = [@sprintf("%.4f", x) for x in arr]

    n = length(Rs_array)
    df = DataFrame(s=collect(0:n-1))
    df[!, run_name] = format4(Rs_array)
    CSV.write(filename, df)
end

function save_legacy_beta(run_name::String, ring_number::Int, beta_array::Vector{Float64})
    base_path = "_data/csv_legacy/"
    filename = base_path * "beta_ring$(ring_number).csv"
    format4(arr) = [@sprintf("%.4f", x) for x in arr]

    n = length(beta_array)
    df = DataFrame(s=collect(0:n-1))
    df[!, run_name] = format4(beta_array)
    CSV.write(filename, df)
end

run_name = "test_legacy"

# Process Ring 1
println("\nProcessing Ring 1 (legacy)...")
rg1_leg, rg2_leg, rg3_leg, rg_leg = create_Rgs_array(coords_ring1_legacy)
save_legacy_rgs(run_name, 1, params.dt, params.logger_steps, rg_leg, rg1_leg, rg2_leg, rg3_leg)

msd_leg = create_msd_array(coords_ring1_legacy)
save_legacy_msd(run_name, 1, params.dt, params.logger_steps, msd_leg)

Rs_leg = Rs(coords_ring1_legacy)
save_legacy_rs(run_name, 1, Rs_leg)

beta_leg = beta(tangents_ring1_legacy)
save_legacy_beta(run_name, 1, beta_leg)

println("  ✓ Ring 1 complete")

# Process Ring 2
println("Processing Ring 2 (legacy)...")
rg1_leg2, rg2_leg2, rg3_leg2, rg_leg2 = create_Rgs_array(coords_ring2_legacy)
save_legacy_rgs(run_name, 2, params.dt, params.logger_steps, rg_leg2, rg1_leg2, rg2_leg2, rg3_leg2)

msd_leg2 = create_msd_array(coords_ring2_legacy)
save_legacy_msd(run_name, 2, params.dt, params.logger_steps, msd_leg2)

Rs_leg2 = Rs(coords_ring2_legacy)
save_legacy_rs(run_name, 2, Rs_leg2)

beta_leg2 = beta(tangents_ring2_legacy)
save_legacy_beta(run_name, 2, beta_leg2)

println("  ✓ Ring 2 complete")

println("\n" * "="^60)
println("RUNNING NEW UNIFIED ANALYSIS")
println("="^60)

# Run new analysis
analyze_simulation(filepath; phase=:active, run_name="test_new", output_dir="_data/csv_new")

println("\n" * "="^60)
println("COMPARING CSV OUTPUTS")
println("="^60)

# Compare all CSV files
metrics = ["Rg", "Rg1", "Rg2", "Rg3", "MSD", "Rs", "beta"]
rings = [1, 2]

all_match = true

for ring in rings
    println("\n--- Ring $ring ---")
    for metric in metrics
        legacy_file = "_data/csv_legacy/$(metric)_ring$(ring).csv"
        new_file = "_data/csv_new/$(metric)_ring$(ring).csv"

        if isfile(legacy_file) && isfile(new_file)
            df_legacy = CSV.read(legacy_file, DataFrame)
            df_new = CSV.read(new_file, DataFrame)

            # Extract data columns (skip 'time' or 's' column)
            # Convert to Float64 in case they're strings
            legacy_col = df_legacy[!, 2]
            new_col = df_new[!, 2]

            legacy_values = legacy_col isa Vector{String} ? parse.(Float64, legacy_col) : Float64.(legacy_col)
            new_values = new_col isa Vector{String} ? parse.(Float64, new_col) : Float64.(new_col)

            # Compare
            match = isapprox(legacy_values, new_values, rtol=1e-10)
            max_diff = maximum(abs.(legacy_values .- new_values))

            status = match ? "✓ MATCH" : "✗ MISMATCH"
            println("  $metric: $status (max diff: $(max_diff))")

            global all_match = all_match && match
        else
            println("  $metric: ⚠ Missing file")
            global all_match = false
        end
    end
end

println("\n" * "="^60)
if all_match
    println("✅ ALL TESTS PASSED - New analysis matches legacy exactly!")
else
    println("❌ SOME TESTS FAILED - Check output above")
end
println("="^60)
