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
println("SINGLE RING ANALYSIS TEST")
println("="^60)

# Test file
filepath = "_data/jld2/single_20_5_0.0_1.0_test_analysis.jld2"

println("\nLoading simulation data from: $filepath")

params = jldopen(filepath, "r") do f
    f["params"]
end

println("  System: $(params.system_type)")
println("  Monomers: $(params.n_monomers), Active: $(params.n_active)")
traj_int = hasproperty(params, :traj_interval) ? params.traj_interval : params.logger_steps
println("  dt = $(params.dt), traj_interval = $(traj_int)")

# Create output directories
mkpath("_data/csv_legacy_single")
mkpath("_data/csv_new_single")

println("\n" * "="^60)
println("RUNNING LEGACY ANALYSIS")
println("="^60)

# Load data
coords_active = jldopen(filepath, "r") do f
    f["active/loggers"]["coords"].history
end

tangents_active = jldopen(filepath, "r") do f
    tang_raw = f["active/loggers"]["tangents"].tang_vecs
    [[Vector{Float64}(t) for t in frame] for frame in tang_raw]
end

# Legacy save functions
function save_legacy_rgs(run_name::String, ring_number::Int, dt::Float64,
    logger_steps::Int, Rg::Vector{Float64}, Rg1::Vector{Float64},
    Rg2::Vector{Float64}, Rg3::Vector{Float64}; base_path="_data/csv_legacy_single/")
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
    logger_steps::Int, msd_array::Vector{Float64}; base_path="_data/csv_legacy_single/")
    filename = base_path * "MSD_ring$(ring_number).csv"
    format4(arr) = [@sprintf("%.4f", x) for x in arr]
    n = length(msd_array)
    df = DataFrame(time=dt * collect(1:n) * logger_steps)
    df[!, run_name] = format4(msd_array)
    CSV.write(filename, df)
end

function save_legacy_rs(run_name::String, ring_number::Int, Rs_array::Vector{Float64}; base_path="_data/csv_legacy_single/")
    filename = base_path * "Rs_ring$(ring_number).csv"
    format4(arr) = [@sprintf("%.4f", x) for x in arr]
    n = length(Rs_array)
    df = DataFrame(s=collect(0:n-1))
    df[!, run_name] = format4(Rs_array)
    CSV.write(filename, df)
end

function save_legacy_beta(run_name::String, ring_number::Int, beta_array::Vector{Float64}; base_path="_data/csv_legacy_single/")
    filename = base_path * "beta_ring$(ring_number).csv"
    format4(arr) = [@sprintf("%.4f", x) for x in arr]
    n = length(beta_array)
    df = DataFrame(s=collect(0:n-1))
    df[!, run_name] = format4(beta_array)
    CSV.write(filename, df)
end

run_name = "test_legacy"

println("Processing single ring (legacy)...")
rg1_leg, rg2_leg, rg3_leg, rg_leg = create_Rgs_array(coords_active)
save_legacy_rgs(run_name, 1, params.dt, traj_int, rg_leg, rg1_leg, rg2_leg, rg3_leg)

msd_leg = create_msd_array(coords_active)
save_legacy_msd(run_name, 1, params.dt, traj_int, msd_leg)

Rs_leg = Rs(coords_active)
save_legacy_rs(run_name, 1, Rs_leg)

beta_leg = beta(tangents_active)
save_legacy_beta(run_name, 1, beta_leg)

println("  ✓ Legacy analysis complete")

println("\n" * "="^60)
println("RUNNING NEW UNIFIED ANALYSIS")
println("="^60)

analyze_simulation(filepath; phase=:active, run_name="test_new", output_dir="_data/csv_new_single")

println("\n" * "="^60)
println("COMPARING CSV OUTPUTS")
println("="^60)

metrics = ["Rg", "Rg1", "Rg2", "Rg3", "MSD", "Rs", "beta"]
all_match = true

println("\n--- Single Ring ---")
for metric in metrics
    legacy_file = "_data/csv_legacy_single/$(metric)_ring1.csv"
    new_file = "_data/csv_new_single/$(metric)_ring1.csv"

    if isfile(legacy_file) && isfile(new_file)
        df_legacy = CSV.read(legacy_file, DataFrame)
        df_new = CSV.read(new_file, DataFrame)

        legacy_col = df_legacy[!, 2]
        new_col = df_new[!, 2]

        legacy_values = legacy_col isa Vector{String} ? parse.(Float64, legacy_col) : Float64.(legacy_col)
        new_values = new_col isa Vector{String} ? parse.(Float64, new_col) : Float64.(new_col)

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

println("\n" * "="^60)
if all_match
    println("✅ ALL TESTS PASSED - Single ring analysis matches legacy!")
else
    println("❌ SOME TESTS FAILED - Check output above")
end
println("="^60)
