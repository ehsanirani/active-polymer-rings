#!/usr/bin/env julia

using ActiveRings
using JLD2
using StaticArrays
using Random
using Statistics

# Load legacy functions
include("../2Rings-Codes/Three_base_function.jl")
include("../2Rings-Codes/Read_Simres.jl")

println("Loading test data...")
filepath = "_data/jld2/single_20_5_0.0_1.0_test_analysis.jld2"

# Load with new unified system
coords_new, tangents_new, params = load_simulation_data(filepath; phase=:active)

# Load with legacy system
coords_legacy = jldopen(filepath, "r") do f
    f["active/loggers"]["coords"].history
end

tangents_legacy = jldopen(filepath, "r") do f
    # Convert SVector to Vector{Float64} for legacy code
    tang_raw = f["active/loggers"]["tangents"].tang_vecs
    [[Vector{Float64}(t) for t in frame] for frame in tang_raw]
end

println("\n=== Testing Rg Calculation ===")
# New implementation
Rg1_new, Rg2_new, Rg3_new, Rg_new = compute_rg_timeseries(coords_new)

# Legacy implementation
rg1_leg, rg2_leg, rg3_leg, rg_leg = create_Rgs_array(coords_legacy)

println("Rg match: ", isapprox(Rg_new, rg_leg, rtol=1e-10))
println("Rg1 match: ", isapprox(Rg1_new, rg1_leg, rtol=1e-10))
println("Max difference: ", maximum(abs.(Rg_new .- rg_leg)))

println("\n=== Testing MSD Calculation ===")
# New implementation
msd_new = compute_msd(coords_new)

# Legacy implementation
msd_leg = create_msd_array(coords_legacy)

println("MSD match: ", isapprox(msd_new, msd_leg, rtol=1e-10))
println("Max difference: ", maximum(abs.(msd_new .- msd_leg)))

println("\n=== Testing Rs Calculation ===")
# New implementation
Rs_new = compute_rs(coords_new)

# Legacy implementation
Rs_leg = Rs(coords_legacy)

println("Rs match: ", isapprox(Rs_new, Rs_leg, rtol=1e-10))
println("Max difference: ", maximum(abs.(Rs_new .- Rs_leg)))

println("\n=== Testing Beta Calculation ===")
# New implementation
beta_new = compute_beta(tangents_new)

# Legacy implementation
beta_leg = beta(tangents_legacy)

println("Beta match: ", isapprox(beta_new, beta_leg, rtol=1e-10))
println("Max difference: ", maximum(abs.(beta_new .- beta_leg)))

println("\n✅ All tests completed!")
