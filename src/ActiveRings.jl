module ActiveRings

using Molly
using Molly: AtomsCalculators
using StaticArrays
using ArgParse
using JLD2
using CellListMap
using LinearAlgebra
using CSV
using DataFrames
using Printf

export AbstractParameters, Parameters
export get_n_particles, get_n_active, get_activity_vector
export AbstractSystemBuilder, SingleRingBuilder, DoubleRingBuilder
export SimBodies, Particle
export create_bodies, create_neighbor_finder, create_initial_system
export simulate_system, create_simulation
export RgLogger, TangentLogger
export make_logging_schedule
export ActiveTangentForce, LangevinThermostat

# State I/O and checkpointing exports
export SimulationState, save_state, load_state, create_bodies_from_state
export CheckpointLogger, get_latest_checkpoint

# Analysis exports
export gyration_tensor_eigenvalues, compute_rg_timeseries
export compute_msd, compute_rs, compute_beta
export compute_autocorrelation, estimate_correlation_time
export load_simulation_data, split_rings_coords, split_rings_tangents
export save_rg_csv, save_msd_csv, save_rs_csv, save_beta_csv
export analyze_simulation, analyze_ring

# Core structures
include("parameters.jl")
include("simulation.jl")

# Forces and interactions
include(joinpath("forces", "active_force.jl"))
include(joinpath("forces", "langevin.jl"))

# Utility functions
include(joinpath("utils", "polymer.jl"))
include(joinpath("utils", "logging_schedule.jl"))
include(joinpath("utils", "loggers.jl"))
include(joinpath("utils", "state_io.jl"))
include(joinpath("utils", "checkpoint_logger.jl"))

# System builders
include(joinpath("systems", "single_ring.jl"))
include(joinpath("systems", "double_ring.jl"))

# Analysis
include(joinpath("analysis", "metrics.jl"))
include(joinpath("analysis", "io.jl"))
include(joinpath("analysis", "analyze.jl"))

end
