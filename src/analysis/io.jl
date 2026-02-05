using JLD2
using CSV
using DataFrames
using Printf
using StaticArrays

export load_simulation_data, load_msd_from_file, split_rings_coords, split_rings_tangents
export save_rg_csv, save_msd_csv, save_rs_csv, save_beta_csv

"""
    load_simulation_data(filepath::String; phase::Symbol=:active)

Load simulation results from JLD2 file.

Returns: (coords_history, tangents_history, params)
"""
function load_simulation_data(filepath::String; phase::Symbol=:active)
    jldopen(filepath, "r") do f
        phase_str = string(phase)

        coords = f["$(phase_str)/loggers"]["coords"].history
        tangents_raw = f["$(phase_str)/loggers"]["tangents"].tang_vecs
        params = f["params"]

        # Convert tangents to SVector format if needed
        tangents = convert_tangents_to_svector(tangents_raw)

        return coords, tangents, params
    end
end

"""
    load_msd_from_file(filepath::String; phase::Symbol=:active)

Load both time-averaged and non-time-averaged MSD from a JLD2 file.

Returns:
- Named tuple with msd_monomer, msd_com, msd_monomer_avg, msd_com_avg, lag_times
"""
function load_msd_from_file(filepath::String; phase::Symbol=:active)
    data = jldopen(filepath, "r")
    phase_str = string(phase)

    # Load non-time-averaged from logger
    msd_logger = data["$(phase_str)/loggers"]["msd"]
    msd_monomer = msd_logger.msd_monomer
    msd_com = msd_logger.msd_com

    # Load coords for time-averaged computation
    coords_history = data["$(phase_str)/loggers"]["coords"].history
    params = data["params"]
    close(data)

    # Compute time-averaged
    msd_monomer_avg = compute_msd(coords_history)
    msd_com_avg = compute_msd_com_timeaveraged(coords_history)

    # Time arrays
    dt = phase == :active ? params.dt : params.dt_thermal
    time_interval = dt * params.logger_steps
    lag_times = collect(1:length(msd_monomer)) .* time_interval

    return (msd_monomer=msd_monomer,
            msd_com=msd_com,
            msd_monomer_avg=msd_monomer_avg,
            msd_com_avg=msd_com_avg,
            lag_times=lag_times)
end

"""
Convert tangent vectors to SVector format for consistency
"""
function convert_tangents_to_svector(tangents_raw)
    if eltype(eltype(tangents_raw)) <: SVector
        return tangents_raw
    else
        # Convert Vector{Vector{Vector{Float64}}} to Vector{Vector{SVector{3,Float64}}}
        return [
            [SVector{3,Float64}(t) for t in frame]
            for frame in tangents_raw
        ]
    end
end

"""
    split_rings_coords(coords_history, n1::Int, n2::Int)

Split coordinates history into two rings for double ring systems.

Returns: (coords_ring1, coords_ring2)
"""
function split_rings_coords(coords_history::Vector{Vector{SVector{3, Float64}}}, n1::Int, n2::Int)
    n_frames = length(coords_history)

    coords_ring1 = [frame[1:n1] for frame in coords_history]
    coords_ring2 = [frame[(n1+1):(n1+n2)] for frame in coords_history]

    return coords_ring1, coords_ring2
end

"""
    split_rings_tangents(tangents_history, n1::Int, n2::Int)

Split tangent vectors history into two rings for double ring systems.

Returns: (tangents_ring1, tangents_ring2)
"""
function split_rings_tangents(tangents_history::Vector{Vector{SVector{3, Float64}}}, n1::Int, n2::Int)
    n_frames = length(tangents_history)

    tangents_ring1 = [frame[1:n1] for frame in tangents_history]
    tangents_ring2 = [frame[(n1+1):(n1+n2)] for frame in tangents_history]

    return tangents_ring1, tangents_ring2
end

"""
    save_rg_csv(run_name::String, ring_number::Int, dt::Float64, logger_steps::Int,
                Rg::Vector{Float64}, Rg1::Vector{Float64}, Rg2::Vector{Float64}, Rg3::Vector{Float64};
                output_dir::String="_data/csv")

Save radius of gyration components to CSV files.
"""
function save_rg_csv(run_name::String, ring_number::Int, dt::Float64, logger_steps::Int,
                     Rg::Vector{Float64}, Rg1::Vector{Float64}, Rg2::Vector{Float64}, Rg3::Vector{Float64};
                     output_dir::String="_data/csv")

    mkpath(output_dir)

    # File paths
    files = Dict(
        :Rg => joinpath(output_dir, "Rg_ring$(ring_number).csv"),
        :Rg1 => joinpath(output_dir, "Rg1_ring$(ring_number).csv"),
        :Rg2 => joinpath(output_dir, "Rg2_ring$(ring_number).csv"),
        :Rg3 => joinpath(output_dir, "Rg3_ring$(ring_number).csv"),
    )

    # Format function
    format_values(arr) = [@sprintf("%.4f", x) for x in arr]

    # Save each component
    for (key, data) in zip([:Rg, :Rg1, :Rg2, :Rg3], [Rg, Rg1, Rg2, Rg3])
        filename = files[key]

        if isfile(filename)
            df = CSV.read(filename, DataFrame)
            df[!, run_name] = format_values(data)
        else
            df = DataFrame(time = dt * collect(1:length(data)) * logger_steps)
            df[!, run_name] = format_values(data)
        end

        CSV.write(filename, df)
    end
end

"""
    save_msd_csv(run_name::String, ring_number::Int, dt::Float64, logger_steps::Int,
                 msd_array::Vector{Float64}; output_dir::String="_data/csv")

Save MSD data to CSV file.
"""
function save_msd_csv(run_name::String, ring_number::Int, dt::Float64, logger_steps::Int,
                      msd_array::Vector{Float64}; output_dir::String="_data/csv")

    mkpath(output_dir)

    filename = joinpath(output_dir, "MSD_ring$(ring_number).csv")
    format_values(arr) = [@sprintf("%.4f", x) for x in arr]

    if isfile(filename)
        df = CSV.read(filename, DataFrame)
        df[!, run_name] = format_values(msd_array)
    else
        n = length(msd_array)
        df = DataFrame(time = dt * collect(1:n) * logger_steps)
        df[!, run_name] = format_values(msd_array)
    end

    CSV.write(filename, df)
end

"""
    save_rs_csv(run_name::String, ring_number::Int, Rs_array::Vector{Float64};
                output_dir::String="_data/csv")

Save Rs (end-to-end distance) data to CSV file.
"""
function save_rs_csv(run_name::String, ring_number::Int, Rs_array::Vector{Float64};
                     output_dir::String="_data/csv")

    mkpath(output_dir)

    filename = joinpath(output_dir, "Rs_ring$(ring_number).csv")
    format_values(arr) = [@sprintf("%.4f", x) for x in arr]

    if isfile(filename)
        df = CSV.read(filename, DataFrame)
        df[!, run_name] = format_values(Rs_array)
    else
        n = length(Rs_array)
        df = DataFrame(s = collect(0:n-1))
        df[!, run_name] = format_values(Rs_array)
    end

    CSV.write(filename, df)
end

"""
    save_beta_csv(run_name::String, ring_number::Int, beta_array::Vector{Float64};
                  output_dir::String="_data/csv")

Save beta (tangent correlation) data to CSV file.
"""
function save_beta_csv(run_name::String, ring_number::Int, beta_array::Vector{Float64};
                       output_dir::String="_data/csv")

    mkpath(output_dir)

    filename = joinpath(output_dir, "beta_ring$(ring_number).csv")
    format_values(arr) = [@sprintf("%.4f", x) for x in arr]

    if isfile(filename)
        df = CSV.read(filename, DataFrame)
        df[!, run_name] = format_values(beta_array)
    else
        n = length(beta_array)
        df = DataFrame(s = collect(0:n-1))
        df[!, run_name] = format_values(beta_array)
    end

    CSV.write(filename, df)
end
