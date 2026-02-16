using JLD2
using CSV
using DataFrames
using Printf
using StaticArrays

export load_simulation_data, load_msd_from_file, split_rings_coords, split_rings_tangents
export save_rg_csv, save_msd_csv, save_rs_csv, save_beta_csv, format_with_precision

"""
    format_with_precision(x::Real, precision::Int)

Format a number with the specified number of decimal places.
"""
function format_with_precision(x::Real, precision::Int)
    fmt = Printf.Format("%.$(precision)f")
    return Printf.format(fmt, x)
end

"""
    load_simulation_data(filepath::String; phase::Symbol=:active)

Load simulation results from JLD2 file.

Returns: (coords_history, tangents_history, params)
Note: coords_history may be `nothing` if coords were not stored (msd_time_averaged disabled).
"""
function load_simulation_data(filepath::String; phase::Symbol=:active)
    jldopen(filepath, "r") do f
        phase_str = string(phase)
        loggers = f["$(phase_str)/loggers"]

        # Coords may not be present if msd_time_averaged was disabled
        coords = haskey(loggers, "coords") ? loggers["coords"].history : nothing
        tangents_raw = loggers["tangents"].tang_vecs
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
    loggers = data["$(phase_str)/loggers"]
    msd_logger = loggers["msd"]
    msd_monomer = msd_logger.msd_monomer
    msd_com = msd_logger.msd_com

    # Coords may not be present if msd_time_averaged was disabled
    has_coords = haskey(loggers, "coords")
    coords_history = has_coords ? loggers["coords"].history : nothing
    params = data["params"]
    close(data)

    # Check flags (with backward-compatible defaults)
    do_com = hasproperty(params, :msd_com) ? params.msd_com : true
    # Time-averaged MSD requires coords - disable if not available
    do_timeavg = hasproperty(params, :msd_time_averaged) ? params.msd_time_averaged : true
    if do_timeavg && !has_coords
        do_timeavg = false
    end

    # Compute time-averaged (if enabled and coords available)
    msd_monomer_avg = do_timeavg ? compute_msd(coords_history) : Float64[]
    msd_com_avg = (do_timeavg && do_com) ? compute_msd_com_timeaveraged(coords_history) : Float64[]

    # Time arrays — use step_indices if available, fall back to traj_interval
    dt = phase == :active ? params.dt : params.dt_thermal
    if hasproperty(msd_logger, :step_indices) && !isempty(msd_logger.step_indices)
        lag_times = msd_logger.step_indices .* dt
    else
        # Fallback for old data files
        traj_int = hasproperty(params, :traj_interval) ? params.traj_interval : params.logger_steps
        time_interval = dt * traj_int
        lag_times = collect(1:length(msd_monomer)) .* time_interval
    end

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
                output_dir::String="_data/csv", precision::Int=4)

Save radius of gyration components to CSV files.
"""
function save_rg_csv(run_name::String, ring_number::Int, dt::Float64, logger_steps::Int,
                     Rg::Vector{Float64}, Rg1::Vector{Float64}, Rg2::Vector{Float64}, Rg3::Vector{Float64};
                     output_dir::String="_data/csv", precision::Int=4)

    mkpath(output_dir)

    # File paths
    files = Dict(
        :Rg => joinpath(output_dir, "Rg_ring$(ring_number).csv"),
        :Rg1 => joinpath(output_dir, "Rg1_ring$(ring_number).csv"),
        :Rg2 => joinpath(output_dir, "Rg2_ring$(ring_number).csv"),
        :Rg3 => joinpath(output_dir, "Rg3_ring$(ring_number).csv"),
    )

    # Format function with configurable precision
    format_values(arr) = [format_with_precision(x, precision) for x in arr]

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
                 msd_array::Vector{Float64}; output_dir::String="_data/csv", precision::Int=4)

Save MSD data to CSV file.
"""
function save_msd_csv(run_name::String, ring_number::Int, dt::Float64, logger_steps::Int,
                      msd_array::Vector{Float64}; output_dir::String="_data/csv", precision::Int=4)

    mkpath(output_dir)

    filename = joinpath(output_dir, "MSD_ring$(ring_number).csv")
    format_values(arr) = [format_with_precision(x, precision) for x in arr]

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
                output_dir::String="_data/csv", precision::Int=4)

Save Rs (end-to-end distance) data to CSV file.
"""
function save_rs_csv(run_name::String, ring_number::Int, Rs_array::Vector{Float64};
                     output_dir::String="_data/csv", precision::Int=4)

    mkpath(output_dir)

    filename = joinpath(output_dir, "Rs_ring$(ring_number).csv")
    format_values(arr) = [format_with_precision(x, precision) for x in arr]

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
                  output_dir::String="_data/csv", precision::Int=4)

Save beta (tangent correlation) data to CSV file.
"""
function save_beta_csv(run_name::String, ring_number::Int, beta_array::Vector{Float64};
                       output_dir::String="_data/csv", precision::Int=4)

    mkpath(output_dir)

    filename = joinpath(output_dir, "beta_ring$(ring_number).csv")
    format_values(arr) = [format_with_precision(x, precision) for x in arr]

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
