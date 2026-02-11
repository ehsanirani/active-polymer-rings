using StaticArrays
using LinearAlgebra
using Statistics
using Random

export gyration_tensor_eigenvalues, compute_rg_timeseries
export compute_msd, compute_msd_com_timeaveraged, compute_msd_com_frame, compute_rs, compute_beta
export compute_autocorrelation, estimate_correlation_time

"""
    gyration_tensor_eigenvalues(positions::Vector{SVector{3, Float64}})

Calculate the square roots of the sorted eigenvalues of the gyration tensor
and the total radius of gyration for a set of 3D positions.

Returns: (Rg1, Rg2, Rg3, Rg_total) where Rg1 >= Rg2 >= Rg3
"""
function gyration_tensor_eigenvalues(positions::Vector{SVector{3, Float64}})
    N = length(positions)

    # Compute center of mass
    r_cm = sum(positions) / N

    # Build gyration tensor
    S = zeros(3, 3)
    for r in positions
        dr = r - r_cm
        S += dr * dr'
    end
    S /= N

    # Compute eigenvalues and sort in descending order
    sorted_vals = sort(eigen(S).values, rev=true)

    # Total squared radius of gyration
    rg2 = sum(sorted_vals)

    return sqrt(sorted_vals[1]), sqrt(sorted_vals[2]), sqrt(sorted_vals[3]), sqrt(rg2)
end

"""
    compute_rg_timeseries(coords_history::Vector{Vector{SVector{3, Float64}}})

Compute radius of gyration components for each frame in the trajectory.

Returns: (Rg1, Rg2, Rg3, Rg_total) arrays
"""
function compute_rg_timeseries(coords_history::Vector{Vector{SVector{3, Float64}}})
    n_frames = length(coords_history)

    Rg1 = Vector{Float64}(undef, n_frames)
    Rg2 = Vector{Float64}(undef, n_frames)
    Rg3 = Vector{Float64}(undef, n_frames)
    Rg_total = Vector{Float64}(undef, n_frames)

    for t in 1:n_frames
        r1, r2, r3, rtot = gyration_tensor_eigenvalues(coords_history[t])
        Rg1[t] = r1
        Rg2[t] = r2
        Rg3[t] = r3
        Rg_total[t] = rtot
    end

    return Rg1, Rg2, Rg3, Rg_total
end

"""
    compute_msd(coords_history::Vector{Vector{SVector{3, Float64}}}, max_samples::Int=1000)

Calculate mean squared displacement (MSD) as a function of time lag.

Uses random sampling of frame pairs to improve performance for long trajectories.
"""
function compute_msd(coords_history::Vector{Vector{SVector{3, Float64}}}, max_samples::Int=1000)
    n_frames = length(coords_history)
    N = length(coords_history[1])

    msd_array = Vector{Float64}(undef, n_frames-1)

    for t in 1:n_frames-1
        total_msd = 0.0

        # Sample frame pairs to avoid O(n_frames^2) complexity
        m = min(max_samples, n_frames - t)
        idxs = randperm(n_frames - t)[1:m]

        for j in idxs
            diff = coords_history[j+t] .- coords_history[j]
            total_msd += sum(d -> dot(d, d), diff)
        end

        msd_array[t] = total_msd / (N * m)
    end

    return msd_array
end

"""
    compute_msd_com_timeaveraged(coords_history::Vector{Vector{SVector{3, Float64}}}, max_samples::Int=1000)

Compute time-averaged MSD for the center of mass.

Uses multiple time origins to improve statistics.

Arguments:
- coords_history: Vector of coordinate frames
- max_samples: Maximum number of time origins to sample

Returns:
- msd_array: COM MSD values for each lag time
"""
function compute_msd_com_timeaveraged(coords_history::Vector{Vector{SVector{3, Float64}}}, max_samples::Int=1000)
    n_frames = length(coords_history)

    # Compute center of mass for each frame
    com_history = [sum(coords) / length(coords) for coords in coords_history]

    msd_array = Vector{Float64}(undef, n_frames-1)

    for t in 1:n_frames-1
        total_msd = 0.0

        # Sample time origins
        m = min(max_samples, n_frames - t)
        idxs = randperm(n_frames - t)[1:m]

        for j in idxs
            diff = com_history[j+t] - com_history[j]
            total_msd += dot(diff, diff)
        end

        msd_array[t] = total_msd / m
    end

    return msd_array
end

"""
    compute_msd_com_frame(coords_history::Vector{Vector{SVector{3, Float64}}}, max_samples::Int=1000)

Compute time-averaged MSD of monomers in the Center of Mass reference frame.

This removes the overall polymer drift, showing only internal fluctuations.

Arguments:
- coords_history: Vector of coordinate frames
- max_samples: Maximum number of time origins to sample

Returns:
- msd_array: COM-frame MSD values for each lag time
"""
function compute_msd_com_frame(coords_history::Vector{Vector{SVector{3, Float64}}}, max_samples::Int=1000)
    n_frames = length(coords_history)
    N = length(coords_history[1])

    # Step 1: Transform all frames to COM reference frame
    coords_com_frame = Vector{Vector{SVector{3,Float64}}}(undef, n_frames)
    for t in 1:n_frames
        com = sum(coords_history[t]) / N
        coords_com_frame[t] = coords_history[t] .- Ref(com)
    end

    # Step 2: Compute time-averaged MSD on transformed coordinates
    msd_array = Vector{Float64}(undef, n_frames-1)
    for τ in 1:n_frames-1
        total_msd = 0.0
        m = min(max_samples, n_frames - τ)
        idxs = randperm(n_frames - τ)[1:m]

        for j in idxs
            diff = coords_com_frame[j+τ] .- coords_com_frame[j]
            total_msd += sum(d -> dot(d, d), diff)
        end
        msd_array[τ] = total_msd / (N * m)
    end
    return msd_array
end

"""
    compute_rs(coords_history::Vector{Vector{SVector{3, Float64}}})

Compute end-to-end distance as a function of contour distance s.

Uses last 20% of frames for equilibrated statistics.
"""
function compute_rs(coords_history::Vector{Vector{SVector{3, Float64}}})
    n_frames = length(coords_history)
    n_monomers = length(coords_history[1])
    max_s = div(n_monomers, 2)
    s_values = collect(0:max_s)

    # Use last 20% of frames for equilibrated statistics
    last_n = div(n_frames, 5)
    frame_range = (n_frames - last_n + 1):n_frames

    Rs = zeros(Float64, length(s_values))

    for (s_idx, s) in enumerate(s_values)
        values = Float64[]

        for t in frame_range
            coords = coords_history[t]

            for i in 1:n_monomers
                j = mod(i + s - 1, n_monomers) + 1
                d = coords[j] - coords[i]
                push!(values, dot(d, d))
            end
        end

        Rs[s_idx] = sqrt(mean(values))
    end

    return Rs
end

"""
    compute_beta(tangents_history::Vector{Vector{SVector{3, Float64}}})

Compute tangent-tangent correlation function as a function of contour distance s.

Uses last 20% of frames for equilibrated statistics.
"""
function compute_beta(tangents_history::Vector{Vector{SVector{3, Float64}}})
    n_frames = length(tangents_history)
    n_monomers = length(tangents_history[1])
    max_s = div(n_monomers, 2)
    s_values = collect(0:max_s)

    # Use last 20% of frames
    last_n = div(n_frames, 5)
    frame_range = (n_frames - last_n + 1):n_frames

    correlations = zeros(Float64, length(s_values))

    for (s_idx, s) in enumerate(s_values)
        values = Float64[]

        for t in frame_range
            tangents = tangents_history[t]

            for i in 1:n_monomers
                j = mod(i + s - 1, n_monomers) + 1
                push!(values, dot(tangents[i], tangents[j]))
            end
        end

        correlations[s_idx] = mean(values)
    end

    return correlations
end

"""
    compute_autocorrelation(signal::Vector{Float64}, max_lag::Int=0) -> Vector{Float64}

Compute normalized autocorrelation function of a time series.

Returns ACF[τ] for τ = 0, 1, ..., max_lag where ACF[0] = 1.0 by definition.

# Arguments
- `signal`: Time series data
- `max_lag`: Maximum lag to compute (default: n ÷ 2)

# Returns
- Vector of ACF values from lag 0 to max_lag
"""
function compute_autocorrelation(signal::Vector{Float64}, max_lag::Int=0)
    n = length(signal)
    max_lag = max_lag > 0 ? min(max_lag, n-1) : n ÷ 2

    # Subtract mean
    μ = mean(signal)
    x = signal .- μ

    # Variance (for normalization)
    var_x = dot(x, x) / n

    acf = Vector{Float64}(undef, max_lag + 1)
    for τ in 0:max_lag
        # ACF(τ) = (1/N) Σ_t x(t) * x(t+τ)
        acf[τ+1] = dot(x[1:n-τ], x[τ+1:n]) / ((n - τ) * var_x)
    end

    return acf
end

"""
    estimate_correlation_time(acf::Vector{Float64}; method::Symbol=:integrated) -> Float64

Estimate correlation time from an autocorrelation function.

# Methods
- `:first_zero` - lag where ACF first crosses zero
- `:e_folding` - lag where ACF ≈ 1/e ≈ 0.368
- `:integrated` - ∫ ACF(τ) dτ from 0 to first zero (most robust, default)

# Arguments
- `acf`: Autocorrelation function values (from compute_autocorrelation)
- `method`: Method for estimating correlation time

# Returns
- Estimated correlation time in units of lag indices
"""
function estimate_correlation_time(acf::Vector{Float64}; method::Symbol=:integrated)
    if method == :first_zero
        idx = findfirst(x -> x <= 0, acf)
        return idx === nothing ? Float64(length(acf) - 1) : Float64(idx - 1)
    elseif method == :e_folding
        idx = findfirst(x -> x <= 1/ℯ, acf)
        return idx === nothing ? Float64(length(acf) - 1) : Float64(idx - 1)
    elseif method == :integrated
        # Integrated autocorrelation time (more robust)
        idx = findfirst(x -> x <= 0, acf)
        cutoff = idx === nothing ? length(acf) : idx
        return sum(acf[1:cutoff])
    else
        error("Unknown method: $method. Use :first_zero, :e_folding, or :integrated")
    end
end
