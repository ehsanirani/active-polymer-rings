using StaticArrays
using LinearAlgebra
using Statistics
using Random

export gyration_tensor_eigenvalues, compute_rg_timeseries
export compute_msd, compute_rs, compute_beta

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
