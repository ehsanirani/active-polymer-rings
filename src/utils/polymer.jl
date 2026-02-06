using Statistics
using LinearAlgebra: norm

export create_bodies, create_initial_system

function create_bodies(params::Parameters)
    L = params.L
    boundary = CubicBoundary(L, L, L)

    n_particles = get_n_particles(params)
    activity = get_activity_vector(params)

    atoms = [Particle(id, 1.0, 1.0, 1.0, activity[id],
                      SVector(0.0, 0.0, 0.0),
                      SVector(0.0, 0.0, 0.0),
                      SVector(0, 0, 0)) for id in 1:n_particles]

    coords = initialize_polymer_coords(params, boundary)
    velocities = [random_velocity(1.0, params.KT) for _ in 1:n_particles]

    neighbor_finder = create_neighbor_finder(coords, params.rcut_nf)

    return SimBodies(atoms, coords, velocities, boundary, activity, neighbor_finder)
end

function initialize_polymer_coords(params::Parameters, boundary)
    if params.system_type == :single
        return initialize_single_ring(params, boundary)
    else
        return initialize_double_ring(params, boundary)
    end
end

function initialize_single_ring(params::Parameters, boundary)
    if params.init_method == :fourier
        return initialize_ring_fourier(params.n_monomers, boundary;
            k_max=params.init_kmax, bond_length=0.97)
    else  # :circle
        _, wrapped = place_polymer_ring(params.n_monomers, boundary; σ=0.97)
        return wrapped
    end
end

function initialize_double_ring(params::Parameters, boundary)
    if params.init_method == :fourier
        # Generate each ring separately with Fourier method
        n1 = params.n_monomers_1
        n2 = params.n_monomers_2

        ring1_coords = initialize_ring_fourier(n1, boundary;
            k_max=params.init_kmax, bond_length=0.97)
        ring2_coords = initialize_ring_fourier(n2, boundary;
            k_max=params.init_kmax, bond_length=0.97)

        # Offset second ring to avoid overlap (simple displacement)
        box_center = boundary.side_lengths ./ 2
        offset = SVector(3.0, 0.0, 0.0)
        ring2_coords = [c + offset for c in ring2_coords]

        return vcat(ring1_coords, ring2_coords)
    else  # :circle - use catenated placement
        return place_two_rings(params.n_monomers_1, params.n_monomers_2, boundary)
    end
end

function place_polymer_ring(N, boundary; σ=1.0, remove_rcm=true, pos_1=1, pos_2=2, pos0=zeros(3))
    box_center = boundary.side_lengths ./ 2
    theta = 2π/N
    radius_of_sphere = σ / sqrt(2*(1-cos(theta)))

    walk = [zeros(3) for _ in 1:N]
    wrapped_walk = [zeros(3) for _ in 1:N]

    for i in 1:N
        # Create mutable array from pos0 (handle both Array and SVector)
        particle_position = collect(pos0)
        # Use (i-1) to start at angle 0 and create a proper ring
        particle_position[pos_1] = cos(theta*(i-1))
        particle_position[pos_2] = sin(theta*(i-1))
        walk[i] = particle_position .* radius_of_sphere
    end

    if remove_rcm
        rcm = zeros(3)
        for j in 1:3
            rcm[j] = mean(walk[i][j] for i in 1:N)
        end
        for i in 1:N
            walk[i] = walk[i] .- rcm .+ box_center .+ pos0
            wrapped_walk[i] = wrap_coords(walk[i], boundary)
        end
    end

    swalk = [SVector(s0...) for s0 in walk]
    swrapped_walk = [SVector(s0...) for s0 in wrapped_walk]
    return swalk, swrapped_walk
end

function place_two_rings(n1::Int64, n2::Int64, boundary; σ=1.0)
    box_center = boundary.side_lengths ./ 2

    # Calculate radii for both rings
    theta1 = 2π / n1
    theta2 = 2π / n2
    radius1 = σ / sqrt(2 * (1 - cos(theta1)))
    radius2 = σ / sqrt(2 * (1 - cos(theta2)))

    # First ring in XY plane (horizontal)
    _, wrapped1 = place_polymer_ring(n1, boundary; σ=σ, pos_1=1, pos_2=2)

    # Second ring in XZ plane (vertical), centered at edge of first ring
    # This creates catenation (rings interlinked)
    center_offset = SVector(radius1, 0.0, 0.0)
    _, wrapped2 = place_polymer_ring(n2, boundary; σ=σ, pos_1=1, pos_2=3, pos0=center_offset)

    return vcat(wrapped1, wrapped2)
end

"""
    initialize_ring_fourier(N, boundary; k_max=10, bond_length=0.97)

Initialize a ring polymer using low-frequency Fourier modes.
This produces smooth, globally curved configurations that are
guaranteed unknotted for k_max ≤ 15.

The ring is constructed by summing Fourier modes with amplitudes
decaying as 1/k (Rouse-like), then rescaled to achieve target bond length.

Arguments:
- `N`: Number of monomers
- `boundary`: Simulation box boundary
- `k_max`: Maximum Fourier mode number (default: 10)
- `bond_length`: Target mean bond length (default: 0.97)

Returns vector of SVector{3,Float64} coordinates centered in the box.
"""
function initialize_ring_fourier(N::Int, boundary;
                                  k_max::Int=10,
                                  bond_length::Float64=0.97)
    # Initialize with zeros
    coords = [zeros(SVector{3,Float64}) for _ in 1:N]

    # Sum Fourier modes k=1 to k_max
    for k in 1:k_max
        σ_k = 1.0 / k  # Rouse-like amplitude decay
        a_k = SVector{3,Float64}(randn(3) .* σ_k)
        b_k = SVector{3,Float64}(randn(3) .* σ_k)

        for i in 1:N
            s = (i - 1) / N
            coords[i] = coords[i] + a_k * cos(2π * k * s) + b_k * sin(2π * k * s)
        end
    end

    # Compute center of mass and box center
    com = sum(coords) / N
    box_center = boundary.side_lengths ./ 2

    # Remove COM and center in box
    coords = [c - com + box_center for c in coords]

    # Compute current bond lengths and rescale to target
    bond_lengths = [norm(coords[mod1(i+1, N)] - coords[i]) for i in 1:N]
    mean_bond = mean(bond_lengths)
    scale = bond_length / mean_bond

    # Rescale around box center
    coords = [box_center + (c - box_center) * scale for c in coords]

    return coords
end

function create_neighbor_finder(coords::Vector, rcut::Float64)
    n_particles = length(coords)
    eligible = trues(n_particles, n_particles)

    neighbor_finder = CellListMapNeighborFinder(
        eligible=eligible,
        n_steps=1,
        dist_cutoff=rcut,
        x0=coords
    )
end

function create_initial_system(params::Parameters, sim_bodies::SimBodies; loggers)
    if params.system_type == :single
        return SingleRingBuilder(params, sim_bodies; loggers=loggers)
    else
        return DoubleRingBuilder(params, sim_bodies; loggers=loggers)
    end
end
