using Statistics

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
        return initialize_single_ring(params.n_monomers, boundary)
    else
        return initialize_double_ring(params.n_monomers_1, params.n_monomers_2, boundary)
    end
end

function initialize_single_ring(n_monomers::Int64, boundary)
    _, wrapped = place_polymer_ring(n_monomers, boundary; σ=0.97)
    return wrapped
end

function initialize_double_ring(n_monomers_1::Int64, n_monomers_2::Int64, boundary)
    return place_two_rings(n_monomers_1, n_monomers_2, boundary)
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