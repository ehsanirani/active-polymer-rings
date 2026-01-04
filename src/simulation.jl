using Molly
using StaticArrays

export SimBodies, Particle

mutable struct SimBodies{A,B,C,D,E,F}
    atoms::A
    coords::B
    velocities::C
    boundary::D
    activity::E
    neighbor_finder::F
end

SimBodies(; atoms, coords, velocities, boundary, activity, neighbor_finder) =
    SimBodies(atoms, coords, velocities, boundary, activity, neighbor_finder)

mutable struct Particle
    id::Int
    mass::Float64
    σ::Float64
    ϵ::Float64
    active::Bool
    tang_vec::SVector{3,Float64}
    rold::SVector{3,Float64}
    image::SVector{3,Int64}
end

function neighbor_finder(params, coords::Vector, neighbor_finder_rcut::Float64)
    n_particles = length(coords)
    nb_matrix = trues(n_particles, n_particles)

    for i in 1:n_particles-1
        nb_matrix[i, i+1] = false
        nb_matrix[i+1, i] = false
    end

    neighbor_finder = CellListMapNeighborFinder(
        nb_matrix=nb_matrix,
        n_steps=1,
        dist_cutoff=neighbor_finder_rcut,
        x0=coords
    )
end
