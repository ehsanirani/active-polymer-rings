using JLD2
using StaticArrays
using Dates

export SimulationState, save_state, load_state, create_bodies_from_state

"""
Captures the complete state of a simulation for saving/loading.

Contains all information needed to resume a simulation from a checkpoint:
- Coordinates and velocities
- MSD tracking state (reference positions and images)
- Activity pattern
- System configuration
- Metadata about when/where saved
"""
struct SimulationState
    # Coordinates and dynamics
    coords::Vector{SVector{3,Float64}}
    velocities::Vector{SVector{3,Float64}}
    L::Float64

    # MSD tracking (for continuing MSD computation)
    rold::Vector{SVector{3,Float64}}
    image::Vector{SVector{3,Int64}}
    tang_vec::Vector{SVector{3,Float64}}

    # Activity pattern
    activity::Vector{Bool}
    n_active::Int64
    activity_pattern::Symbol  # :random or :block

    # System info
    system_type::Symbol
    n_monomers::Int64
    n_monomers_1::Int64
    n_monomers_2::Int64

    # Metadata
    original_params::Parameters
    save_timestamp::String
    save_phase::Symbol          # :thermalization or :active
    current_step::Int64         # Step number when saved
    total_steps::Int64          # Total steps planned
end

"""
    save_state(filepath, sim_bodies, sys, params; phase, step, total_steps)

Save simulation state to a JLD2 file.

Arguments:
- `filepath`: Path to save the state file
- `sim_bodies`: SimBodies struct with current simulation state
- `sys`: Molly System (to extract tangent vectors from atoms)
- `params`: Parameters struct
- `phase`: Current phase (:thermalization or :active)
- `step`: Current step number
- `total_steps`: Total steps planned for this phase
"""
function save_state(filepath::String, sim_bodies::SimBodies, sys, params::Parameters;
                    phase::Symbol, step::Int64, total_steps::Int64)
    # Extract state from system
    n_particles = get_n_particles(params)

    # Extract tangent vectors and MSD tracking data from atoms
    tang_vec = [sys.atoms[i].tang_vec for i in 1:n_particles]
    rold = [sys.atoms[i].rold for i in 1:n_particles]
    image = [sys.atoms[i].image for i in 1:n_particles]

    state = SimulationState(
        # Coordinates and dynamics
        Vector(sys.coords),
        Vector(sys.velocities),
        sim_bodies.boundary.side_lengths[1],
        # MSD tracking
        rold,
        image,
        tang_vec,
        # Activity pattern
        Vector(sim_bodies.activity),
        get_n_active(params),
        params.activity_pattern,
        # System info
        params.system_type,
        params.n_monomers,
        params.n_monomers_1,
        params.n_monomers_2,
        # Metadata
        params,
        Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"),
        phase,
        step,
        total_steps
    )

    # Ensure directory exists
    mkpath(dirname(filepath))

    # Save to JLD2
    jldopen(filepath, "w"; compress=true) do file
        file["state"] = state
    end

    return state
end

"""
    load_state(filepath) -> SimulationState

Load simulation state from a JLD2 file.
"""
function load_state(filepath::String)
    if !isfile(filepath)
        error("State file not found: $filepath")
    end

    state = jldopen(filepath, "r") do file
        file["state"]
    end

    return state
end

"""
    create_bodies_from_state(state, params; load_activity=:keep) -> SimBodies

Create SimBodies from a loaded simulation state.

Arguments:
- `state`: SimulationState loaded from file
- `params`: Parameters for the new simulation (may differ from original)
- `load_activity`: How to handle activity pattern
  - `:keep` - Use the activity pattern from the saved state
  - `:new` - Generate new activity pattern from current params
"""
function create_bodies_from_state(state::SimulationState, params::Parameters;
                                   load_activity::Symbol=:keep)
    L = params.L > 0 ? params.L : state.L
    boundary = CubicBoundary(L, L, L)

    n_particles = length(state.coords)

    # Determine activity pattern
    if load_activity == :keep
        activity = state.activity
    else  # :new
        activity = get_activity_vector(params)
    end

    # Validate particle counts match
    expected_n = get_n_particles(params)
    if n_particles != expected_n
        error("Particle count mismatch: state has $n_particles, params expects $expected_n")
    end

    # Create atoms with restored state
    atoms = [Particle(
        id,
        1.0,                          # mass
        1.0,                          # σ
        1.0,                          # ϵ
        activity[id],                 # active (from chosen activity pattern)
        state.tang_vec[id],           # tang_vec
        state.rold[id],               # rold
        state.image[id]               # image
    ) for id in 1:n_particles]

    neighbor_finder = create_neighbor_finder(state.coords, params.rcut_nf)

    return SimBodies(
        atoms,
        state.coords,
        state.velocities,
        boundary,
        activity,
        neighbor_finder
    )
end
