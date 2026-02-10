using Molly
using Printf

export CheckpointLogger, get_latest_checkpoint

"""
A Molly logger that saves simulation checkpoints at regular intervals.

Checkpoints are saved with filenames like:
  checkpoint_thermal_step100000.jld2
  checkpoint_active_step200000.jld2

Only the last `keep_last` checkpoints are retained to save disk space.
"""
mutable struct CheckpointLogger
    interval::Int64              # Save every N steps
    directory::String            # Directory for checkpoint files
    params::Parameters
    phase::Symbol                # :thermalization or :active
    keep_last::Int               # Number of checkpoints to keep
    total_steps::Int64           # Total steps planned for this phase
    sim_bodies::SimBodies        # Reference to sim_bodies for saving
    saved_files::Vector{String}  # Track saved checkpoint files
end

"""
    CheckpointLogger(interval, directory, params, phase, total_steps, sim_bodies; keep_last=2)

Create a checkpoint logger.

Arguments:
- `interval`: Save checkpoint every N steps (0 = disabled)
- `directory`: Directory to save checkpoint files
- `params`: Simulation parameters
- `phase`: Current phase (:thermalization or :active)
- `total_steps`: Total steps planned for this phase
- `sim_bodies`: SimBodies reference for state saving
- `keep_last`: Number of checkpoints to keep (default: 2)
"""
function CheckpointLogger(interval::Int64, directory::String, params::Parameters,
                          phase::Symbol, total_steps::Int64, sim_bodies::SimBodies;
                          keep_last::Int=2)
    # Create directory if it doesn't exist
    if interval > 0
        mkpath(directory)
    end

    CheckpointLogger(
        interval,
        directory,
        params,
        phase,
        keep_last,
        total_steps,
        sim_bodies,
        String[]
    )
end

function Molly.log_property!(logger::CheckpointLogger, sys, buffers, neighbors=nothing,
                              step_n::Integer=0; n_threads=1, kwargs...)
    # Skip if checkpointing is disabled
    if logger.interval <= 0
        return
    end

    # Only checkpoint at specified intervals
    if step_n % logger.interval != 0 || step_n == 0
        return
    end

    # Generate checkpoint filename
    phase_str = string(logger.phase)
    filename = @sprintf("checkpoint_%s_step%d.jld2", phase_str, step_n)
    filepath = joinpath(logger.directory, filename)

    # Save state
    save_state(filepath, logger.sim_bodies, sys, logger.params;
               phase=logger.phase, step=Int64(step_n), total_steps=logger.total_steps)

    # Track this checkpoint
    push!(logger.saved_files, filepath)

    # Clean up old checkpoints (keep only the last N)
    cleanup_old_checkpoints!(logger)
end

"""
Remove old checkpoint files, keeping only the last `keep_last` files.
"""
function cleanup_old_checkpoints!(logger::CheckpointLogger)
    while length(logger.saved_files) > logger.keep_last
        old_file = popfirst!(logger.saved_files)
        if isfile(old_file)
            rm(old_file)
        end
    end
end

"""
    get_latest_checkpoint(directory; phase=nothing) -> Union{String, Nothing}

Find the latest checkpoint file in a directory.

Arguments:
- `directory`: Directory to search for checkpoints
- `phase`: Optional phase filter (:thermalization or :active)

Returns the path to the latest checkpoint, or nothing if none found.
"""
function get_latest_checkpoint(directory::String; phase::Union{Symbol,Nothing}=nothing)
    if !isdir(directory)
        return nothing
    end

    # Find all checkpoint files
    pattern = if phase === nothing
        r"checkpoint_(thermal|active)_step(\d+)\.jld2"
    elseif phase == :thermalization
        r"checkpoint_thermal_step(\d+)\.jld2"
    else  # :active
        r"checkpoint_active_step(\d+)\.jld2"
    end

    checkpoints = String[]
    for file in readdir(directory)
        if occursin(pattern, file)
            push!(checkpoints, joinpath(directory, file))
        end
    end

    if isempty(checkpoints)
        return nothing
    end

    # Sort by step number (extract from filename) and return latest
    function extract_step(filepath)
        m = match(r"step(\d+)\.jld2", filepath)
        return m === nothing ? 0 : parse(Int, m.captures[1])
    end

    sort!(checkpoints, by=extract_step, rev=true)
    return first(checkpoints)
end
