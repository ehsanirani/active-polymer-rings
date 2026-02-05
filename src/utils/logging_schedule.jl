export make_logging_schedule

"""
    make_logging_schedule(mode::Symbol, total_steps::Int;
                          fixed_interval::Int=500, n_points::Int=1000)

Generate a vector of step numbers at which to log metrics.

# Modes
- `:fixed` — evenly spaced at `fixed_interval`: `[interval, 2*interval, ..., total_steps]`
- `:logspaced` — logarithmically spaced with `n_points` samples from step 1 to `total_steps`

Returns a sorted `Vector{Int}` of unique step indices.
"""
function make_logging_schedule(mode::Symbol, total_steps::Int;
                               fixed_interval::Int=500, n_points::Int=1000)
    if mode == :fixed
        indices = collect(fixed_interval:fixed_interval:total_steps)
    elseif mode == :logspaced
        raw = 10 .^ range(0, log10(total_steps), length=n_points)
        indices = sort(unique(round.(Int, raw)))
        if isempty(indices) || indices[end] != total_steps
            push!(indices, total_steps)
        end
    else
        error("Unknown logging mode: $mode. Use :fixed or :logspaced")
    end
    return indices
end
