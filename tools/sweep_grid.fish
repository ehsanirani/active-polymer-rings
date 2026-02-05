#!/usr/bin/env fish
# =============================================================================
# sweep_grid.fish - Run batch simulations sweeping multiple parameters (grid search)
# =============================================================================
#
# Usage:
#   ./tools/sweep_grid.fish <config_file> [overrides...] --sweep --<p1> "vals" --sweep --<p2> "vals" [options]
#
# Examples:
#   # Sweep n-active and fact (3x3 = 9 simulations)
#   ./tools/sweep_grid.fish config/single_ring.toml \
#       --sweep --n-active "10 30 50" --sweep --fact "1.0 3.0 5.0"
#
#   # Override some params and sweep others
#   ./tools/sweep_grid.fish config/single_ring.toml \
#       --n-monomers 200 --kangle 5.0 \
#       --sweep --n-active "10 30 50" --sweep --fact "1.0 3.0"
#
#   # With parallel execution
#   ./tools/sweep_grid.fish config/single_ring.toml \
#       --n-monomers 200 \
#       --sweep --n-active "10 20 30" --sweep --fact "1.0 2.0" --parallel 4
#
#   # Dry run to preview all combinations
#   ./tools/sweep_grid.fish config/single_ring.toml \
#       --sweep --n-active "10 20" --sweep --fact "1.0 2.0" --dry-run
#
# =============================================================================

# Default values
set -g PARALLEL_JOBS 1
set -g DRY_RUN false
set -g PREFIX ""
set -g CONFIG_FILE ""
set -g FIXED_OVERRIDES ""

# Arrays to store parameter names and their values
set -g PARAM_NAMES
set -g PARAM_VALUES

# Generated combinations
set -g COMBINATIONS

# Print usage
function usage
    echo "Usage: "(status basename)" <config_file> [overrides...] --sweep --<p1> \"vals\" --sweep --<p2> \"vals\" [options]"
    echo ""
    echo "Options:"
    echo "  --sweep         Mark the next parameter as one to sweep"
    echo "  --parallel N    Run N simulations in parallel (default: 1)"
    echo "  --dry-run       Print commands without executing"
    echo "  --prefix STR    Prefix for simulation IDs"
    echo "  --help          Show this help message"
    echo ""
    echo "Example:"
    echo "  "(status basename)" config/single_ring.toml --n-monomers 200 --sweep --n-active \"10 30\" --sweep --fact \"1.0 3.0\" --parallel 4"
    exit 1
end

# Parse arguments
function parse_args
    if test (count $argv) -lt 2
        usage
    end

    set -g CONFIG_FILE $argv[1]
    set -e argv[1]

    if not test -f $CONFIG_FILE
        set_color red
        echo "Error: Config file not found: $CONFIG_FILE"
        set_color normal
        exit 1
    end

    set -l in_sweep false

    while test (count $argv) -gt 0
        switch $argv[1]
            case --parallel
                set -g PARALLEL_JOBS $argv[2]
                set -e argv[1..2]
            case --dry-run
                set -g DRY_RUN true
                set -e argv[1]
            case --prefix
                set -g PREFIX $argv[2]
                set -e argv[1..2]
            case --help
                usage
            case --sweep
                set in_sweep true
                set -e argv[1]
            case '--*'
                if test "$in_sweep" = true
                    # This is a sweep parameter
                    set -a PARAM_NAMES $argv[1]
                    set -a PARAM_VALUES $argv[2]
                    set in_sweep false
                    set -e argv[1..2]
                else
                    # This is a fixed override
                    set -g FIXED_OVERRIDES "$FIXED_OVERRIDES $argv[1] $argv[2]"
                    set -e argv[1..2]
                end
            case '*'
                set_color red
                echo "Error: Unknown argument: $argv[1]"
                set_color normal
                usage
        end
    end

    if test (count $PARAM_NAMES) -eq 0
        set_color red
        echo "Error: At least one sweep parameter is required (use --sweep --<param> \"values\")"
        set_color normal
        usage
    end
end

# Generate all combinations of parameter values
function generate_combinations
    set -g COMBINATIONS

    # Start with empty combinations
    set -l current_combos ""

    for param_idx in (seq (count $PARAM_NAMES))
        set -l values (string split " " $PARAM_VALUES[$param_idx])
        set -l new_combos

        if test -z "$current_combos"
            # First parameter - just use its values
            for v in $values
                set -a new_combos $v
            end
        else
            # Subsequent parameters - combine with existing
            for combo in $current_combos
                for v in $values
                    set -a new_combos "$combo,$v"
                end
            end
        end

        set current_combos $new_combos
    end

    set -g COMBINATIONS $current_combos
end

# Run a single simulation with given parameter values
function run_simulation
    set -l combo_str $argv[1]
    set -l sim_idx $argv[2]

    # Parse combination string back to array
    set -l values (string split "," $combo_str)

    # Build command
    set -l cmd "julia --project=. scripts/simulate.jl --config $CONFIG_FILE$FIXED_OVERRIDES"
    set -l param_desc ""

    for i in (seq (count $PARAM_NAMES))
        set cmd "$cmd $PARAM_NAMES[$i] $values[$i]"
        if test -n "$param_desc"
            set param_desc "$param_desc, "
        end
        set param_desc "$param_desc$PARAM_NAMES[$i]=$values[$i]"
    end

    # Add simulation ID
    set -l simid ""
    if test -n "$PREFIX"
        set simid {$PREFIX}_{$sim_idx}
    else
        set simid grid_{$sim_idx}
    end
    set cmd "$cmd --simid $simid"

    if test "$DRY_RUN" = true
        set_color blue
        echo -n "[DRY-RUN] "
        set_color normal
        echo $cmd
        return 0
    else
        set_color green
        echo -n "[RUNNING] "
        set_color normal
        echo "[$sim_idx] $param_desc"

        if eval $cmd
            set_color green
            echo -n "[DONE] "
            set_color normal
            echo "[$sim_idx] $param_desc"
            return 0
        else
            set_color red
            echo -n "[FAILED] "
            set_color normal
            echo "[$sim_idx] $param_desc"
            return 1
        end
    end
end

# Main function
function main
    parse_args $argv

    # Generate all combinations
    generate_combinations
    set -l total (count $COMBINATIONS)

    set_color yellow
    echo "=== Grid Parameter Sweep ==="
    set_color normal
    echo "Config:    $CONFIG_FILE"
    if test -n "$FIXED_OVERRIDES"
        echo "Overrides:$FIXED_OVERRIDES"
    end
    echo "Sweep parameters:"
    for i in (seq (count $PARAM_NAMES))
        echo "  $PARAM_NAMES[$i]: $PARAM_VALUES[$i]"
    end
    echo "Total:     $total simulations"
    echo "Parallel:  $PARALLEL_JOBS jobs"
    if test -n "$PREFIX"
        echo "Prefix:    $PREFIX"
    end
    echo ""

    if test "$DRY_RUN" = true
        set_color yellow
        echo "Dry run mode - commands will be printed but not executed"
        set_color normal
        echo ""
    end

    # Run simulations
    set -l completed 0
    set -l failed 0
    set -l jobs

    for i in (seq (count $COMBINATIONS))
        set -l combo $COMBINATIONS[$i]
        set -l idx (math $i - 1)

        if test $PARALLEL_JOBS -eq 1
            # Sequential execution
            if run_simulation $combo $idx
                set completed (math $completed + 1)
            else
                set failed (math $failed + 1)
            end
        else
            # Parallel execution
            run_simulation $combo $idx &
            set -a jobs $last_pid

            # Wait if we've reached the parallel limit
            if test (count $jobs) -ge $PARALLEL_JOBS
                wait $jobs[1]
                if test $status -eq 0
                    set completed (math $completed + 1)
                else
                    set failed (math $failed + 1)
                end
                set -e jobs[1]
            end
        end
    end

    # Wait for remaining parallel jobs
    for pid in $jobs
        wait $pid
        if test $status -eq 0
            set completed (math $completed + 1)
        else
            set failed (math $failed + 1)
        end
    end

    # Summary
    if test "$DRY_RUN" != true
        echo ""
        set_color yellow
        echo "=== Summary ==="
        set_color normal
        set_color green
        echo "Completed: $completed"
        set_color normal
        if test $failed -gt 0
            set_color red
            echo "Failed:    $failed"
            set_color normal
        end
    end
end

main $argv
