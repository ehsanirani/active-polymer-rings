#!/usr/bin/env fish
# =============================================================================
# sweep_parameter.fish - Run batch simulations sweeping a single parameter
# =============================================================================
#
# Usage:
#   ./tools/sweep_parameter.fish <config_file> --<param> "value1 value2 ..." [options]
#
# Examples:
#   # Sweep active force
#   ./tools/sweep_parameter.fish config/single_ring.toml --fact "1.0 2.0 3.0 4.0 5.0"
#
#   # Sweep with 4 parallel jobs
#   ./tools/sweep_parameter.fish config/single_ring.toml --n-active "10 20 30 40" --parallel 4
#
#   # Dry run to preview commands
#   ./tools/sweep_parameter.fish config/single_ring.toml --fact "1.0 2.0" --dry-run
#
#   # Add prefix to simulation IDs
#   ./tools/sweep_parameter.fish config/single_ring.toml --fact "1.0 2.0" --prefix sweep01
#
# =============================================================================

# Default values
set -g PARALLEL_JOBS 1
set -g DRY_RUN false
set -g PREFIX ""
set -g CONFIG_FILE ""
set -g PARAM_NAME ""
set -g PARAM_VALUES

# Print usage
function usage
    echo "Usage: "(status basename)" <config_file> --<param> \"value1 value2 ...\" [options]"
    echo ""
    echo "Options:"
    echo "  --parallel N    Run N simulations in parallel (default: 1)"
    echo "  --dry-run       Print commands without executing"
    echo "  --prefix STR    Prefix for simulation IDs"
    echo "  --help          Show this help message"
    echo ""
    echo "Example:"
    echo "  "(status basename)" config/single_ring.toml --fact \"1.0 2.0 3.0\" --parallel 4"
    exit 1
end

# Parse arguments
function parse_args
    if test (count $argv) -lt 3
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
            case '--*'
                if test -z "$PARAM_NAME"
                    set -g PARAM_NAME $argv[1]
                    set -g PARAM_VALUES (string split " " $argv[2])
                    set -e argv[1..2]
                else
                    set_color red
                    echo "Error: Only one parameter can be swept at a time"
                    set_color normal
                    echo "Use sweep_grid.fish for multiple parameters"
                    exit 1
                end
            case '*'
                set_color red
                echo "Error: Unknown argument: $argv[1]"
                set_color normal
                usage
        end
    end

    if test -z "$PARAM_NAME" -o (count $PARAM_VALUES) -eq 0
        set_color red
        echo "Error: Parameter name and values are required"
        set_color normal
        usage
    end
end

# Run a single simulation
function run_simulation
    set -l value $argv[1]
    set -l simid $argv[2]
    set -l cmd "julia --project=. scripts/simulate.jl --config $CONFIG_FILE $PARAM_NAME $value"

    if test -n "$simid"
        set cmd "$cmd --simid $simid"
    end

    if test "$DRY_RUN" = true
        set_color blue
        echo -n "[DRY-RUN] "
        set_color normal
        echo $cmd
    else
        set_color green
        echo -n "[RUNNING] "
        set_color normal
        echo "$PARAM_NAME=$value"

        if eval $cmd
            set_color green
            echo -n "[DONE] "
            set_color normal
            echo "$PARAM_NAME=$value"
            return 0
        else
            set_color red
            echo -n "[FAILED] "
            set_color normal
            echo "$PARAM_NAME=$value"
            return 1
        end
    end
end

# Main function
function main
    parse_args $argv

    set -l total (count $PARAM_VALUES)

    set_color yellow
    echo "=== Parameter Sweep ==="
    set_color normal
    echo "Config:    $CONFIG_FILE"
    echo "Parameter: $PARAM_NAME"
    echo "Values:    $PARAM_VALUES"
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

    for i in (seq (count $PARAM_VALUES))
        set -l value $PARAM_VALUES[$i]
        set -l simid ""

        if test -n "$PREFIX"
            set simid {$PREFIX}_{$i}
        end

        if test $PARALLEL_JOBS -eq 1
            # Sequential execution
            if run_simulation $value $simid
                set completed (math $completed + 1)
            else
                set failed (math $failed + 1)
            end
        else
            # Parallel execution using fish's job control
            run_simulation $value $simid &
            set jobs $jobs $last_pid

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
