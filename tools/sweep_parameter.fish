#!/usr/bin/env fish
# =============================================================================
# sweep_parameter.fish - Run batch simulations sweeping a single parameter
# =============================================================================
#
# Usage:
#   ./tools/sweep_parameter.fish <config_file> [overrides...] [--sweep --<param> "values"] [options]
#
# Examples:
#   # Sweep active force
#   ./tools/sweep_parameter.fish config/single_ring.toml \
#       --sweep --fact "1.0 2.0 3.0 4.0 5.0"
#
#   # Override some params and sweep another
#   ./tools/sweep_parameter.fish config/single_ring.toml \
#       --n-monomers 200 --kangle 5.0 \
#       --sweep --fact "1.0 2.0 3.0"
#
#   # Run 10 replicates with same parameters
#   ./tools/sweep_parameter.fish config/single_ring.toml --runs 10
#
#   # Sweep parameter with 5 replicates each (3 values × 5 runs = 15 sims)
#   ./tools/sweep_parameter.fish config/single_ring.toml \
#       --sweep --fact "1.0 2.0 3.0" --runs 5
#
#   # With parallel execution
#   ./tools/sweep_parameter.fish config/single_ring.toml \
#       --sweep --n-active "10 20 30 40" --runs 3 --parallel 4
#
#   # Dry run to preview commands
#   ./tools/sweep_parameter.fish config/single_ring.toml \
#       --sweep --fact "1.0 2.0" --runs 3 --dry-run
#
# =============================================================================

# Default values
set -g PARALLEL_JOBS 1
set -g DRY_RUN false
set -g PREFIX "run"
set -g CONFIG_FILE ""
set -g PARAM_NAME ""
set -g PARAM_VALUES
set -g FIXED_OVERRIDES ""
set -g NUM_RUNS 1

# Print usage
function usage
    echo "Usage: "(status basename)" <config_file> [overrides...] [--sweep --<param> \"values\"] [options]"
    echo ""
    echo "Options:"
    echo "  --sweep         Mark the next parameter as the one to sweep"
    echo "  --runs N        Number of replicates per parameter combination (default: 1)"
    echo "  --parallel N    Run N simulations in parallel (default: 1)"
    echo "  --dry-run       Print commands without executing"
    echo "  --prefix STR    Prefix for simulation IDs (default: run)"
    echo "  --help          Show this help message"
    echo ""
    echo "Examples:"
    echo "  # Sweep a parameter"
    echo "  "(status basename)" config/single_ring.toml --sweep --fact \"1.0 2.0 3.0\""
    echo ""
    echo "  # Run 10 replicates with same parameters"
    echo "  "(status basename)" config/single_ring.toml --runs 10"
    echo ""
    echo "  # Sweep with 5 replicates each"
    echo "  "(status basename)" config/single_ring.toml --sweep --fact \"1.0 2.0 3.0\" --runs 5"
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
            case --runs
                set -g NUM_RUNS $argv[2]
                set -e argv[1..2]
            case --help
                usage
            case --sweep
                set in_sweep true
                set -e argv[1]
            case '--*'
                if test "$in_sweep" = true
                    # This is the sweep parameter
                    if test -n "$PARAM_NAME"
                        set_color red
                        echo "Error: Only one parameter can be swept at a time"
                        set_color normal
                        echo "Use sweep_grid.fish for multiple parameters"
                        exit 1
                    end
                    set -g PARAM_NAME $argv[1]
                    set -g PARAM_VALUES (string split " " $argv[2])
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

    # If no sweep parameter, we're just doing replicates
    if test -z "$PARAM_NAME"
        if test $NUM_RUNS -lt 2
            set_color red
            echo "Error: Must specify --sweep or --runs N (with N >= 2)"
            set_color normal
            usage
        end
    end
end

# Run a single simulation
function run_simulation
    set -l value $argv[1]
    set -l simid $argv[2]
    set -l cmd "julia --project=. scripts/simulate.jl --config $CONFIG_FILE$FIXED_OVERRIDES"

    if test -n "$value"
        set cmd "$cmd $PARAM_NAME $value"
    end

    if test -n "$simid"
        set cmd "$cmd --simid $simid"
    end

    if test "$DRY_RUN" = true
        set_color blue
        echo -n "[DRY-RUN] "
        set_color normal
        echo $cmd
        return 0
    else
        set -l desc ""
        if test -n "$value"
            set desc "$PARAM_NAME=$value, "
        end
        set desc "{$desc}simid=$simid"
        set_color green
        echo -n "[RUNNING] "
        set_color normal
        echo $desc

        if eval $cmd
            set_color green
            echo -n "[DONE] "
            set_color normal
            echo $desc
            return 0
        else
            set_color red
            echo -n "[FAILED] "
            set_color normal
            echo $desc
            return 1
        end
    end
end

# Main function
function main
    parse_args $argv

    # Determine what we're doing
    set -l has_sweep false
    set -l n_values 1

    if test -n "$PARAM_NAME"
        set has_sweep true
        set n_values (count $PARAM_VALUES)
    end

    set -l total (math "$n_values * $NUM_RUNS")

    set_color yellow
    echo "=== Parameter Sweep ==="
    set_color normal
    echo "Config:    $CONFIG_FILE"
    if test -n "$FIXED_OVERRIDES"
        echo "Overrides:$FIXED_OVERRIDES"
    end
    if test "$has_sweep" = true
        echo "Sweep:     $PARAM_NAME = $PARAM_VALUES"
    end
    if test $NUM_RUNS -gt 1
        echo "Runs:      $NUM_RUNS replicates per combination"
    end
    echo "Total:     $total simulations"
    echo "Parallel:  $PARALLEL_JOBS jobs"
    echo "Prefix:    $PREFIX"
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

    # If no sweep, just do runs
    if test "$has_sweep" = false
        for run in (seq 0 (math "$NUM_RUNS - 1"))
            set -l simid "$PREFIX"_"$run"

            if test $PARALLEL_JOBS -eq 1
                if run_simulation "" $simid
                    set completed (math $completed + 1)
                else
                    set failed (math $failed + 1)
                end
            else
                run_simulation "" $simid &
                set -a jobs $last_pid

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
    else
        # Sweep with optional runs
        for i in (seq (count $PARAM_VALUES))
            set -l value $PARAM_VALUES[$i]
            set -l idx (math $i - 1)

            for run in (seq 0 (math "$NUM_RUNS - 1"))
                set -l simid ""
                if test $NUM_RUNS -gt 1
                    set simid "$PREFIX"_"$idx"_"$run"
                else
                    set simid "$PREFIX"_"$idx"
                end

                if test $PARALLEL_JOBS -eq 1
                    if run_simulation $value $simid
                        set completed (math $completed + 1)
                    else
                        set failed (math $failed + 1)
                    end
                else
                    run_simulation $value $simid &
                    set -a jobs $last_pid

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
