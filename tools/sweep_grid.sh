#!/usr/bin/env bash
# =============================================================================
# sweep_grid.sh - Run batch simulations sweeping multiple parameters (grid search)
# =============================================================================
#
# Usage:
#   ./tools/sweep_grid.sh <config_file> --<param1> "values" --<param2> "values" [options]
#
# Examples:
#   # Sweep n-active and fact (3x3 = 9 simulations)
#   ./tools/sweep_grid.sh config/single_ring.toml \
#       --n-active "10 30 50" --fact "1.0 3.0 5.0"
#
#   # With parallel execution
#   ./tools/sweep_grid.sh config/single_ring.toml \
#       --n-active "10 20 30" --fact "1.0 2.0" --parallel 4
#
#   # Dry run to preview all combinations
#   ./tools/sweep_grid.sh config/single_ring.toml \
#       --n-active "10 20" --fact "1.0 2.0" --dry-run
#
# =============================================================================

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default values
PARALLEL_JOBS=1
DRY_RUN=false
PREFIX=""
CONFIG_FILE=""

# Arrays to store parameter names and their values
declare -a PARAM_NAMES=()
declare -a PARAM_VALUES=()

# Print usage
usage() {
    echo "Usage: $0 <config_file> --<param1> \"values\" --<param2> \"values\" [options]"
    echo ""
    echo "Options:"
    echo "  --parallel N    Run N simulations in parallel (default: 1)"
    echo "  --dry-run       Print commands without executing"
    echo "  --prefix STR    Prefix for simulation IDs"
    echo "  --help          Show this help message"
    echo ""
    echo "Example:"
    echo "  $0 config/single_ring.toml --n-active \"10 30 50\" --fact \"1.0 3.0 5.0\" --parallel 4"
    exit 1
}

# Parse arguments
parse_args() {
    if [[ $# -lt 3 ]]; then
        usage
    fi

    CONFIG_FILE="$1"
    shift

    if [[ ! -f "$CONFIG_FILE" ]]; then
        echo -e "${RED}Error: Config file not found: $CONFIG_FILE${NC}"
        exit 1
    fi

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --parallel)
                PARALLEL_JOBS="$2"
                shift 2
                ;;
            --dry-run)
                DRY_RUN=true
                shift
                ;;
            --prefix)
                PREFIX="$2"
                shift 2
                ;;
            --help)
                usage
                ;;
            --*)
                PARAM_NAMES+=("$1")
                PARAM_VALUES+=("$2")
                shift 2
                ;;
            *)
                echo -e "${RED}Error: Unknown argument: $1${NC}"
                usage
                ;;
        esac
    done

    if [[ ${#PARAM_NAMES[@]} -eq 0 ]]; then
        echo -e "${RED}Error: At least one parameter is required${NC}"
        usage
    fi
}

# Generate all combinations of parameter values
# Stores results in global COMBINATIONS array
declare -a COMBINATIONS=()

generate_combinations() {
    # Start with a single empty combination
    COMBINATIONS=("")

    for param_idx in "${!PARAM_NAMES[@]}"; do
        local values_str="${PARAM_VALUES[$param_idx]}"
        read -ra values_array <<< "$values_str"

        local -a new_combinations=()
        for combo in "${COMBINATIONS[@]}"; do
            for value in "${values_array[@]}"; do
                if [[ -z "$combo" ]]; then
                    new_combinations+=("$value")
                else
                    new_combinations+=("$combo,$value")
                fi
            done
        done
        COMBINATIONS=("${new_combinations[@]}")
    done
}

# Run a single simulation with given parameter values
run_simulation() {
    local combo_str="$1"
    local sim_idx="$2"

    # Parse combination string back to array
    IFS=',' read -ra values <<< "$combo_str"

    # Build command
    local cmd="julia --project=. scripts/simulate.jl --config $CONFIG_FILE"
    local param_desc=""

    for i in "${!PARAM_NAMES[@]}"; do
        cmd="$cmd ${PARAM_NAMES[$i]} ${values[$i]}"
        if [[ -n "$param_desc" ]]; then
            param_desc="$param_desc, "
        fi
        param_desc="$param_desc${PARAM_NAMES[$i]}=${values[$i]}"
    done

    # Add simulation ID
    local simid=""
    if [[ -n "$PREFIX" ]]; then
        simid="${PREFIX}_${sim_idx}"
    else
        simid="grid_${sim_idx}"
    fi
    cmd="$cmd --simid $simid"

    if [[ "$DRY_RUN" == true ]]; then
        echo -e "${BLUE}[DRY-RUN]${NC} $cmd"
        return 0
    else
        echo -e "${GREEN}[RUNNING]${NC} [$sim_idx] $param_desc"
        if eval "$cmd"; then
            echo -e "${GREEN}[DONE]${NC} [$sim_idx] $param_desc"
        else
            echo -e "${RED}[FAILED]${NC} [$sim_idx] $param_desc"
            return 1
        fi
    fi
}

# Main function
main() {
    parse_args "$@"

    # Generate all combinations
    generate_combinations
    local total=${#COMBINATIONS[@]}

    echo -e "${YELLOW}=== Grid Parameter Sweep ===${NC}"
    echo -e "Config:    $CONFIG_FILE"
    echo -e "Parameters:"
    for i in "${!PARAM_NAMES[@]}"; do
        echo -e "  ${PARAM_NAMES[$i]}: ${PARAM_VALUES[$i]}"
    done
    echo -e "Total:     $total simulations"
    echo -e "Parallel:  $PARALLEL_JOBS jobs"
    if [[ -n "$PREFIX" ]]; then
        echo -e "Prefix:    $PREFIX"
    fi
    echo ""

    if [[ "$DRY_RUN" == true ]]; then
        echo -e "${YELLOW}Dry run mode - commands will be printed but not executed${NC}"
        echo ""
    fi

    # Run simulations
    local running=0
    local completed=0
    local failed=0
    local pids=()

    for i in "${!COMBINATIONS[@]}"; do
        local combo="${COMBINATIONS[$i]}"

        if [[ "$PARALLEL_JOBS" -eq 1 ]]; then
            # Sequential execution
            if run_simulation "$combo" "$i"; then
                completed=$((completed + 1))
            else
                failed=$((failed + 1))
            fi
        else
            # Parallel execution
            run_simulation "$combo" "$i" &
            pids+=($!)
            running=$((running + 1))

            # Wait if we've reached the parallel limit
            if [[ $running -ge $PARALLEL_JOBS ]]; then
                wait "${pids[0]}"
                if [[ $? -eq 0 ]]; then
                    completed=$((completed + 1))
                else
                    failed=$((failed + 1))
                fi
                pids=("${pids[@]:1}")
                running=$((running - 1))
            fi
        fi
    done

    # Wait for remaining parallel jobs
    for pid in "${pids[@]}"; do
        wait "$pid"
        if [[ $? -eq 0 ]]; then
            completed=$((completed + 1))
        else
            failed=$((failed + 1))
        fi
    done

    # Summary
    if [[ "$DRY_RUN" != true ]]; then
        echo ""
        echo -e "${YELLOW}=== Summary ===${NC}"
        echo -e "Completed: ${GREEN}$completed${NC}"
        if [[ $failed -gt 0 ]]; then
            echo -e "Failed:    ${RED}$failed${NC}"
        fi
    fi
}

main "$@"
