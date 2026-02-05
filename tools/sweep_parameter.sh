#!/usr/bin/env bash
# =============================================================================
# sweep_parameter.sh - Run batch simulations sweeping a single parameter
# =============================================================================
#
# Usage:
#   ./tools/sweep_parameter.sh <config_file> --<param> "value1 value2 ..." [options]
#
# Examples:
#   # Sweep active force
#   ./tools/sweep_parameter.sh config/single_ring.toml --fact "1.0 2.0 3.0 4.0 5.0"
#
#   # Sweep with 4 parallel jobs
#   ./tools/sweep_parameter.sh config/single_ring.toml --n-active "10 20 30 40" --parallel 4
#
#   # Dry run to preview commands
#   ./tools/sweep_parameter.sh config/single_ring.toml --fact "1.0 2.0" --dry-run
#
#   # Add prefix to simulation IDs
#   ./tools/sweep_parameter.sh config/single_ring.toml --fact "1.0 2.0" --prefix sweep01
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
PARAM_NAME=""
PARAM_VALUES=""

# Print usage
usage() {
    echo "Usage: $0 <config_file> --<param> \"value1 value2 ...\" [options]"
    echo ""
    echo "Options:"
    echo "  --parallel N    Run N simulations in parallel (default: 1)"
    echo "  --dry-run       Print commands without executing"
    echo "  --prefix STR    Prefix for simulation IDs"
    echo "  --help          Show this help message"
    echo ""
    echo "Example:"
    echo "  $0 config/single_ring.toml --fact \"1.0 2.0 3.0\" --parallel 4"
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
                if [[ -z "$PARAM_NAME" ]]; then
                    PARAM_NAME="$1"
                    PARAM_VALUES="$2"
                    shift 2
                else
                    echo -e "${RED}Error: Only one parameter can be swept at a time${NC}"
                    echo "Use sweep_grid.sh for multiple parameters"
                    exit 1
                fi
                ;;
            *)
                echo -e "${RED}Error: Unknown argument: $1${NC}"
                usage
                ;;
        esac
    done

    if [[ -z "$PARAM_NAME" || -z "$PARAM_VALUES" ]]; then
        echo -e "${RED}Error: Parameter name and values are required${NC}"
        usage
    fi
}

# Run a single simulation
run_simulation() {
    local value="$1"
    local simid="$2"
    local cmd="julia --project=. scripts/simulate.jl --config $CONFIG_FILE $PARAM_NAME $value"

    if [[ -n "$simid" ]]; then
        cmd="$cmd --simid $simid"
    fi

    if [[ "$DRY_RUN" == true ]]; then
        echo -e "${BLUE}[DRY-RUN]${NC} $cmd"
        return 0
    else
        echo -e "${GREEN}[RUNNING]${NC} $PARAM_NAME=$value"
        if eval "$cmd"; then
            echo -e "${GREEN}[DONE]${NC} $PARAM_NAME=$value"
        else
            echo -e "${RED}[FAILED]${NC} $PARAM_NAME=$value"
            return 1
        fi
    fi
}

# Main function
main() {
    parse_args "$@"

    # Convert values string to array
    read -ra VALUES_ARRAY <<< "$PARAM_VALUES"
    local total=${#VALUES_ARRAY[@]}

    echo -e "${YELLOW}=== Parameter Sweep ===${NC}"
    echo -e "Config:    $CONFIG_FILE"
    echo -e "Parameter: $PARAM_NAME"
    echo -e "Values:    ${VALUES_ARRAY[*]}"
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

    for i in "${!VALUES_ARRAY[@]}"; do
        local value="${VALUES_ARRAY[$i]}"
        local simid=""

        if [[ -n "$PREFIX" ]]; then
            simid="${PREFIX}_${i}"
        fi

        if [[ "$PARALLEL_JOBS" -eq 1 ]]; then
            # Sequential execution
            if run_simulation "$value" "$simid"; then
                completed=$((completed + 1))
            else
                failed=$((failed + 1))
            fi
        else
            # Parallel execution
            run_simulation "$value" "$simid" &
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
