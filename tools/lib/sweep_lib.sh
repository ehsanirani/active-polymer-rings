#!/usr/bin/env bash
# =============================================================================
# sweep_lib.sh - Shared library for sweep scripts
# =============================================================================
#
# This library provides common functions for parameter sweep scripts:
# - Color definitions and logging utilities
# - Simulation runner with parallel execution support
# - Base state management (passive/active equilibration)
# - Range notation expansion (e.g., "0:25:200")
# - Cartesian product generation for grid sweeps
#
# Usage:
#   source "$(dirname "${BASH_SOURCE[0]}")/lib/sweep_lib.sh"
#
# =============================================================================

# Prevent multiple sourcing
[[ -n "${_SWEEP_LIB_LOADED:-}" ]] && return 0
_SWEEP_LIB_LOADED=1

# =============================================================================
# Color definitions
# =============================================================================

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# =============================================================================
# Logging utilities
# =============================================================================

log_info() {
    echo -e "${GREEN}[INFO]${NC} $*"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $*"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $*"
}

log_exists() {
    echo -e "${GREEN}[EXISTS]${NC} $*"
}

log_skip() {
    echo -e "${BLUE}[SKIP]${NC} $*"
}

log_running() {
    echo -e "${GREEN}[RUNNING]${NC} $*"
}

log_done() {
    echo -e "${GREEN}[DONE]${NC} $*"
}

log_failed() {
    echo -e "${RED}[FAILED]${NC} $*"
}

log_dry_run() {
    echo -e "${BLUE}[DRY-RUN]${NC} $*"
}

# =============================================================================
# Header and summary utilities
# =============================================================================

print_header() {
    local title="$1"
    echo -e "${YELLOW}=== $title ===${NC}"
}

print_summary() {
    local completed="$1"
    local failed="$2"

    echo ""
    print_header "Summary"
    echo -e "Completed: ${GREEN}$completed${NC}"
    if [[ $failed -gt 0 ]]; then
        echo -e "Failed:    ${RED}$failed${NC}"
    fi
}

# =============================================================================
# Value expansion utilities
# =============================================================================

# Expand range notation (start:step:end) to space-separated values
# Examples:
#   "0:25:100" -> "0 25 50 75 100"
#   "1.0:0.5:3.0" -> "1.0 1.5 2.0 2.5 3.0"
#   "10 20 30" -> "10 20 30" (passthrough)
expand_values() {
    local input="$1"

    # Check if it's range notation (contains exactly 2 colons)
    if [[ "$input" =~ ^([0-9.]+):([0-9.]+):([0-9.]+)$ ]]; then
        local start="${BASH_REMATCH[1]}"
        local step="${BASH_REMATCH[2]}"
        local end="${BASH_REMATCH[3]}"

        # Use awk for floating point support
        awk -v s="$start" -v step="$step" -v e="$end" 'BEGIN {
            for (i = s; i <= e + step/1000; i += step) {
                if (i == int(i)) printf "%d ", int(i)
                else printf "%g ", i
            }
        }'
    else
        # Return as-is (space-separated values)
        echo "$input"
    fi
}

# =============================================================================
# Cartesian product generation
# =============================================================================

# Generate all combinations of parameter values
# Arguments:
#   PARAM_NAMES array (global)
#   PARAM_VALUES array (global)
# Sets:
#   COMBINATIONS array (global) - comma-separated value strings
generate_cartesian_product() {
    # Start with a single empty combination
    COMBINATIONS=("")

    for param_idx in "${!PARAM_NAMES[@]}"; do
        local values_str="${PARAM_VALUES[$param_idx]}"
        # Expand range notation if present
        values_str="$(expand_values "$values_str")"
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

# =============================================================================
# Base state management
# =============================================================================

# Compute base state path for a given type
# Arguments:
#   $1 - type: "passive", "active", or file path
#   $2 - n_monomers
#   $3 - thermal_steps (optional, default: 2000000)
#   $4 - data_dir (optional, default: _data)
#   $5 - run_suffix (optional, e.g., "_run0", "_run1")
# Output:
#   Path to base state file
get_base_state_path() {
    local type="$1"
    local n_monomers="$2"
    local thermal_steps="${3:-2000000}"
    local data_dir="${4:-_data}"
    local run_suffix="${5:-}"

    case "$type" in
        passive|active)
            local thermal_suffix
            if [[ "$thermal_steps" -ge 1000000 ]]; then
                thermal_suffix="$((thermal_steps / 1000000))M"
            else
                thermal_suffix="$thermal_steps"
            fi
            echo "${data_dir}/base_states/${type}_N${n_monomers}_thermal${thermal_suffix}${run_suffix}.jld2"
            ;;
        *)
            # Assume it's a file path
            echo "$type"
            ;;
    esac
}

# Ensure base state exists, creating it if necessary
# Arguments:
#   $1 - base_state_type: "passive", "active", or file path
#   $2 - config_file
#   $3 - n_monomers
#   $4 - thermal_steps (optional)
#   $5 - fact (optional, default: 5.0)
#   $6 - dry_run: true/false
#   $7 - data_dir (optional, default: _data)
#   $8 - run_suffix (optional, e.g., "_run0", "_run1")
# Output:
#   Sets BASE_STATE_PATH global variable
# Returns:
#   0 on success, 1 on failure
ensure_base_state() {
    local base_state_type="$1"
    local config_file="$2"
    local n_monomers="$3"
    local thermal_steps="${4:-2000000}"
    local fact="${5:-5.0}"
    local dry_run="${6:-false}"
    local data_dir="${7:-_data}"
    local run_suffix="${8:-}"

    # Compute path
    BASE_STATE_PATH="$(get_base_state_path "$base_state_type" "$n_monomers" "$thermal_steps" "$data_dir" "$run_suffix")"

    # If it's a custom path, just validate it exists
    if [[ "$base_state_type" != "passive" && "$base_state_type" != "active" ]]; then
        if [[ ! -f "$BASE_STATE_PATH" ]]; then
            log_error "Base state file not found: $BASE_STATE_PATH"
            return 1
        fi
        log_exists "Using base state: $BASE_STATE_PATH"
        return 0
    fi

    # Check if file already exists
    if [[ -f "$BASE_STATE_PATH" ]]; then
        log_exists "Using base state: $BASE_STATE_PATH"
        return 0
    fi

    # Create directory if needed
    local base_dir
    base_dir="$(dirname "$BASE_STATE_PATH")"
    if [[ ! -d "$base_dir" ]]; then
        if [[ "$dry_run" == true ]]; then
            log_dry_run "mkdir -p $base_dir"
        else
            mkdir -p "$base_dir"
        fi
    fi

    # Determine n_active for equilibration
    local n_active
    case "$base_state_type" in
        passive) n_active=0 ;;
        active)  n_active="$n_monomers" ;;
    esac

    # Build command
    local cmd="julia --project=. scripts/simulate.jl"
    if [[ -n "$config_file" ]]; then
        cmd="$cmd --config $config_file"
    fi
    cmd="$cmd --data-dir $data_dir"
    cmd="$cmd --n-monomers $n_monomers"
    cmd="$cmd --activity-pattern block"
    cmd="$cmd --n-active $n_active"
    cmd="$cmd --thermal-steps $thermal_steps"
    cmd="$cmd --n-steps 0"
    cmd="$cmd --fact $fact"
    cmd="$cmd --save-state $BASE_STATE_PATH"
    cmd="$cmd --simid base_state_creation"

    if [[ "$dry_run" == true ]]; then
        log_dry_run "Creating base state ($base_state_type):"
        log_dry_run "$cmd"
        return 0
    fi

    print_header "Creating Base State ($base_state_type)"
    echo -e "Output: $BASE_STATE_PATH"
    echo -e "Running $base_state_type equilibration (n-active=$n_active, thermal-steps=$thermal_steps)..."
    echo ""

    if eval "$cmd"; then
        log_done "Base state created: $BASE_STATE_PATH"
        return 0
    else
        log_failed "Failed to create base state"
        return 1
    fi
}

# Ensure base states exist for each run, creating them if necessary
# Arguments:
#   $1 - base_state_type: "passive", "active", or file path
#   $2 - config_file
#   $3 - n_monomers
#   $4 - thermal_steps
#   $5 - fact
#   $6 - dry_run: true/false
#   $7 - data_dir
#   $8 - run_start
#   $9 - num_runs
# Output:
#   Sets BASE_STATE_PATHS global array
# Returns:
#   0 on success, 1 on failure
ensure_base_states_per_run() {
    local base_state_type="$1"
    local config_file="$2"
    local n_monomers="$3"
    local thermal_steps="$4"
    local fact="$5"
    local dry_run="$6"
    local data_dir="$7"
    local run_start="$8"
    local num_runs="$9"

    # Initialize global array
    declare -ga BASE_STATE_PATHS=()

    for run in $(seq "$run_start" $((run_start + num_runs - 1))); do
        local run_suffix="_run${run}"
        if ! ensure_base_state "$base_state_type" "$config_file" "$n_monomers" \
            "$thermal_steps" "$fact" "$dry_run" "$data_dir" "$run_suffix"; then
            return 1
        fi
        BASE_STATE_PATHS+=("$BASE_STATE_PATH")
    done

    return 0
}

# =============================================================================
# Simulation execution
# =============================================================================

# Build simulation command from parameters
# Arguments:
#   $1 - config_file
#   $2 - fixed_overrides
#   $3 - param_args
#   $4 - simid
#   $5 - base_state_path (optional)
#   $6 - data_dir (optional, default: _data)
# Output:
#   Command string
build_simulation_command() {
    local config_file="$1"
    local fixed_overrides="$2"
    local param_args="$3"
    local simid="$4"
    local base_state_path="${5:-}"
    local data_dir="${6:-_data}"

    local cmd="julia --project=. scripts/simulate.jl"

    # Add config file if specified
    if [[ -n "$config_file" ]]; then
        cmd="$cmd --config $config_file"
    fi

    # Add data directory
    cmd="$cmd --data-dir $data_dir"

    # Add fixed overrides
    if [[ -n "$fixed_overrides" ]]; then
        cmd="$cmd$fixed_overrides"
    fi

    # Add sweep parameters
    if [[ -n "$param_args" ]]; then
        cmd="$cmd $param_args"
    fi

    # Add base state loading if specified
    if [[ -n "$base_state_path" ]]; then
        cmd="$cmd --thermal-steps 0"
        cmd="$cmd --load-state $base_state_path"
        cmd="$cmd --load-activity new"
    fi

    # Add simulation ID
    cmd="$cmd --simid $simid"

    echo "$cmd"
}

# Run a single simulation
# Arguments:
#   $1 - command to run
#   $2 - description
#   $3 - dry_run: true/false
# Returns:
#   0 on success, 1 on failure
run_simulation() {
    local cmd="$1"
    local desc="$2"
    local dry_run="${3:-false}"

    if [[ "$dry_run" == true ]]; then
        log_dry_run "$cmd"
        return 0
    fi

    log_running "$desc"
    if eval "$cmd"; then
        log_done "$desc"
        return 0
    else
        log_failed "$desc"
        return 1
    fi
}

# =============================================================================
# Parallel execution engine
# =============================================================================

# Run commands in parallel with job limiting
# Arguments:
#   $1 - parallel_jobs: number of parallel jobs
#   $2 - dry_run: true/false
#   Commands are read from global PENDING_COMMANDS and PENDING_DESCRIPTIONS arrays
# Returns:
#   Sets COMPLETED and FAILED counts
run_parallel() {
    local parallel_jobs="$1"
    local dry_run="${2:-false}"

    COMPLETED=0
    FAILED=0

    local running=0
    local pids=()

    for i in "${!PENDING_COMMANDS[@]}"; do
        local cmd="${PENDING_COMMANDS[$i]}"
        local desc="${PENDING_DESCRIPTIONS[$i]}"

        if [[ "$parallel_jobs" -eq 1 ]]; then
            # Sequential execution
            if run_simulation "$cmd" "$desc" "$dry_run"; then
                COMPLETED=$((COMPLETED + 1))
            else
                FAILED=$((FAILED + 1))
            fi
        else
            # Parallel execution
            run_simulation "$cmd" "$desc" "$dry_run" &
            pids+=($!)
            running=$((running + 1))

            # Wait if we've reached the parallel limit
            if [[ $running -ge $parallel_jobs ]]; then
                wait "${pids[0]}"
                if [[ $? -eq 0 ]]; then
                    COMPLETED=$((COMPLETED + 1))
                else
                    FAILED=$((FAILED + 1))
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
            COMPLETED=$((COMPLETED + 1))
        else
            FAILED=$((FAILED + 1))
        fi
    done
}
