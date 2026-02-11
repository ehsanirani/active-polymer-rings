#!/usr/bin/env bash
# =============================================================================
# sweep.sh - Unified parameter sweep tool
# =============================================================================
#
# Subcommands:
#   param - Single parameter sweep (or replicates)
#   grid  - Multi-parameter grid sweep
#
# Features:
#   - General base state support (passive/active equilibration)
#   - Per-run base states (independent equilibrated states per replica)
#   - Range notation (e.g., "0:25:200")
#   - Parallel execution
#   - Dry-run mode
#
# Usage:
#   ./tools/sweep.sh param <config> [options] --sweep --<param> "values"
#   ./tools/sweep.sh grid <config> [options] --sweep --<p1> "v1" --sweep --<p2> "v2"
#
# Examples:
#   # Simple sweep (no base state)
#   ./tools/sweep.sh param config/single_ring.toml --sweep --fact "1 2 3"
#
#   # Sweep with passive base state
#   ./tools/sweep.sh param config/single_ring.toml \
#       --base-state passive --n-monomers 200 \
#       --sweep --fact "1 2 3 4 5" --runs 5 --parallel 4
#
#   # Sweep n-active with range notation
#   ./tools/sweep.sh param config/single_ring.toml \
#       --base-state active --sweep --n-active "0:25:200" --runs 3
#
#   # Grid sweep with passive base
#   ./tools/sweep.sh grid config/single_ring.toml \
#       --base-state passive \
#       --sweep --n-active "50 100" --sweep --fact "3 5" --parallel 8
#
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

# Source the shared library
source "$SCRIPT_DIR/lib/sweep_lib.sh"

# =============================================================================
# Global variables
# =============================================================================

SUBCOMMAND=""
CONFIG_FILE=""
FIXED_OVERRIDES=""

# Sweep parameters
declare -a PARAM_NAMES=()
declare -a PARAM_VALUES=()

# Options
PARALLEL_JOBS=1
DRY_RUN=false
SKIP_EXISTING=false
PREFIX="run"
NUM_RUNS=1
RUN_START=0
SWEEP_ID=""

# Base state options
BASE_STATE_TYPE=""
BASE_STATE_PER_RUN=false
N_MONOMERS=""
THERMAL_STEPS=2000000
FACT=5.0

# Data directory
DATA_DIR="_data"

# =============================================================================
# Usage functions
# =============================================================================

usage() {
    echo "Usage: $0 <subcommand> [options]"
    echo ""
    echo "Subcommands:"
    echo "  param   Single parameter sweep (or replicates)"
    echo "  grid    Multi-parameter grid sweep"
    echo ""
    echo "Run '$0 <subcommand> --help' for subcommand-specific help."
    exit 1
}

usage_param() {
    echo "Usage: $0 param [options] --sweep --<param> \"values\""
    echo ""
    echo "Single parameter sweep or replicates."
    echo ""
    echo "Options:"
    echo "  --config FILE              Config file (optional)"
    echo "  --data-dir DIR             Base data directory (default: _data)"
    echo "  --sweep --<param> \"vals\"  Parameter to sweep (supports range notation)"
    echo "  --runs N                   Number of replicates per value (default: 1)"
    echo "  --run-start N              Starting run index (default: 0)"
    echo "  --parallel N               Run N simulations in parallel (default: 1)"
    echo "  --prefix STR               Prefix for simulation IDs (default: run)"
    echo "  --sweep-id ID              Identifier for this sweep (used in metadata filename)"
    echo "  --dry-run                  Print commands without executing"
    echo "  --skip-existing            Skip simulations whose output files already exist"
    echo "                             (checks by simid only - best for resuming same sweep)"
    echo ""
    echo "Base state options:"
    echo "  --base-state TYPE          passive, active, or file path"
    echo "  --base-state-per-run       Create separate base state for each replica"
    echo "  --n-monomers N             Number of monomers (required with passive/active)"
    echo "  --thermal-steps N          Thermal steps for base state (default: 2000000)"
    echo "  --fact F                   Active force for base state (default: 5.0)"
    echo ""
    echo "Parameter overrides:"
    echo "  Any --<param> value not after --sweep is passed as a fixed override"
    echo "  Example: --n-monomers 150 --kangle 5.0"
    echo ""
    echo "Range notation:"
    echo "  \"0:25:200\"               Expands to \"0 25 50 75 100 125 150 175 200\""
    echo "  \"1.0:0.5:3.0\"             Expands to \"1.0 1.5 2.0 2.5 3.0\""
    echo ""
    echo "Examples:"
    echo "  # Simple sweep"
    echo "  $0 param --config config/single_ring.toml --sweep --fact \"1.0 2.0 3.0\""
    echo ""
    echo "  # With parameter overrides"
    echo "  $0 param --config config/single_ring.toml --n-monomers 150 --kangle 5.0 \\"
    echo "      --sweep --fact \"1 2 3\""
    echo ""
    echo "  # With passive base state"
    echo "  $0 param --config config/single_ring.toml --base-state passive --n-monomers 200 \\"
    echo "      --sweep --n-active \"0:25:200\""
    echo ""
    echo "  # With per-run base states (independent equilibrated states per replica)"
    echo "  $0 param --config config/single_ring.toml --base-state passive --base-state-per-run \\"
    echo "      --n-monomers 200 --sweep --n-active \"0 50 100\" --runs 3"
    exit 1
}

usage_grid() {
    echo "Usage: $0 grid [options] --sweep --<p1> \"v1\" --sweep --<p2> \"v2\""
    echo ""
    echo "Multi-parameter grid sweep."
    echo ""
    echo "Options:"
    echo "  --config FILE              Config file (optional)"
    echo "  --data-dir DIR             Base data directory (default: _data)"
    echo "  --sweep --<param> \"vals\"  Parameter to sweep (can specify multiple)"
    echo "  --runs N                   Number of replicates per combination (default: 1)"
    echo "  --parallel N               Run N simulations in parallel (default: 1)"
    echo "  --prefix STR               Prefix for simulation IDs (default: run)"
    echo "  --sweep-id ID              Identifier for this sweep (used in metadata filename)"
    echo "  --dry-run                  Print commands without executing"
    echo "  --skip-existing            Skip simulations whose output files already exist"
    echo "                             (checks by simid only - best for resuming same sweep)"
    echo ""
    echo "Base state options:"
    echo "  --base-state TYPE          passive, active, or file path"
    echo "  --base-state-per-run       Create separate base state for each replica"
    echo "  --n-monomers N             Number of monomers (required with passive/active)"
    echo "  --thermal-steps N          Thermal steps for base state (default: 2000000)"
    echo "  --fact F                   Active force for base state (default: 5.0)"
    echo ""
    echo "Parameter overrides:"
    echo "  Any --<param> value not after --sweep is passed as a fixed override"
    echo ""
    echo "Examples:"
    echo "  # Grid sweep (3x3 = 9 simulations)"
    echo "  $0 grid --config config/single_ring.toml \\"
    echo "      --sweep --n-active \"10 30 50\" --sweep --fact \"1.0 3.0 5.0\""
    echo ""
    echo "  # With base state (2x3x5 = 30 simulations)"
    echo "  $0 grid --config config/single_ring.toml \\"
    echo "      --base-state passive --n-monomers 200 \\"
    echo "      --sweep --n-active \"50 100\" --sweep --fact \"1.0 2.0 3.0\" --runs 5"
    echo ""
    echo "  # With per-run base states (independent equilibrated states per replica)"
    echo "  $0 grid --config config/single_ring.toml --base-state passive --base-state-per-run \\"
    echo "      --n-monomers 200 --sweep --n-active \"50 100\" --sweep --fact \"1.0 3.0\" --runs 3"
    exit 1
}

# =============================================================================
# Argument parsing
# =============================================================================

parse_common_args() {
    local in_sweep=false

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --config)
                CONFIG_FILE="$2"
                shift 2
                ;;
            --data-dir)
                DATA_DIR="$2"
                shift 2
                ;;
            --parallel)
                PARALLEL_JOBS="$2"
                shift 2
                ;;
            --dry-run)
                DRY_RUN=true
                shift
                ;;
            --skip-existing)
                SKIP_EXISTING=true
                shift
                ;;
            --prefix)
                PREFIX="$2"
                shift 2
                ;;
            --runs)
                NUM_RUNS="$2"
                shift 2
                ;;
            --run-start)
                RUN_START="$2"
                shift 2
                ;;
            --sweep-id)
                SWEEP_ID="$2"
                shift 2
                ;;
            --base-state)
                BASE_STATE_TYPE="$2"
                shift 2
                ;;
            --base-state-per-run)
                BASE_STATE_PER_RUN=true
                shift
                ;;
            --n-monomers)
                N_MONOMERS="$2"
                # Also add to fixed overrides
                FIXED_OVERRIDES="$FIXED_OVERRIDES --n-monomers $2"
                shift 2
                ;;
            --thermal-steps)
                THERMAL_STEPS="$2"
                shift 2
                ;;
            --fact)
                if [[ "$in_sweep" == true ]]; then
                    PARAM_NAMES+=("--fact")
                    PARAM_VALUES+=("$2")
                    in_sweep=false
                else
                    FACT="$2"
                    FIXED_OVERRIDES="$FIXED_OVERRIDES --fact $2"
                fi
                shift 2
                ;;
            --help|-h)
                case "$SUBCOMMAND" in
                    param) usage_param ;;
                    grid) usage_grid ;;
                    *) usage ;;
                esac
                ;;
            --sweep)
                in_sweep=true
                shift
                ;;
            --*)
                if [[ "$in_sweep" == true ]]; then
                    # This is a sweep parameter
                    PARAM_NAMES+=("$1")
                    PARAM_VALUES+=("$2")
                    in_sweep=false
                    shift 2
                else
                    # This is a fixed override
                    FIXED_OVERRIDES="$FIXED_OVERRIDES $1 $2"
                    shift 2
                fi
                ;;
            *)
                log_error "Unknown argument: $1"
                case "$SUBCOMMAND" in
                    param) usage_param ;;
                    grid) usage_grid ;;
                    *) usage ;;
                esac
                ;;
        esac
    done
}

parse_args() {
    if [[ $# -lt 1 ]]; then
        usage
    fi

    SUBCOMMAND="$1"
    shift

    case "$SUBCOMMAND" in
        param|grid)
            parse_common_args "$@"

            # Validate config file if specified
            if [[ -n "$CONFIG_FILE" && ! -f "$CONFIG_FILE" ]]; then
                log_error "Config file not found: $CONFIG_FILE"
                exit 1
            fi
            ;;
        --help|-h)
            usage
            ;;
        *)
            log_error "Unknown subcommand: $SUBCOMMAND"
            usage
            ;;
    esac
}

# =============================================================================
# Validation
# =============================================================================

validate_param_args() {
    # For param subcommand, require either --sweep or --runs >= 2
    if [[ ${#PARAM_NAMES[@]} -eq 0 && "$NUM_RUNS" -lt 2 ]]; then
        log_error "Must specify --sweep or --runs N (with N >= 2)"
        usage_param
    fi

    # Check that at most one sweep parameter is specified
    if [[ ${#PARAM_NAMES[@]} -gt 1 ]]; then
        log_error "Only one parameter can be swept with 'param' subcommand"
        echo "Use 'grid' subcommand for multiple parameters"
        exit 1
    fi

    # If base state is passive or active, require n_monomers
    if [[ "$BASE_STATE_TYPE" == "passive" || "$BASE_STATE_TYPE" == "active" ]]; then
        if [[ -z "$N_MONOMERS" ]]; then
            log_error "--n-monomers is required when using --base-state passive/active"
            exit 1
        fi
    fi
}

validate_grid_args() {
    # For grid subcommand, require at least one sweep parameter
    if [[ ${#PARAM_NAMES[@]} -eq 0 ]]; then
        log_error "At least one sweep parameter is required (use --sweep --<param> \"values\")"
        usage_grid
    fi

    # If base state is passive or active, require n_monomers
    if [[ "$BASE_STATE_TYPE" == "passive" || "$BASE_STATE_TYPE" == "active" ]]; then
        if [[ -z "$N_MONOMERS" ]]; then
            log_error "--n-monomers is required when using --base-state passive/active"
            exit 1
        fi
    fi
}

# =============================================================================
# Generate simid from parameters
# =============================================================================

# Generate simid for distinguishing replicates
# Arguments:
#   $1 - run_index: replicate index
# Always includes run index for consistency and resumability
generate_simid() {
    local run_index="$1"
    echo "${PREFIX}${run_index}"
}

# =============================================================================
# Build parameter description
# =============================================================================

build_param_desc() {
    local combo_str="$1"

    IFS=',' read -ra values <<< "$combo_str"

    local desc=""
    for idx in "${!PARAM_NAMES[@]}"; do
        if [[ -n "$desc" ]]; then
            desc="$desc, "
        fi
        desc="$desc${PARAM_NAMES[$idx]}=${values[$idx]}"
    done

    echo "$desc"
}

# =============================================================================
# Main command execution
# =============================================================================

cmd_param() {
    validate_param_args

    cd "$PROJECT_DIR"

    # Expand values if sweep parameter specified
    local has_sweep=false
    local expanded_values=""
    local n_values=1

    if [[ ${#PARAM_NAMES[@]} -gt 0 ]]; then
        has_sweep=true
        expanded_values="$(expand_values "${PARAM_VALUES[0]}")"
        read -ra VALUES_ARRAY <<< "$expanded_values"
        n_values=${#VALUES_ARRAY[@]}
    fi

    local total=$((n_values * NUM_RUNS))

    # Print header
    print_header "Parameter Sweep"
    if [[ -n "$CONFIG_FILE" ]]; then
        echo -e "Config:    $CONFIG_FILE"
    fi
    echo -e "Data dir:  $DATA_DIR"
    if [[ -n "$FIXED_OVERRIDES" ]]; then
        echo -e "Overrides:$FIXED_OVERRIDES"
    fi
    if [[ "$has_sweep" == true ]]; then
        echo -e "Sweep:     ${PARAM_NAMES[0]} = ${VALUES_ARRAY[*]}"
    fi
    if [[ -n "$BASE_STATE_TYPE" ]]; then
        if [[ "$BASE_STATE_PER_RUN" == true ]]; then
            echo -e "Base:      $BASE_STATE_TYPE (per-run)"
        else
            echo -e "Base:      $BASE_STATE_TYPE"
        fi
    fi
    if [[ "$NUM_RUNS" -gt 1 ]]; then
        echo -e "Runs:      $NUM_RUNS replicates"
    fi
    if [[ "$RUN_START" -gt 0 ]]; then
        echo -e "Run range: run${RUN_START} to run$((RUN_START + NUM_RUNS - 1))"
    fi
    echo -e "Total:     $total simulations"
    echo -e "Parallel:  $PARALLEL_JOBS jobs"
    echo -e "Prefix:    $PREFIX"
    echo ""

    if [[ "$DRY_RUN" == true ]]; then
        echo -e "${YELLOW}Dry run mode - commands will be printed but not executed${NC}"
        echo ""
    fi

    # Save sweep metadata
    save_sweep_metadata "$DATA_DIR" "$SWEEP_ID" "param" "$CONFIG_FILE" "$FIXED_OVERRIDES" \
        "$DRY_RUN" "$N_MONOMERS" "$BASE_STATE_TYPE" "$BASE_STATE_PER_RUN" "$THERMAL_STEPS" \
        "$FACT" "$NUM_RUNS" "$RUN_START" "$PREFIX" "$PARALLEL_JOBS" "$total"

    # Handle base state if specified
    local base_state_path=""
    if [[ -n "$BASE_STATE_TYPE" ]]; then
        if [[ "$BASE_STATE_PER_RUN" == true ]]; then
            # Create one base state per run
            if ! ensure_base_states_per_run "$BASE_STATE_TYPE" "$CONFIG_FILE" "$N_MONOMERS" \
                "$THERMAL_STEPS" "$FACT" "$DRY_RUN" "$DATA_DIR" "$RUN_START" "$NUM_RUNS"; then
                exit 1
            fi
            # BASE_STATE_PATHS array now populated
        else
            # Shared base state for all runs
            if ! ensure_base_state "$BASE_STATE_TYPE" "$CONFIG_FILE" "$N_MONOMERS" "$THERMAL_STEPS" "$FACT" "$DRY_RUN" "$DATA_DIR"; then
                exit 1
            fi
            base_state_path="$BASE_STATE_PATH"
        fi
        echo ""
    fi

    # Build commands
    declare -a PENDING_COMMANDS=()
    declare -a PENDING_DESCRIPTIONS=()
    local SKIPPED=0

    if [[ "$DRY_RUN" != true ]]; then
        print_header "Running Sweep Simulations"
        echo ""
    fi

    if [[ "$has_sweep" == false ]]; then
        # Just replicates
        for run in $(seq "$RUN_START" $((RUN_START + NUM_RUNS - 1))); do
            local simid="${PREFIX}_${run}"
            # Check if output already exists
            if [[ "$SKIP_EXISTING" == true ]] && simulation_output_exists "$simid" "$DATA_DIR"; then
                log_skip "simid=$simid (output exists)"
                SKIPPED=$((SKIPPED + 1))
                continue
            fi
            # Select appropriate base state
            if [[ "$BASE_STATE_PER_RUN" == true && -n "$BASE_STATE_TYPE" ]]; then
                local run_idx=$((run - RUN_START))
                base_state_path="${BASE_STATE_PATHS[$run_idx]}"
            fi
            local cmd
            cmd="$(build_simulation_command "$CONFIG_FILE" "$FIXED_OVERRIDES" "" "$simid" "$base_state_path" "$DATA_DIR")"
            PENDING_COMMANDS+=("$cmd")
            PENDING_DESCRIPTIONS+=("simid=$simid")
        done
    else
        # Sweep with optional replicates
        for i in "${!VALUES_ARRAY[@]}"; do
            local value="${VALUES_ARRAY[$i]}"
            local param_args="${PARAM_NAMES[0]} $value"

            for run in $(seq "$RUN_START" $((RUN_START + NUM_RUNS - 1))); do
                local simid
                simid="$(generate_simid "$run")"
                local desc="${PARAM_NAMES[0]}=$value, simid=$simid"
                # Check if output already exists
                if [[ "$SKIP_EXISTING" == true ]] && simulation_output_exists "$simid" "$DATA_DIR"; then
                    log_skip "$desc (output exists)"
                    SKIPPED=$((SKIPPED + 1))
                    continue
                fi
                # Select appropriate base state
                if [[ "$BASE_STATE_PER_RUN" == true && -n "$BASE_STATE_TYPE" ]]; then
                    local run_idx=$((run - RUN_START))
                    base_state_path="${BASE_STATE_PATHS[$run_idx]}"
                fi
                local cmd
                cmd="$(build_simulation_command "$CONFIG_FILE" "$FIXED_OVERRIDES" "$param_args" "$simid" "$base_state_path" "$DATA_DIR")"

                PENDING_COMMANDS+=("$cmd")
                PENDING_DESCRIPTIONS+=("$desc")
            done
        done
    fi

    # Run simulations
    run_parallel "$PARALLEL_JOBS" "$DRY_RUN"

    # Print summary
    if [[ "$DRY_RUN" != true ]]; then
        print_summary "$COMPLETED" "$FAILED" "$SKIPPED"
    fi
}

cmd_grid() {
    validate_grid_args

    cd "$PROJECT_DIR"

    # Generate all combinations
    generate_cartesian_product
    local n_combos=${#COMBINATIONS[@]}
    local total=$((n_combos * NUM_RUNS))

    # Print header
    print_header "Grid Parameter Sweep"
    if [[ -n "$CONFIG_FILE" ]]; then
        echo -e "Config:    $CONFIG_FILE"
    fi
    echo -e "Data dir:  $DATA_DIR"
    if [[ -n "$FIXED_OVERRIDES" ]]; then
        echo -e "Overrides:$FIXED_OVERRIDES"
    fi
    echo -e "Sweep parameters:"
    for i in "${!PARAM_NAMES[@]}"; do
        local expanded
        expanded="$(expand_values "${PARAM_VALUES[$i]}")"
        echo -e "  ${PARAM_NAMES[$i]}: $expanded"
    done
    echo -e "Combinations: $n_combos"
    if [[ -n "$BASE_STATE_TYPE" ]]; then
        if [[ "$BASE_STATE_PER_RUN" == true ]]; then
            echo -e "Base:      $BASE_STATE_TYPE (per-run)"
        else
            echo -e "Base:      $BASE_STATE_TYPE"
        fi
    fi
    if [[ "$NUM_RUNS" -gt 1 ]]; then
        echo -e "Runs:      $NUM_RUNS replicates per combination"
    fi
    echo -e "Total:     $total simulations"
    echo -e "Parallel:  $PARALLEL_JOBS jobs"
    echo -e "Prefix:    $PREFIX"
    echo ""

    if [[ "$DRY_RUN" == true ]]; then
        echo -e "${YELLOW}Dry run mode - commands will be printed but not executed${NC}"
        echo ""
    fi

    # Save sweep metadata
    save_sweep_metadata "$DATA_DIR" "$SWEEP_ID" "grid" "$CONFIG_FILE" "$FIXED_OVERRIDES" \
        "$DRY_RUN" "$N_MONOMERS" "$BASE_STATE_TYPE" "$BASE_STATE_PER_RUN" "$THERMAL_STEPS" \
        "$FACT" "$NUM_RUNS" 0 "$PREFIX" "$PARALLEL_JOBS" "$total"

    # Handle base state if specified
    local base_state_path=""
    if [[ -n "$BASE_STATE_TYPE" ]]; then
        if [[ "$BASE_STATE_PER_RUN" == true ]]; then
            # Create one base state per run
            if ! ensure_base_states_per_run "$BASE_STATE_TYPE" "$CONFIG_FILE" "$N_MONOMERS" \
                "$THERMAL_STEPS" "$FACT" "$DRY_RUN" "$DATA_DIR" 0 "$NUM_RUNS"; then
                exit 1
            fi
            # BASE_STATE_PATHS array now populated
        else
            # Shared base state for all runs
            if ! ensure_base_state "$BASE_STATE_TYPE" "$CONFIG_FILE" "$N_MONOMERS" "$THERMAL_STEPS" "$FACT" "$DRY_RUN" "$DATA_DIR"; then
                exit 1
            fi
            base_state_path="$BASE_STATE_PATH"
        fi
        echo ""
    fi

    # Build commands
    declare -a PENDING_COMMANDS=()
    declare -a PENDING_DESCRIPTIONS=()
    local SKIPPED=0

    if [[ "$DRY_RUN" != true ]]; then
        print_header "Running Sweep Simulations"
        echo ""
    fi

    for i in "${!COMBINATIONS[@]}"; do
        local combo="${COMBINATIONS[$i]}"
        IFS=',' read -ra values <<< "$combo"

        # Build parameter args
        local param_args=""
        for idx in "${!PARAM_NAMES[@]}"; do
            param_args="$param_args ${PARAM_NAMES[$idx]} ${values[$idx]}"
        done

        for run in $(seq 0 $((NUM_RUNS - 1))); do
            local simid
            simid="$(generate_simid "$run")"
            local desc
            desc="$(build_param_desc "$combo"), simid=$simid"
            # Check if output already exists
            if [[ "$SKIP_EXISTING" == true ]] && simulation_output_exists "$simid" "$DATA_DIR"; then
                log_skip "$desc (output exists)"
                SKIPPED=$((SKIPPED + 1))
                continue
            fi
            # Select appropriate base state
            if [[ "$BASE_STATE_PER_RUN" == true && -n "$BASE_STATE_TYPE" ]]; then
                base_state_path="${BASE_STATE_PATHS[$run]}"
            fi
            local cmd
            cmd="$(build_simulation_command "$CONFIG_FILE" "$FIXED_OVERRIDES" "$param_args" "$simid" "$base_state_path" "$DATA_DIR")"

            PENDING_COMMANDS+=("$cmd")
            PENDING_DESCRIPTIONS+=("$desc")
        done
    done

    # Run simulations
    run_parallel "$PARALLEL_JOBS" "$DRY_RUN"

    # Print summary
    if [[ "$DRY_RUN" != true ]]; then
        print_summary "$COMPLETED" "$FAILED" "$SKIPPED"
    fi
}

# =============================================================================
# Main
# =============================================================================

main() {
    parse_args "$@"

    case "$SUBCOMMAND" in
        param)
            cmd_param
            ;;
        grid)
            cmd_grid
            ;;
    esac
}

main "$@"
