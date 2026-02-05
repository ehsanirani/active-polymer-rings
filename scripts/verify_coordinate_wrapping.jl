#!/usr/bin/env julia

using JLD2
using LinearAlgebra
using StaticArrays

"""
Detect boundary crossing events where coordinates exceed box boundaries.
This is THE definitive test for unwrapped coordinates.
"""
function detect_boundary_crossings(coords_history, box_size)
    crossings = []

    for particle_idx in 1:length(coords_history[1])
        trajectory = [frame[particle_idx] for frame in coords_history]

        for t in 2:length(trajectory)
            prev_pos = trajectory[t-1]
            curr_pos = trajectory[t]

            # Check each dimension
            for dim in 1:3
                # Wrapped version of current position
                wrapped_curr = mod(curr_pos[dim], box_size)

                # If wrapped position differs significantly from actual position,
                # we've crossed a boundary AND coordinates are unwrapped
                if abs(curr_pos[dim] - wrapped_curr) > 0.1
                    push!(crossings, (
                        particle=particle_idx,
                        time=t,
                        dim=dim,
                        position=curr_pos[dim],
                        wrapped_position=wrapped_curr,
                        exceeded_boundary=abs(curr_pos[dim]) > box_size
                    ))
                end
            end
        end
    end

    return crossings
end

"""
Check for coordinates that directly violate box boundaries.
Only possible if coordinates are unwrapped.
"""
function check_boundary_violations(coords_history, box_size)
    violations = 0
    max_coord = 0.0
    min_coord = 0.0
    violation_examples = []

    for (frame_idx, frame) in enumerate(coords_history)
        for (particle_idx, coord) in enumerate(frame)
            for dim in 1:3
                val = coord[dim]
                max_coord = max(max_coord, val)
                min_coord = min(min_coord, val)

                if abs(val) > box_size
                    violations += 1
                    if length(violation_examples) < 5
                        push!(violation_examples, (
                            frame=frame_idx,
                            particle=particle_idx,
                            dim=dim,
                            value=val
                        ))
                    end
                end
            end
        end
    end

    return violations, max_coord, min_coord, violation_examples
end

"""
Find large discontinuous jumps in trajectories.
Wrapped coordinates show jumps ~= box_size, unwrapped don't.
"""
function find_large_jumps(coords_history, box_size)
    large_jumps = []
    threshold = box_size / 2

    for particle_idx in 1:length(coords_history[1])
        trajectory = [frame[particle_idx] for frame in coords_history]

        for t in 2:length(trajectory)
            displacement = trajectory[t] - trajectory[t-1]
            jump_size = norm(displacement)

            # A jump > box_size/2 suggests wrapping
            if jump_size > threshold
                push!(large_jumps, (
                    particle=particle_idx,
                    time=t,
                    size=jump_size,
                    displacement=displacement
                ))
            end
        end
    end

    return large_jumps
end

"""
Calculate center of mass trajectory and its drift.
Unwrapped coords allow COM to drift freely.
"""
function analyze_com_drift(coords_history, box_size)
    com_trajectory = []

    for frame in coords_history
        com = sum(frame) / length(frame)
        push!(com_trajectory, com)
    end

    initial_com = com_trajectory[1]
    final_com = com_trajectory[end]
    total_drift = norm(final_com - initial_com)
    max_drift_from_origin = maximum(norm(com) for com in com_trajectory)

    return com_trajectory, total_drift, max_drift_from_origin
end

"""
Main verification function.
"""
function verify_coordinate_wrapping(jld2_file::String)
    println("\n" * "="^70)
    println("COORDINATE WRAPPING VERIFICATION")
    println("="^70)
    println("File: $jld2_file\n")

    # Load data
    if !isfile(jld2_file)
        println("ERROR: File not found: $jld2_file")
        return
    end

    data = jldopen(jld2_file, "r")
    coords_history = data["active/loggers"]["coords"].history
    params = data["params"]
    box_size = params.L
    close(data)

    println("Box size: $box_size")
    println("Number of frames: $(length(coords_history))")
    println("Number of particles: $(length(coords_history[1]))")
    println()

    # Test 1: Boundary Violation Check
    println("TEST 1: Direct Boundary Violation Check")
    println("-"^70)
    violations, max_coord, min_coord, examples = check_boundary_violations(coords_history, box_size)

    println("Coordinate range: [$min_coord, $max_coord]")
    println("Box boundaries: [0.0, $box_size]")

    if violations > 0
        println("✓ UNWRAPPED: Found $violations coordinates exceeding box boundaries")
        println("\nExample violations:")
        for ex in examples
            println("  Frame $(ex.frame), Particle $(ex.particle), Dim $(ex.dim): $(ex.value)")
        end
    else
        println("⚠ All coordinates within box - possibly wrapped OR no crossings yet")
    end

    # Test 2: Boundary Crossing Event Detection
    println("\n\nTEST 2: Boundary Crossing Event Detection (CRITICAL)")
    println("-"^70)
    crossings = detect_boundary_crossings(coords_history, box_size)

    if length(crossings) > 0
        println("✓ FOUND $(length(crossings)) boundary crossing events")
        println("\nFirst 5 crossing events:")
        for (i, crossing) in enumerate(crossings[1:min(5, length(crossings))])
            dim_name = ["x", "y", "z"][crossing.dim]
            println("  $(i). Particle $(crossing.particle) at time $(crossing.time) ($dim_name-axis)")
            println("      Actual position: $(round(crossing.position, digits=3))")
            println("      Wrapped would be: $(round(crossing.wrapped_position, digits=3))")
            println("      Exceeded boundary: $(crossing.exceeded_boundary)")
        end
        println("  ... and $(max(0, length(crossings) - 5)) more")
    else
        println("✗ No boundary crossings detected - test inconclusive")
        println("  Need longer simulation or higher activity")
    end

    # Test 3: Trajectory Continuity (Jump Detection)
    println("\n\nTEST 3: Trajectory Continuity (Large Jump Detection)")
    println("-"^70)
    jumps = find_large_jumps(coords_history, box_size)

    if length(jumps) == 0
        println("✓ CONTINUOUS: No large discontinuous jumps (unwrapped behavior)")
        println("  All displacements < $(box_size/2) (threshold for wrapping jumps)")
    else
        println("✗ DISCONTINUOUS: Found $(length(jumps)) large jumps")
        println("  This suggests WRAPPED coordinates (unexpected!)")
        println("\nFirst 3 jumps:")
        for (i, jump) in enumerate(jumps[1:min(3, length(jumps))])
            println("  $(i). Particle $(jump.particle) at time $(jump.time)")
            println("      Jump size: $(round(jump.size, digits=3))")
            println("      Displacement: $(round.(jump.displacement, digits=3))")
        end
    end

    # Test 4: Center of Mass Drift
    println("\n\nTEST 4: Center of Mass Drift Analysis")
    println("-"^70)
    com_trajectory, total_drift, max_drift = analyze_com_drift(coords_history, box_size)

    println("Initial COM: $(round.(com_trajectory[1], digits=3))")
    println("Final COM: $(round.(com_trajectory[end], digits=3))")
    println("Total COM drift: $(round(total_drift, digits=3))")
    println("Max distance from origin: $(round(max_drift, digits=3))")

    if total_drift > box_size * 0.1
        println("✓ SIGNIFICANT DRIFT: COM moved freely (unwrapped behavior)")
    else
        println("⚠ Small drift - may need longer simulation")
    end

    # Final Conclusion
    println("\n" * "="^70)
    println("CONCLUSION:")
    println("="^70)

    unwrapped_evidence = 0
    wrapped_evidence = 0

    if violations > 0
        println("✓ Coordinates exceed box boundaries")
        unwrapped_evidence += 1
    end

    if length(crossings) > 0
        println("✓ Boundary crossing events detected")
        unwrapped_evidence += 1
    end

    if length(jumps) == 0
        println("✓ No discontinuous jumps in trajectories")
        unwrapped_evidence += 1
    else
        println("✗ Discontinuous jumps detected")
        wrapped_evidence += 1
    end

    if total_drift > box_size * 0.1
        println("✓ Significant COM drift observed")
        unwrapped_evidence += 1
    end

    println()
    if unwrapped_evidence >= 2 && wrapped_evidence == 0
        println("🎉 VERDICT: Coordinates are UNWRAPPED (continuous)")
        println("   The comment in the code is CORRECT.")
    elseif wrapped_evidence > 0
        println("⚠️  VERDICT: Evidence suggests WRAPPED coordinates!")
        println("   The comment in the code may be INCORRECT.")
    else
        println("🤷 VERDICT: Inconclusive - need longer/more active simulation")
    end

    println("="^70 * "\n")
end

# Main execution
if length(ARGS) > 0
    verify_coordinate_wrapping(ARGS[1])
else
    println("Usage: julia verify_coordinate_wrapping.jl <path_to_jld2_file>")
    println("\nExample:")
    println("  julia scripts/verify_coordinate_wrapping.jl _data/jld2/single_20_10_0.0_5.0_verify_test.jld2")
end
