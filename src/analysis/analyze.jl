using Printf

export analyze_simulation, analyze_ring

"""
    analyze_ring(run_name::String, ring_number::Int, coords_history, tangents_history,
                 dt::Float64, logger_steps::Int; output_dir::String="_data/csv")

Perform all analysis for a single ring and save results to CSV.

Computes:
- Radius of gyration (Rg) and principal components (Rg1, Rg2, Rg3)
- Mean Squared Displacement (MSD)
- End-to-end distance (Rs)
- Tangent correlation function (beta)
"""
function analyze_ring(run_name::String, ring_number::Int, coords_history, tangents_history,
                      dt::Float64, logger_steps::Int; output_dir::String="_data/csv")

    println("Analyzing ring $ring_number...")

    # Compute Rg components
    println("  Computing radius of gyration...")
    Rg1, Rg2, Rg3, Rg = compute_rg_timeseries(coords_history)

    # Compute MSD
    println("  Computing MSD...")
    msd = compute_msd(coords_history)

    # Compute Rs
    println("  Computing Rs...")
    Rs = compute_rs(coords_history)

    # Compute beta
    println("  Computing beta (tangent correlation)...")
    beta_vals = compute_beta(tangents_history)

    # Save all results
    println("  Saving results...")
    save_rg_csv(run_name, ring_number, dt, logger_steps, Rg, Rg1, Rg2, Rg3; output_dir=output_dir)
    save_msd_csv(run_name, ring_number, dt, logger_steps, msd; output_dir=output_dir)
    save_rs_csv(run_name, ring_number, Rs; output_dir=output_dir)
    save_beta_csv(run_name, ring_number, beta_vals; output_dir=output_dir)

    println("  ✓ Ring $ring_number analysis complete")

    return (Rg=Rg, Rg1=Rg1, Rg2=Rg2, Rg3=Rg3, msd=msd, Rs=Rs, beta=beta_vals)
end

"""
    analyze_simulation(filepath::String; phase::Symbol=:active, run_name::String="",
                      output_dir::String="_data/csv")

Analyze a simulation from JLD2 file and save all results to CSV.

Automatically handles both single and double ring systems.
"""
function analyze_simulation(filepath::String; phase::Symbol=:active, run_name::String="",
                           output_dir::String="_data/csv")

    println("Loading simulation data from: $filepath")

    # Load data
    coords, tangents, params = load_simulation_data(filepath; phase=phase)

    # Generate run name if not provided
    if isempty(run_name)
        run_name = generate_run_name(params)
    end

    println("Run name: $run_name")
    println("System type: $(params.system_type)")

    results = Dict()

    if params.system_type == :single
        # Single ring system
        println("\nAnalyzing single ring system...")

        results[:ring1] = analyze_ring(run_name, 1, coords, tangents,
                                       params.dt, params.logger_steps;
                                       output_dir=output_dir)

    else  # :double
        # Double ring system
        println("\nAnalyzing double ring system...")

        n1 = params.n_monomers_1
        n2 = params.n_monomers_2

        # Split data
        coords1, coords2 = split_rings_coords(coords, n1, n2)
        tangents1, tangents2 = split_rings_tangents(tangents, n1, n2)

        # Analyze each ring
        results[:ring1] = analyze_ring(run_name, 1, coords1, tangents1,
                                       params.dt, params.logger_steps;
                                       output_dir=output_dir)

        results[:ring2] = analyze_ring(run_name, 2, coords2, tangents2,
                                       params.dt, params.logger_steps;
                                       output_dir=output_dir)
    end

    println("\n✅ Analysis complete!")
    println("Results saved to: $output_dir/")

    return results
end

"""
Generate a run name from parameters
"""
function generate_run_name(params)
    if params.system_type == :single
        return "single_$(params.n_monomers)_$(params.n_active)_$(params.kangle)_$(params.factive)"
    else
        return "double_$(params.n_monomers_1)_$(params.n_monomers_2)_$(params.n_active_1)_$(params.n_active_2)_$(params.kangle)_$(params.factive)"
    end
end
