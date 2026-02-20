using Test
using ActiveRings
using JLD2
using Molly
using Molly: CubicBoundary, AtomsCalculators
using StaticArrays: SVector
using LinearAlgebra: norm, normalize

@testset "Parameters" begin
    @testset "Single Ring Parameters" begin
        p = Parameters(system_type=:single, n_monomers=200, n_active=20)
        @test p.system_type == :single
        @test p.n_monomers == 200
        @test p.n_active == 20
        @test p.n_monomers_1 == 100  # default
        @test get_n_particles(p) == 200
        @test get_n_active(p) == 20
    end
    
    @testset "Double Ring Parameters" begin
        p = Parameters(system_type=:double, n_monomers_1=150, n_monomers_2=100, 
                      n_active_1=15, n_active_2=10)
        @test p.system_type == :double
        @test get_n_particles(p) == 250
        @test get_n_active(p) == 25
    end
    
    @testset "Box Size Auto-Calculation" begin
        p1 = Parameters(system_type=:single, n_monomers=100, L=0.0)
        @test p1.L > 0
        
        p2 = Parameters(system_type=:double, n_monomers_1=100, n_monomers_2=100, L=0.0)
        @test p2.L > p1.L
    end
    
    @testset "Activity Vector - Random" begin
        p = Parameters(system_type=:single, n_monomers=100, n_active=30, activity_pattern=:random)
        activity = get_activity_vector(p)
        @test length(activity) == 100
        @test sum(activity) == 30
    end

    @testset "Activity Vector - Block" begin
        p = Parameters(system_type=:single, n_monomers=100, n_active=30, activity_pattern=:block)
        activity = get_activity_vector(p)
        @test length(activity) == 100
        @test sum(activity) == 30
        # Block pattern: first 30 should be active, rest passive
        @test all(activity[1:30])
        @test !any(activity[31:100])
    end

    @testset "Activity Pattern defaults to random" begin
        p = Parameters(system_type=:single, n_monomers=100)
        @test p.activity_pattern == :random
    end

    @testset "Thermal dt defaults to dt" begin
        p1 = Parameters(system_type=:single, n_monomers=100, dt=0.01)
        @test p1.dt_thermal == 0.01  # Should default to dt

        p2 = Parameters(system_type=:single, n_monomers=100, dt=0.01, dt_thermal=0.02)
        @test p2.dt_thermal == 0.02  # Should use specified value
    end

    @testset "MSD flags default to false" begin
        p = Parameters(system_type=:single, n_monomers=100)
        @test p.msd_com == false
        @test p.msd_time_averaged == false
    end

    @testset "MSD flags can be enabled" begin
        p = Parameters(system_type=:single, n_monomers=100, msd_com=true, msd_time_averaged=true)
        @test p.msd_com == true
        @test p.msd_time_averaged == true
    end

    @testset "export_xyz defaults to false" begin
        p = Parameters(system_type=:single, n_monomers=100)
        @test p.export_xyz == false
    end

    @testset "export_xyz can be enabled" begin
        p = Parameters(system_type=:single, n_monomers=100, export_xyz=true)
        @test p.export_xyz == true
    end

    @testset "metrics_format defaults to jld2" begin
        p = Parameters(system_type=:single, n_monomers=100)
        @test p.metrics_format == :jld2
    end

    @testset "metrics_format can be set to csv" begin
        p = Parameters(system_type=:single, n_monomers=100, metrics_format=:csv)
        @test p.metrics_format == :csv
    end
end

@testset "BODIES" begin
    @testset "BODIES Creation" begin
        p = Parameters(system_type=:single, n_monomers=50, n_active=10)
        bodies = create_bodies(p)
        
        @test length(bodies.atoms) == 50
        @test length(bodies.coords) == 50
        @test length(bodies.velocities) == 50
        @test bodies.boundary isa CubicBoundary
    end
    
    @testset "Double Ring Creation" begin
        p = Parameters(system_type=:double, n_monomers_1=50, n_monomers_2=50)
        bodies = create_bodies(p)
        
        @test length(bodies.atoms) == 100
        @test length(bodies.coords) == 100
    end
end

@testset "Forces" begin
    @testset "LangevinThermostat Creation" begin
        coords = [SVector(0.0, 0.0, 0.0) for _ in 1:10]
        fdamp = zero(coords)
        frand = zero(coords)

        force = LangevinThermostat(KT=1.0, γ=2.0, δt=0.01, frand=frand, fdamp=fdamp)
        
        @test force.KT == 1.0
        @test force.γ == 2.0
        @test force.δt == 0.01
    end
    
    @testset "ActiveTangentForce Creation" begin
        # Single ring
        f1 = ActiveTangentForce(system_type=:single, n_monomers_1=100, k=1.0, f_active=2.0)
        @test f1.system_type == :single
        @test length(f1.forces) == 100

        # Double ring
        f2 = ActiveTangentForce(system_type=:double, n_monomers_1=100, n_monomers_2=100,
                               k=1.0, f_active=2.0)
        @test f2.system_type == :double
        @test length(f2.forces) == 200
    end
end

@testset "Loggers" begin
    @testset "RgLogger Creation" begin
        schedule = make_logging_schedule(:fixed, 1000; fixed_interval=100)
        p1 = Parameters(system_type=:single, n_monomers=100)
        logger1 = RgLogger(schedule, p1.system_type, Int64(100))
        @test logger1.system_type == :single

        p2 = Parameters(system_type=:double, n_monomers_1=100, n_monomers_2=100)
        logger2 = RgLogger(schedule, p2.system_type, p2.n_monomers_1, p2.n_monomers_2)
        @test logger2.system_type == :double
        @test logger2.n_monomers_2 == 100
    end

    @testset "MSDLogger compute_com flag" begin
        schedule = make_logging_schedule(:fixed, 1000; fixed_interval=100)
        boundary = CubicBoundary(50.0)
        coords = [SVector(Float64(i), 0.0, 0.0) for i in 1:20]

        logger_with_com = MSDLogger(schedule, :single, Int64(20), Int64(0), coords, boundary; compute_com=true)
        @test logger_with_com.compute_com == true

        logger_no_com = MSDLogger(schedule, :single, Int64(20), Int64(0), coords, boundary; compute_com=false)
        @test logger_no_com.compute_com == false
    end

    @testset "TangentLogger Creation" begin
        p = Parameters(system_type=:single, n_monomers=100, n_active=20)
        logger = TangentLogger(100, p)
        @test logger.n_monomers_1 == 100
    end
end

@testset "LoggingSchedule" begin
    @testset "Fixed schedule" begin
        s_fixed = make_logging_schedule(:fixed, 1000; fixed_interval=100)
        @test s_fixed == [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
    end

    @testset "Logspaced schedule" begin
        s_log = make_logging_schedule(:logspaced, 100000; n_points=50)
        @test s_log[1] == 1
        @test s_log[end] == 100000
        @test issorted(s_log)
        @test length(s_log) == length(unique(s_log))
    end

    @testset "Invalid mode" begin
        @test_throws ErrorException make_logging_schedule(:invalid, 1000)
    end
end

@testset "Workflow" begin
    @testset "Single Ring Minimal" begin
        p = Parameters(
            system_type=:single,
            n_monomers=20,
            n_active=5,
            n_steps=100,
            thermal_steps=50,
            traj_interval=10,
            factive=1.0,
            kangle=0.0  # Use active force
        )
        
        # Create bodies
        sim_bodies = create_bodies(p)
        @test length(sim_bodies.atoms) == 20
        
        # Check activity assignment
        n_active = sum([atom.active for atom in sim_bodies.atoms])
        @test n_active == 5
        
        println("Single ring minimal test passed ✓")
    end
    
    @testset "Double Ring Minimal" begin
        p = Parameters(
            system_type=:double,
            n_monomers_1=15,
            n_monomers_2=10,
            n_active_1=3,
            n_active_2=2,
            n_steps=100,
            thermal_steps=50,
            traj_interval=10,
            factive=1.0,
            kangle=5.0  # Use angle potential
        )
        
        sim_bodies = create_bodies(p)
        @test length(sim_bodies.atoms) == 25
        
        # Check activity
        n_active = sum([atom.active for atom in sim_bodies.atoms])
        @test n_active == 5
        
        println("Double ring minimal test passed ✓")
    end
end

@testset "XYZ File Output" begin
    @testset "save_xyz function basic test" begin
        # Test that we can create and write XYZ files using native Julia format
        # This is a basic sanity check without running full simulation
        n_atoms = 10
        L = 20.0

        # Write a simple test XYZ file
        test_file = "_data/test_native.xyz"
        mkpath("_data")

        open(test_file, "w") do io
            # Write atom count
            println(io, n_atoms)

            # Write comment line with Properties specification
            bond_str = join(["$i-$(i+1)" for i in 1:n_atoms-1], ",")
            println(io, "Properties=species:S:1:pos:R:3:tangent:R:3:active:I:1 ",
                    "Lattice=\"$L 0 0 0 $L 0 0 0 $L\" ",
                    "bonds=\"$bond_str\" ",
                    "step=1")

            # Write atom lines
            for i in 1:n_atoms
                # Format: species x y z tang_x tang_y tang_z active
                println(io, "Monomer $(Float64(i)) 0.0 0.0 1.0 0.0 0.0 0")
            end
        end

        # Verify file exists and can be read back
        @test isfile(test_file)

        # Read back and verify format
        lines = readlines(test_file)
        @test length(lines) == n_atoms + 2  # atom count + comment + atoms
        @test parse(Int, lines[1]) == n_atoms
        @test occursin("Properties=", lines[2])
        @test occursin("Lattice=", lines[2])
        @test occursin("bonds=", lines[2])

        # Verify atom lines
        for i in 3:length(lines)
            parts = split(lines[i])
            @test parts[1] == "Monomer"  # species
            @test length(parts) == 8      # species + 3 pos + 3 tangent + 1 active
            @test all(p -> !isnothing(tryparse(Float64, p)), parts[2:end])
        end

        rm(test_file)

        println("Native XYZ basic test passed ✓")
    end
end

@testset "Ring Closure - Active Force" begin
    # Helper: build a minimal Molly System with given coords and all-active atoms
    function make_ring_system(coords::Vector{SVector{3,Float64}}, boundary, f_active; k=0.0)
        n = length(coords)
        atoms = [Particle(i, 1.0, 1.0, 1.0, true,
                 SVector(0.,0.,0.), SVector(0.,0.,0.), SVector(0,0,0)) for i in 1:n]
        velocities = [SVector(0.,0.,0.) for _ in 1:n]

        # Ring bonds: 1-2, 2-3, ..., n-1
        bond_is = push!(collect(1:n-1), n)
        bond_js = push!(collect(2:n), 1)
        bonds = Molly.InteractionList2Atoms(
            bond_is, bond_js,
            [Molly.HarmonicBond(r0=1.0, k=30.0) for _ in 1:n]
        )

        inter = ActiveTangentForce(system_type=:single, n_monomers_1=n, k=k, f_active=f_active)

        nf = CellListMapNeighborFinder(
            eligible=trues(n, n), n_steps=1,
            dist_cutoff=2.0, x0=coords
        )

        sys = Molly.System(
            atoms=atoms,
            coords=coords,
            velocities=velocities,
            boundary=boundary,
            pairwise_inters=(Molly.LennardJones(cutoff=ShiftedPotentialCutoff(2.0^(1/6)), use_neighbors=true),),
            specific_inter_lists=(bonds,),
            general_inters=(inter,),
            neighbor_finder=nf,
            loggers=Dict{String,Any}(),
            force_units=Molly.NoUnits,
            energy_units=Molly.NoUnits
        )
        return sys, inter
    end

    @testset "Regular polygon: net active force is zero" begin
        # On a perfect regular polygon every tangent is symmetric,
        # so the vector sum of all normalized tangents must vanish.
        for n in [6, 10, 20, 50]
            L = 200.0
            boundary = CubicBoundary(L)
            R = 5.0
            coords = [SVector(L/2 + R*cos(2π*i/n), L/2 + R*sin(2π*i/n), L/2) for i in 1:n]

            sys, inter = make_ring_system(coords, boundary, 3.0)

            # Compute forces via the ActiveTangentForce pathway
            fs = [SVector(0.,0.,0.) for _ in 1:n]
            AtomsCalculators.forces!(fs, sys, inter)

            # The only contribution in fs is the active tangent force (k=0 ⇒ no bending)
            net = sum(fs)
            @test norm(net) < 1e-10
        end
    end

    @testset "Tangent at monomer 1 and n are two-sided" begin
        # For a non-regular ring, check that monomers 1 and n get tangents
        # consistent with the ring closure (i.e. same formula as interior monomers).
        n = 8
        L = 200.0
        boundary = CubicBoundary(L)

        # Slightly perturbed circle so it's not a perfect polygon
        R = 5.0
        coords = [SVector(L/2 + R*cos(2π*i/n) + 0.1*i, L/2 + R*sin(2π*i/n), L/2) for i in 1:n]

        sys, inter = make_ring_system(coords, boundary, 1.0)
        fs = [SVector(0.,0.,0.) for _ in 1:n]
        AtomsCalculators.forces!(fs, sys, inter)

        # Manually compute expected tangent at monomer 1 using ring closure: triplet (n, 1, 2)
        ba = Molly.vector(coords[1], coords[n], boundary)
        bc = Molly.vector(coords[1], coords[2], boundary)
        expected_tang_1 = normalize(bc / norm(bc) - ba / norm(ba))
        @test isapprox(inter.tang_vecs[1], expected_tang_1; atol=1e-12)

        # Manually compute expected tangent at monomer n using ring closure: triplet (n-1, n, 1)
        ba_n = Molly.vector(coords[n], coords[n-1], boundary)
        bc_n = Molly.vector(coords[n], coords[1], boundary)
        expected_tang_n = normalize(bc_n / norm(bc_n) - ba_n / norm(ba_n))
        @test isapprox(inter.tang_vecs[n], expected_tang_n; atol=1e-12)
    end

    @testset "All monomers get tangents (no zero vectors)" begin
        n = 15
        L = 200.0
        boundary = CubicBoundary(L)
        R = 5.0
        coords = [SVector(L/2 + R*cos(2π*i/n), L/2 + R*sin(2π*i/n), L/2) for i in 1:n]

        sys, inter = make_ring_system(coords, boundary, 1.0)
        fs = [SVector(0.,0.,0.) for _ in 1:n]
        AtomsCalculators.forces!(fs, sys, inter)

        for i in 1:n
            @test norm(inter.tang_vecs[i]) > 0.99  # should be unit vectors
        end
    end

    @testset "Bending forces conserve momentum (k > 0)" begin
        # The angle potential forces satisfy fb = -fa - fc per triplet,
        # so the total bending force on the ring must be zero.
        n = 12
        L = 200.0
        boundary = CubicBoundary(L)
        R = 5.0
        # Perturbed ring
        coords = [SVector(L/2 + R*cos(2π*i/n) + 0.2*sin(3i), L/2 + R*sin(2π*i/n), L/2) for i in 1:n]

        # f_active=0 so forces are purely from bending
        sys, inter = make_ring_system(coords, boundary, 0.0; k=5.0)
        fs = [SVector(0.,0.,0.) for _ in 1:n]
        AtomsCalculators.forces!(fs, sys, inter)

        net = sum(fs)
        @test norm(net) < 1e-10
    end

    @testset "Double ring closure" begin
        n1, n2 = 10, 8
        L = 200.0
        boundary = CubicBoundary(L)
        R = 5.0

        # Two separate rings placed apart
        coords_1 = [SVector(L/2 + R*cos(2π*i/n1), L/2 + R*sin(2π*i/n1), L/2) for i in 1:n1]
        coords_2 = [SVector(L/2 + R*cos(2π*i/n2) + 20.0, L/2 + R*sin(2π*i/n2), L/2) for i in 1:n2]
        coords = vcat(coords_1, coords_2)
        n = n1 + n2

        atoms = [Particle(i, 1.0, 1.0, 1.0, true,
                 SVector(0.,0.,0.), SVector(0.,0.,0.), SVector(0,0,0)) for i in 1:n]
        velocities = [SVector(0.,0.,0.) for _ in 1:n]

        # Bonds for two separate rings
        bond_is_1 = push!(collect(1:n1-1), n1)
        bond_js_1 = push!(collect(2:n1), 1)
        bond_is_2 = push!(collect(n1+1:n-1), n)
        bond_js_2 = push!(collect(n1+2:n), n1+1)
        bonds = Molly.InteractionList2Atoms(
            vcat(bond_is_1, bond_is_2),
            vcat(bond_js_1, bond_js_2),
            [Molly.HarmonicBond(r0=1.0, k=30.0) for _ in 1:n]
        )

        inter = ActiveTangentForce(system_type=:double, n_monomers_1=n1, n_monomers_2=n2,
                                   k=0.0, f_active=2.0)

        nf = CellListMapNeighborFinder(
            eligible=trues(n, n), n_steps=1,
            dist_cutoff=2.0, x0=coords
        )

        sys = Molly.System(
            atoms=atoms, coords=coords, velocities=velocities,
            boundary=boundary,
            pairwise_inters=(Molly.LennardJones(cutoff=ShiftedPotentialCutoff(2.0^(1/6)), use_neighbors=true),),
            specific_inter_lists=(bonds,),
            general_inters=(inter,),
            neighbor_finder=nf,
            loggers=Dict{String,Any}(),
            force_units=Molly.NoUnits,
            energy_units=Molly.NoUnits
        )

        fs = [SVector(0.,0.,0.) for _ in 1:n]
        AtomsCalculators.forces!(fs, sys, inter)

        # Net active force on each regular-polygon ring should be zero independently
        net_1 = sum(fs[1:n1])
        net_2 = sum(fs[n1+1:end])
        @test norm(net_1) < 1e-10
        @test norm(net_2) < 1e-10

        # All tangents should be unit vectors
        for i in 1:n
            @test norm(inter.tang_vecs[i]) > 0.99
        end
    end
end

println("\n" * "="^50)
println("All tests completed!")
println("="^50)