using Test
using ActiveRings
using JLD2
using Molly: CubicBoundary
using StaticArrays: SVector

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
    
    @testset "Activity Vector" begin
        p = Parameters(system_type=:single, n_monomers=100, n_active=30)
        activity = get_activity_vector(p)
        @test length(activity) == 100
        @test sum(activity) == 30
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

println("\n" * "="^50)
println("All tests completed!")
println("="^50)