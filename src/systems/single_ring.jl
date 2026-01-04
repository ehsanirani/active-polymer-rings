export SingleRingBuilder, create_single_ring_system, create_single_ring_active_system

function SingleRingBuilder(params::Parameters, sim_bodies::SimBodies; loggers)
    if params.kangle > 0.0
        return create_single_ring_hard_system(params, sim_bodies; loggers=loggers)
    else
        return create_single_ring_active_system(params, sim_bodies; loggers=loggers)
    end
end

function create_single_ring_hard_system(params::Parameters, sim_bodies::SimBodies; loggers)
    n_particles = params.n_monomers
    n_monomers = params.n_monomers

    polymer_bonds = InteractionList2Atoms(
        collect(1:n_monomers),
        push!(collect(2:n_monomers), 1),
        [FENEBond(r0=1.6, k=params.kbond, σ=1.0, ϵ=1.0) for _ in 1:n_particles]
    )
    
    polymer_angles = InteractionList3Atoms(
        collect(1:n_monomers),
        push!(collect(2:n_monomers), 1),
        push!(collect(3:n_monomers), 1, 2),
        repeat(["cangle"], n_particles),
        [CosineAngle(; k=params.kangle, θ0=0) for _ in 1:n_particles]
    )
    
    pair_inters = (LennardJones(
        cutoff=ShiftedPotentialCutoff(2.0^(1/6)),
        use_neighbors=true
    ),)
    
    specific_inters = (polymer_bonds, polymer_angles)
    
    general_inters = (
        LangevinThermostat(
            KT=params.KT, γ=params.γ, δt=params.dt,
            frand=zero(sim_bodies.coords),
            fdamp=zero(sim_bodies.coords)
        ),
    )
    
    return Molly.System(
        atoms=sim_bodies.atoms,
        coords=sim_bodies.coords,
        velocities=sim_bodies.velocities,
        boundary=sim_bodies.boundary,
        pairwise_inters=pair_inters,
        specific_inter_lists=specific_inters,
        general_inters=general_inters,
        neighbor_finder=sim_bodies.neighbor_finder,
        loggers=loggers,
        force_units=Molly.NoUnits,
        energy_units=Molly.NoUnits
    )
end

function create_single_ring_active_system(params::Parameters, sim_bodies::SimBodies; loggers)
    n_particles = params.n_monomers
    n_monomers = params.n_monomers
    
    polymer_bonds = InteractionList2Atoms(
        collect(1:n_monomers),
        push!(collect(2:n_monomers), 1),
        [FENEBond(r0=1.6, k=params.kbond, σ=1.0, ϵ=1.0) for _ in 1:n_particles]
    )

    pair_inters = (LennardJones(
        cutoff=ShiftedPotentialCutoff(2.0^(1/6)),
        use_neighbors=true
    ),)
    
    specific_inters = (polymer_bonds,)
    
    general_inters = (
        ActiveTangentForce(
            system_type=:single,
            n_monomers_1=n_monomers,
            k=params.kangle,
            f_active=params.factive
        ),
        LangevinThermostat(
            KT=params.KT, γ=params.γ, δt=params.dt,
            frand=zero(sim_bodies.coords),
            fdamp=zero(sim_bodies.coords)
        ),
    )
    
    return Molly.System(
        atoms=sim_bodies.atoms,
        coords=sim_bodies.coords,
        velocities=sim_bodies.velocities,
        boundary=sim_bodies.boundary,
        pairwise_inters=pair_inters,
        specific_inter_lists=specific_inters,
        general_inters=general_inters,
        neighbor_finder=sim_bodies.neighbor_finder,
        loggers=loggers,
        force_units=Molly.NoUnits,
        energy_units=Molly.NoUnits
    )
end