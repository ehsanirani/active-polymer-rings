export DoubleRingBuilder, create_double_ring_system, create_double_ring_active_system

function DoubleRingBuilder(params::Parameters, sim_bodies::SimBodies; loggers)
    if params.kangle > 0.0
        return create_double_ring_hard_system(params, sim_bodies; loggers=loggers)
    else
        return create_double_ring_active_system(params, sim_bodies; loggers=loggers)
    end
end

function create_double_ring_hard_system(params::Parameters, sim_bodies::SimBodies; loggers)
    n1 = params.n_monomers_1
    n2 = params.n_monomers_2
    n_total = n1 + n2
    
    polymer_bonds = InteractionList2Atoms(
        collect(1:n_total),
        append!(push!(collect(2:n1), 1),
                push!(collect(n1+2:n_total), n1+1)),
        [FENEBond(r0=1.6, k=params.kbond, σ=1.0, ϵ=1.0) for _ in 1:n_total]
    )
    
    polymer_angles = InteractionList3Atoms(
        collect(1:n_total),
        append!(push!(collect(2:n1), 1),
                push!(collect(n1+2:n_total), n1+1)),
        append!(push!(collect(3:n1), 1, 2),
                push!(collect(n1+3:n_total), n1+1, n1+2)),
        [CosineAngle(; k=params.kangle, θ0=0) for _ in 1:n_total],
        repeat(["cangle"], n_total)
    )
    
    pair_inters = (LennardJones(
        cutoff=ShiftedPotentialCutoff(2.0^(1/6)),
        use_neighbors=true
    ),)
    
    specific_inters = (polymer_bonds, polymer_angles)

    # Include ActiveTangentForce for activity (k=0 since bending handled by CosineAngle)
    general_inters = (
        ActiveTangentForce(
            system_type=:double,
            n_monomers_1=n1,
            n_monomers_2=n2,
            k=0.0,  # Bending already handled by CosineAngle
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

function create_double_ring_active_system(params::Parameters, sim_bodies::SimBodies; loggers)
    n1 = params.n_monomers_1
    n2 = params.n_monomers_2
    n_total = n1 + n2
    
    polymer_bonds = InteractionList2Atoms(
        collect(1:n_total),
        append!(push!(collect(2:n1), 1),
                push!(collect(n1+2:n_total), n1+1)),
        [FENEBond(r0=1.6, k=params.kbond, σ=1.0, ϵ=1.0) for _ in 1:n_total]
    )
    
    pair_inters = (LennardJones(
        cutoff=ShiftedPotentialCutoff(2.0^(1/6)),
        use_neighbors=true
    ),)
    
    specific_inters = (polymer_bonds,)
    
    general_inters = (
        ActiveTangentForce(
            system_type=:double,
            n_monomers_1=n1,
            n_monomers_2=n2,
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