export LangevinThermostat

mutable struct LangevinThermostat
    KT::Float64
    γ::Float64
    δt::Float64
    sqfr2::Float64
    frand::Vector{SVector{3,Float64}}
    fdamp::Vector{SVector{3,Float64}}
end

function LangevinThermostat(; frand, fdamp, KT::Float64=1.0, γ::Float64=1.0, δt::Float64=0.01)
    sqfr2 = sqrt(24 * KT * γ / δt)
    return LangevinThermostat(KT, γ, δt, sqfr2, frand, fdamp)
end

function AtomsCalculators.forces!(fs, sys, inter::LangevinThermostat;
                                   neighbors=nothing, step_n=0, n_threads=Threads.nthreads(),
                                   buffers=nothing, needs_vir=false)
    inter.fdamp = inter.γ .* sys.velocities
    inter.frand = (rand!(inter.frand) .- Ref([0.5, 0.5, 0.5])) .* inter.sqfr2

    # Add to force accumulator
    for i in 1:length(fs)
        fs[i] += inter.frand[i] - inter.fdamp[i]
    end

    return nothing
end