export RgLogger, CoordinateLogger, TangentLogger, MSDLogger

struct RgLogger
    n_steps::Int64
    system_type::Symbol
    n_monomers_1::Int64
    n_monomers_2::Int64
    Rg_total::Vector{Float64}
    Rg1::Vector{Float64}
    Rg2::Vector{Float64}
    Rg3::Vector{Float64}
end

function RgLogger(n_steps::Int64, system_type::Symbol, n_monomers_1::Int64, n_monomers_2::Int64=0)
    RgLogger(n_steps, system_type, n_monomers_1, n_monomers_2,
             Float64[], Float64[], Float64[], Float64[])
end

function Molly.log_property!(logger::RgLogger, sys, buffers, neighbors=nothing, step_n::Integer=0; n_threads=Threads.nthreads(), kwargs...)
    if step_n % logger.n_steps != 0
        return
    end
    
    if logger.system_type == :single
        unwrap_coords = unwrap_polymer(sys.coords[1:logger.n_monomers_1], sys.boundary)
        Rg = radius_gyration(unwrap_coords[1:logger.n_monomers_1], 
                            sys.atoms[1:logger.n_monomers_1])
        push!(logger.Rg_total, Rg)
        
        # For single ring, all three values are the same
        push!(logger.Rg1, Rg)
        push!(logger.Rg2, Rg)
        push!(logger.Rg3, Rg)
    else
        n_total = logger.n_monomers_1 + logger.n_monomers_2
        
        unwrap_coords1 = unwrap_polymer(sys.coords[1:logger.n_monomers_1], sys.boundary)
        unwrap_coords2 = unwrap_polymer(
            sys.coords[logger.n_monomers_1+1:n_total], sys.boundary)
        
        # Radius of gyration for each ring
        Rg1 = radius_gyration(unwrap_coords1, sys.atoms[1:logger.n_monomers_1])
        Rg2 = radius_gyration(unwrap_coords2,
                             sys.atoms[logger.n_monomers_1+1:n_total])
        
        # Combined Rg for both rings
        all_coords = vcat(unwrap_coords1, unwrap_coords2)
        all_atoms = sys.atoms[1:n_total]
        Rg_total = radius_gyration(all_coords, all_atoms)
        
        push!(logger.Rg_total, Rg_total)
        push!(logger.Rg1, Rg1)
        push!(logger.Rg2, Rg2)
        push!(logger.Rg3, Rg2)  # Using Rg2 again as placeholder
    end
end

function unwrap_polymer(coords::Vector{SVector{3,Float64}}, boundary::CubicBoundary)
    L = boundary.side_lengths[1]
    Lhalf = L / 2
    bonds = coords[2:end] - coords[1:end-1]
    image = zero(coords)
    
    for i in 1:length(bonds)
        dr = bonds[i]
        passed_pbc = abs.(dr) .> Lhalf
        image0 = any(passed_pbc) ? Int.(sign.(dr)) .* passed_pbc : SVector{3,Int64}(0, 0, 0)
        image[i+1:end] = image[i+1:end] .+ Ref(image0)
    end
    
    unwrapped = coords .- (image * L)
    return unwrapped
end

"""
Logger for computing non-time-averaged Mean-Squared Displacement during simulation.

Stores reference coordinates at initialization and computes MSD at each logging step.
"""
struct MSDLogger
    n_steps::Int64
    system_type::Symbol
    n_monomers_1::Int64
    n_monomers_2::Int64
    reference_coords::Vector{SVector{3,Float64}}  # Reference coordinates (t=0)
    reference_com::SVector{3,Float64}              # Reference center of mass
    msd_monomer::Vector{Float64}                   # Monomer MSD timeseries
    msd_com::Vector{Float64}                       # COM MSD timeseries
end

"""
    MSDLogger(n_steps, system_type, n_monomers_1, n_monomers_2, initial_coords, boundary)

Create MSDLogger with reference coordinates from the first frame.
"""
function MSDLogger(n_steps::Int64, system_type::Symbol, n_monomers_1::Int64,
                   n_monomers_2::Int64, initial_coords::Vector{SVector{3,Float64}},
                   boundary)
    # Store reference coordinates (unwrapped)
    n_total = system_type == :single ? n_monomers_1 : n_monomers_1 + n_monomers_2
    ref_coords = unwrap_polymer(initial_coords[1:n_total], boundary)

    # Compute reference COM
    ref_com = sum(ref_coords) / length(ref_coords)

    MSDLogger(n_steps, system_type, n_monomers_1, n_monomers_2,
              ref_coords, ref_com, Float64[], Float64[])
end

"""
Log MSD values at each step.
"""
function Molly.log_property!(logger::MSDLogger, sys, buffers, neighbors=nothing,
                             step_n::Integer=0; n_threads=1, kwargs...)
    if step_n % logger.n_steps != 0
        return
    end

    n_total = logger.system_type == :single ? logger.n_monomers_1 :
              logger.n_monomers_1 + logger.n_monomers_2

    # Get current coordinates (unwrapped)
    current_coords = unwrap_polymer(sys.coords[1:n_total], sys.boundary)

    # Compute monomer MSD
    diff = current_coords .- logger.reference_coords
    msd_mon = sum(d -> dot(d, d), diff) / length(diff)
    push!(logger.msd_monomer, msd_mon)

    # Compute COM MSD
    current_com = sum(current_coords) / length(current_coords)
    com_diff = current_com - logger.reference_com
    msd_com_val = dot(com_diff, com_diff)
    push!(logger.msd_com, msd_com_val)
end

struct TangentLogger
    n_steps::Int64
    system_type::Symbol
    n_monomers_1::Int64
    n_monomers_2::Int64
    tang_vecs::Vector{Vector{SVector{3,Float64}}}
    active::Vector{Vector{Bool}}
end

function TangentLogger(n_steps::Int64, params::Parameters)
    system_type = params.system_type
    n_monomers_1 = params.system_type == :single ? params.n_monomers : params.n_monomers_1
    n_monomers_2 = params.system_type == :single ? 0 : params.n_monomers_2

    TangentLogger(n_steps, system_type, n_monomers_1, n_monomers_2, [], [])
end

function Molly.log_property!(logger::TangentLogger, sys, buffers, neighbors=nothing, step_n::Integer=0; n_threads=Threads.nthreads(), kwargs...)
    if step_n % logger.n_steps == 0
        n_total = logger.n_monomers_1 + logger.n_monomers_2
        
        tang_vecs0 = [sys.atoms[id].tang_vec for id in 1:n_total]
        active0 = [sys.atoms[id].active for id in 1:length(sys.atoms)]
        
        push!(logger.tang_vecs, tang_vecs0)
        push!(logger.active, active0)
    end
end