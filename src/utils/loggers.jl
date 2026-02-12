export RgLogger, CoordinateLogger, TangentLogger, MSDLogger

struct RgLogger
    log_steps::Set{Int64}
    step_indices::Vector{Int64}
    system_type::Symbol
    n_monomers_1::Int64
    n_monomers_2::Int64
    Rg_total::Vector{Float64}
    Rg1::Vector{Float64}
    Rg2::Vector{Float64}
    Rg3::Vector{Float64}
end

function RgLogger(schedule::Vector{Int}, system_type::Symbol, n_monomers_1::Int64, n_monomers_2::Int64=0)
    RgLogger(Set{Int64}(schedule), Int64[], system_type, n_monomers_1, n_monomers_2,
             Float64[], Float64[], Float64[], Float64[])
end

function Molly.log_property!(logger::RgLogger, sys, buffers, neighbors=nothing, step_n::Integer=0; n_threads=Threads.nthreads(), kwargs...)
    if !(step_n in logger.log_steps)
        return
    end
    push!(logger.step_indices, step_n)
    
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
For double ring systems, computes MSD separately for each ring in addition to combined MSD.
"""
struct MSDLogger
    log_steps::Set{Int64}
    step_indices::Vector{Int64}
    compute_com::Bool                              # whether to compute COM MSD
    compute_com_frame::Bool                        # whether to compute COM-frame MSD
    system_type::Symbol
    n_monomers_1::Int64
    n_monomers_2::Int64
    # Combined reference data (all monomers)
    reference_coords::Vector{SVector{3,Float64}}  # Reference coordinates (t=0)
    reference_com::SVector{3,Float64}              # Reference center of mass
    reference_coords_com_frame::Vector{SVector{3,Float64}}  # Reference positions in COM frame
    # Per-ring reference data (for double ring systems)
    reference_coords_1::Vector{SVector{3,Float64}}
    reference_coords_2::Vector{SVector{3,Float64}}
    reference_com_1::SVector{3,Float64}
    reference_com_2::SVector{3,Float64}
    reference_coords_com_frame_1::Vector{SVector{3,Float64}}
    reference_coords_com_frame_2::Vector{SVector{3,Float64}}
    # Combined MSD timeseries (all monomers)
    msd_monomer::Vector{Float64}                   # Monomer MSD timeseries
    msd_com::Vector{Float64}                       # COM MSD timeseries
    msd_com_frame::Vector{Float64}                 # COM-frame MSD timeseries
    # Per-ring MSD timeseries (for double ring systems)
    msd_monomer_1::Vector{Float64}
    msd_monomer_2::Vector{Float64}
    msd_com_1::Vector{Float64}
    msd_com_2::Vector{Float64}
    msd_com_frame_1::Vector{Float64}
    msd_com_frame_2::Vector{Float64}
end

"""
    MSDLogger(schedule, system_type, n_monomers_1, n_monomers_2, initial_coords, boundary;
              compute_com=true, compute_com_frame=false)

Create MSDLogger with reference coordinates from the first frame.
For double ring systems, also initializes per-ring reference data.
"""
function MSDLogger(schedule::Vector{Int}, system_type::Symbol, n_monomers_1::Int64,
                   n_monomers_2::Int64, initial_coords::Vector{SVector{3,Float64}},
                   boundary; compute_com::Bool=true, compute_com_frame::Bool=false)
    # Store reference coordinates (unwrapped)
    n_total = system_type == :single ? n_monomers_1 : n_monomers_1 + n_monomers_2
    ref_coords = unwrap_polymer(initial_coords[1:n_total], boundary)

    # Compute reference COM
    ref_com = (compute_com || compute_com_frame) ? sum(ref_coords) / length(ref_coords) : SVector{3,Float64}(0, 0, 0)

    # Compute reference coordinates in COM frame
    ref_coords_com_frame = compute_com_frame ? ref_coords .- Ref(ref_com) : SVector{3,Float64}[]

    # Per-ring reference data (for double ring systems)
    if system_type == :double
        ref_coords_1 = unwrap_polymer(initial_coords[1:n_monomers_1], boundary)
        ref_coords_2 = unwrap_polymer(initial_coords[n_monomers_1+1:n_total], boundary)

        ref_com_1 = (compute_com || compute_com_frame) ? sum(ref_coords_1) / length(ref_coords_1) : SVector{3,Float64}(0, 0, 0)
        ref_com_2 = (compute_com || compute_com_frame) ? sum(ref_coords_2) / length(ref_coords_2) : SVector{3,Float64}(0, 0, 0)

        ref_coords_com_frame_1 = compute_com_frame ? ref_coords_1 .- Ref(ref_com_1) : SVector{3,Float64}[]
        ref_coords_com_frame_2 = compute_com_frame ? ref_coords_2 .- Ref(ref_com_2) : SVector{3,Float64}[]
    else
        ref_coords_1 = SVector{3,Float64}[]
        ref_coords_2 = SVector{3,Float64}[]
        ref_com_1 = SVector{3,Float64}(0, 0, 0)
        ref_com_2 = SVector{3,Float64}(0, 0, 0)
        ref_coords_com_frame_1 = SVector{3,Float64}[]
        ref_coords_com_frame_2 = SVector{3,Float64}[]
    end

    MSDLogger(
        Set{Int64}(schedule), Int64[], compute_com, compute_com_frame, system_type, n_monomers_1, n_monomers_2,
        # Combined reference data
        ref_coords, ref_com, ref_coords_com_frame,
        # Per-ring reference data
        ref_coords_1, ref_coords_2, ref_com_1, ref_com_2, ref_coords_com_frame_1, ref_coords_com_frame_2,
        # Combined MSD timeseries
        Float64[], Float64[], Float64[],
        # Per-ring MSD timeseries
        Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]
    )
end

"""
Log MSD values at each step.
For double ring systems, computes MSD separately for each ring in addition to combined MSD.
"""
function Molly.log_property!(logger::MSDLogger, sys, buffers, neighbors=nothing,
                             step_n::Integer=0; n_threads=1, kwargs...)
    if !(step_n in logger.log_steps)
        return
    end
    push!(logger.step_indices, step_n)

    n_total = logger.system_type == :single ? logger.n_monomers_1 :
              logger.n_monomers_1 + logger.n_monomers_2

    # Get current coordinates (unwrapped)
    current_coords = unwrap_polymer(sys.coords[1:n_total], sys.boundary)

    # Compute monomer MSD (combined)
    diff = current_coords .- logger.reference_coords
    msd_mon = sum(d -> dot(d, d), diff) / length(diff)
    push!(logger.msd_monomer, msd_mon)

    # Compute current COM (needed for COM MSD and COM-frame MSD)
    current_com = (logger.compute_com || logger.compute_com_frame) ?
                  sum(current_coords) / length(current_coords) : SVector{3,Float64}(0, 0, 0)

    # Compute COM MSD (skip if disabled)
    if logger.compute_com
        com_diff = current_com - logger.reference_com
        msd_com_val = dot(com_diff, com_diff)
        push!(logger.msd_com, msd_com_val)
    end

    # Compute COM-frame MSD (skip if disabled)
    if logger.compute_com_frame
        # Transform current coordinates to COM frame
        current_coords_com_frame = current_coords .- Ref(current_com)
        # Compute MSD in COM frame
        diff_com_frame = current_coords_com_frame .- logger.reference_coords_com_frame
        msd_com_frame_val = sum(d -> dot(d, d), diff_com_frame) / length(diff_com_frame)
        push!(logger.msd_com_frame, msd_com_frame_val)
    end

    # Per-ring MSD for double ring systems
    if logger.system_type == :double
        # Get per-ring coordinates
        current_coords_1 = unwrap_polymer(sys.coords[1:logger.n_monomers_1], sys.boundary)
        current_coords_2 = unwrap_polymer(sys.coords[logger.n_monomers_1+1:n_total], sys.boundary)

        # Monomer MSD per ring
        diff_1 = current_coords_1 .- logger.reference_coords_1
        diff_2 = current_coords_2 .- logger.reference_coords_2
        msd_mon_1 = sum(d -> dot(d, d), diff_1) / length(diff_1)
        msd_mon_2 = sum(d -> dot(d, d), diff_2) / length(diff_2)
        push!(logger.msd_monomer_1, msd_mon_1)
        push!(logger.msd_monomer_2, msd_mon_2)

        # Per-ring COM and COM MSD
        if logger.compute_com
            current_com_1 = sum(current_coords_1) / length(current_coords_1)
            current_com_2 = sum(current_coords_2) / length(current_coords_2)
            com_diff_1 = current_com_1 - logger.reference_com_1
            com_diff_2 = current_com_2 - logger.reference_com_2
            push!(logger.msd_com_1, dot(com_diff_1, com_diff_1))
            push!(logger.msd_com_2, dot(com_diff_2, com_diff_2))
        end

        # Per-ring COM-frame MSD
        if logger.compute_com_frame
            current_com_1 = sum(current_coords_1) / length(current_coords_1)
            current_com_2 = sum(current_coords_2) / length(current_coords_2)
            current_coords_com_frame_1 = current_coords_1 .- Ref(current_com_1)
            current_coords_com_frame_2 = current_coords_2 .- Ref(current_com_2)
            diff_com_frame_1 = current_coords_com_frame_1 .- logger.reference_coords_com_frame_1
            diff_com_frame_2 = current_coords_com_frame_2 .- logger.reference_coords_com_frame_2
            push!(logger.msd_com_frame_1, sum(d -> dot(d, d), diff_com_frame_1) / length(diff_com_frame_1))
            push!(logger.msd_com_frame_2, sum(d -> dot(d, d), diff_com_frame_2) / length(diff_com_frame_2))
        end
    end
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