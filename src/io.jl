# 1D profiles
function interpolate_1d_profile(psi_norm::AbstractVector{<:Real}, ori_1D::AbstractVector{<:Real}, N::Int)
    if length(ori_1D) == N
        return ori_1D
    else
        return pchip_interp(psi_norm, ori_1D, range(0, 1, N))
    end
end

# ── Data structures ───────────────────────────────────────────────────────────

struct PFileProfile
    psinorm::Vector{Float64}
    data::Vector{Float64}
    derivative::Vector{Float64}
    units::String
    description::String
end

struct PFileIonSpecies
    N::Vector{Float64}   # number of nucleons
    Z::Vector{Float64}   # atomic number
    A::Vector{Float64}   # atomic mass
end

mutable struct PFile
    file::String
    time::Float64
    profiles::Dict{String,PFileProfile}
    ion_species::Union{PFileIonSpecies,Nothing}
end

function Base.show(io::IO, p::PFile)
    print(io, "PFile: \"", p.file, "\" t=", p.time, "s profiles: ", join(sort(collect(keys(p.profiles))), ", "))
end

# ── Reader — faithful translation of OMFITpFile.load() ───────────────────────

"""
    readp(pfile::String; set_time=nothing)

Read an Osborne pfile and return a PFile struct.

Units as stored in the pfile (same as OMFIT):
 - ne, ni, nb, nz : 10^20/m^3
 - te, ti          : keV
 - pb, ptot        : kPa
 - omeg, omegp, omgvb, omgpp, omgeb, omghb : kRad/s
 - er              : kV/m
 - vtor, vpol      : km/s
 - kpol            : km/s/T
"""
function readp(pfile::String; set_time=nothing)
    isfile(pfile) || error("$(pfile) does not exist")

    profiles = Dict{String,PFileProfile}()
    ion_species = nothing

    # parse time from filename  e.g. p190904.01906_324 -> 1.906 s
    time = if !isnothing(set_time)
        Float64(set_time)
    else
        try
            bn = basename(pfile)                       # "p190904.01906_324"
            no_ext = split(bn, ".")[end]                # "01906_324"
            time_str = split(no_ext, "_")[1]            # "01906"
            parse(Float64, time_str) / 1000.0
        catch
            0.0
        end
    end

    lines = readlines(pfile)
    ind = 1

    while ind <= length(lines)
        line = lines[ind]
        cur = split(strip(line))
        isempty(cur) && (ind += 1; continue)

        # first token must be an integer (number of data rows)
        num = tryparse(Int, cur[1])
        isnothing(num) && (ind += 1; continue)

        # ── N Z A of ION SPECIES block ────────────────────────────────────────
        if occursin("N Z A of ION SPECIES", line)
            Ns = Float64[]
            Zs = Float64[]
            As = Float64[]
            for j in ind+1:ind+num
                j > length(lines) && break
                row = parse.(Float64, split(strip(lines[j])))
                length(row) >= 3 || continue
                push!(Ns, row[1])
                push!(Zs, row[2])
                push!(As, row[3])
            end
            ion_species = PFileIonSpecies(Ns, Zs, As)
            ind += num + 1
            continue
        end

        # ── profile block ─────────────────────────────────────────────────────
        # header format (from OMFIT regex):
        #   "256 psinorm ne(10^20/m^3) dne/dpsiN"
        #    num  xkey    key(units)   der
        m = match(r"^(\d+)\s+(\S+)\s+(\w+)\(([^)]*)\)\s+(\S+)", line)
        if isnothing(m)
            ind += 1
            continue
        end

        # xkey = m[2]  (always "psinorm")
        key = m[3]
        units = m[4]
        # der = m[5]

        psinorm = Float64[]
        data = Float64[]
        deriv = Float64[]
        for j in ind+1:ind+num
            j > length(lines) && break
            row = tryparse.(Float64, split(strip(lines[j])))
            any(isnothing, row) || length(row) < 3 && continue
            push!(psinorm, row[1])
            push!(data, row[2])
            push!(deriv, row[3])
        end

        if !isempty(data)
            desc_map = Dict(
                "ne" => "Electron density",
                "te" => "Electron temperature",
                "ni" => "Ion density",
                "ti" => "Ion temperature",
                "nb" => "Fast ion density",
                "pb" => "Fast ion pressure",
                "ptot" => "Total pressure",
                "omeg" => "Toroidal rotation VTOR/R",
                "omegp" => "Poloidal rotation Bt*VPOL/(RBp)",
                "omgvb" => "VxB rotation term OMEG+OMEGP",
                "omgpp" => "Diamagnetic term in ExB",
                "omgeb" => "ExB rotation frequency",
                "er" => "Radial electric field",
                "ommvb" => "Main ion VxB term",
                "ommpp" => "Main ion pressure term",
                "omevb" => "Electron VxB term",
                "omepp" => "Electron pressure term",
                "kpol" => "KPOL = VPOL/Bp",
                "omghb" => "Hahm-Burrell ExB shearing rate",
                "vtor1" => "Toroidal velocity 1st impurity",
                "vpol1" => "Poloidal velocity 1st impurity",
                "nz1" => "Density 1st impurity",
                "vtor2" => "Toroidal velocity 2nd impurity",
                "vpol2" => "Poloidal velocity 2nd impurity",
                "nz2" => "Density 2nd impurity",
            )
            desc = get(desc_map, key, key)
            profiles[key] = PFileProfile(psinorm, data, deriv, units, desc)
        end

        ind += num + 1
    end

    return PFile(basename(pfile), time, profiles, ion_species)
end
