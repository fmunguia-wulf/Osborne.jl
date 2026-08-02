"""
    pfile2imas!(p::PFile, dd::IMASdd.dd; gfile::Union{EFIT.GEQDSKFile,Nothing}=nothing, rho_target::Union{Nothing,AbstractVector{<:Real}}=nothing)

Populate dd.core_profiles from a PFile (Julia translation of OMFIT's
OMFITpFile.to_omas()), mapping psinorm -> rho_tor_norm via gfile.rhovn.

If `gfile` is omitted, falls back to dd.equilibrium.time_slice[]'s own
psi_norm/rho_tor_norm mapping.

`rho_target` sets the output grid. If omitted: reuses dd's own pre-existing
core_profiles grid if present (so the result matches ActorReplay.replay_dd's
length requirement automatically), else defaults to a 129-point grid.
"""
function pfile2imas!(p::PFile, dd::IMASdd.dd; gfile::Union{EFIT.GEQDSKFile,Nothing}=nothing, rho_target::Union{Nothing,AbstractVector{<:Real}}=nothing)

    if gfile !== nothing
        sign_Ip = sign(gfile.current)
        psin_eq = collect(range(0.0, 1.0, length=gfile.nw))
        rhotn_eq = collect(gfile.rhovn)
    else
        @info "No gfile provided. Falling back to equilibrium information stored inside of dd.equilibrium.time_slice[]"
        eqt = dd.equilibrium.time_slice[]
        psin_eq = collect(eqt.profiles_1d.psi_norm)
        rhotn_eq = collect(eqt.profiles_1d.rho_tor_norm)
        sign_Ip = sign(eqt.global_quantities.ip)
    end

    # helper: interpolate pfile psinorm onto equilibrium rho_tor_norm grid
    function psin2rho(psin_arr::Vector{Float64})
        clamped = clamp.(psin_arr, minimum(psin_eq), maximum(psin_eq))
        return linear_interp(psin_eq, rhotn_eq, clamped; extrap=FastInterpolations.ClampExtrap())
    end

    # Regrid in rho-space, not psinorm-space (OMFITpFile.remap()'s own
    # approach made edge Te collapse toward zero on this data).
    rho_target = if rho_target !== nothing
        collect(Float64.(rho_target))
    elseif !isempty(dd.core_profiles.profiles_1d)
        collect(dd.core_profiles.profiles_1d[end].grid.rho_tor_norm)
    else
        collect(range(0.0, 1.0, 129))
    end

    function get_prof(key::String, fac::Float64)
        haskey(p.profiles, key) || return nothing, nothing
        prof = p.profiles[key]
        # drop SOL-extension points (psinorm > 1) — otherwise they all
        # collapse onto the single point rho_tor_norm=1, corrupting the edge
        idx = findall(<=(1.0), prof.psinorm)
        rho = psin2rho(prof.psinorm[idx])
        vals = fac .* prof.data[idx]
        clamped_target = clamp.(rho_target, minimum(rho), maximum(rho))
        return rho_target, linear_interp(rho, vals, clamped_target; extrap=FastInterpolations.ClampExtrap())
    end

    # set target time from the user input dd
    time_target = !isempty(dd.equilibrium.time) ? dd.equilibrium.time[end] : p.time

    dd.global_time = time_target

    resize!(dd.core_profiles.profiles1d)
    cp1d = dd.core_profiles.profiles1d[]
    cp1d.time = time_target

    cp1d.grid.rho_tor_norm = rho_target

    if haskey(p.profiles, "ne")
        _, vals = get_prof("ne", 1e20)
        cp1d.electrons.density_thermal = vals
    end

    if haskey(p.profiles, "omgeb")
        _, vals = get_prof("omgeb", sign_Ip * 1e3)
        cp1d.rotation_frequency_tor_sonic = vals
    end

    if haskey(p.profiles, "te")
        _, vals = get_prof("te", 1e3)
        cp1d.electrons.temperature = vals
    end

    if haskey(p.profiles, "ptot")
        _, vals = get_prof("ptot", 1e3 / 3.0)
        cp1d.pressure_perpendicular = vals
        cp1d.pressure_parallel = vals
    end

    if haskey(p.profiles, "er")
        _, vals = get_prof("er", 1e3)
        cp1d.e_field.radial = vals
    end

    # ── ion species ────────────────────────────────────────────────────────
    if !isnothing(p.ion_species) && length(p.ion_species.A) >= 2
        # p-file's N Z A block orders [impurities..., main, beam]; output
        # cp1d.ion[] must put the main ion FIRST (FUSE's own convention) —
        # ActorPedestal's replay blend zips ion arrays positionally, not by
        # label. Main+beam are consolidated onto one ion[] entry.
        n_ions = length(p.ion_species.A)
        n_imp = n_ions - 2
        mk = 1

        resize!(cp1d.ion, n_imp + 1)

        main_Z = p.ion_species.Z[end-1]
        main_A = p.ion_species.A[end-1]
        cp1d.ion[mk].label = main_Z == 1.0 ? "D" : "H"
        resize!(cp1d.ion[mk].element, 1)
        cp1d.ion[mk].element[1].a = main_A
        cp1d.ion[mk].element[1].z_n = main_Z
        cp1d.ion[mk].multiple_states_flag = 0

        if haskey(p.profiles, "ni")
            _, vals = get_prof("ni", 1e20)
            cp1d.ion[mk].density_thermal = vals
        end
        if haskey(p.profiles, "ti")
            _, vals = get_prof("ti", 1e3)
            cp1d.ion[mk].temperature = vals
        end
        if haskey(p.profiles, "vtor1")
            _, vals = get_prof("vtor1", 1e3)
            cp1d.ion[mk].velocity.toroidal = vals
        end
        if haskey(p.profiles, "vpol1")
            _, vals = get_prof("vpol1", 1e3)
            cp1d.ion[mk].velocity.poloidal = vals
        end
        # sonic2ωtor! applies the diamagnetic correction (not equal to
        # rotation_frequency_tor_sonic); populated explicitly since
        # ActorFluxMatcher's replay step requires hasdata(), which a
        # dynamic expression alone doesn't satisfy.
        if haskey(p.profiles, "omgeb")
            IMAS.sonic2ωtor!(cp1d, cp1d.ion[mk])
        end

        if haskey(p.profiles, "nb")
            _, vals = get_prof("nb", 1e20)
            cp1d.ion[mk].density_fast = vals
        end
        if haskey(p.profiles, "pb")
            _, vals = get_prof("pb", 1e3 / 3.0)
            cp1d.ion[mk].pressure_fast_perpendicular = vals
            cp1d.ion[mk].pressure_fast_parallel = vals
        end

        # impurity ions: indices 2..n_imp+1 (main ion occupies index 1)
        for i in 1:n_imp
            imp_Z = p.ion_species.Z[i]
            imp_A = p.ion_species.A[i]
            k = i + 1
            label = if imp_Z == 6.0
                "C12"
            elseif imp_Z == 1.0
                imp_A <= 1.5 ? "H" : "D"
            else
                "Z$(round(Int, imp_Z))"
            end
            cp1d.ion[k].label = label
            resize!(cp1d.ion[k].element, 1)
            cp1d.ion[k].element[1].a = imp_A
            cp1d.ion[k].element[1].z_n = imp_Z
            cp1d.ion[k].multiple_states_flag = 0

            # density: nz1, nz2, ...
            nz_key = "nz$i"
            if haskey(p.profiles, nz_key)
                _, vals = get_prof(nz_key, 1e20)
                cp1d.ion[k].density_thermal = vals
            end
            # all ions assumed same temperature (OMFIT comment)
            if haskey(p.profiles, "ti")
                _, vals = get_prof("ti", 1e3)
                cp1d.ion[k].temperature = vals
            end
        end

    else
        # ── fallback: no N Z A block — assume D + C (OMFIT fallback) ─────────
        @warn "N Z A missing from pFile; assuming deuterium main ion and carbon impurity"
        resize!(cp1d.ion, 2)

        # main ion: D
        cp1d.ion[1].label = "D"
        resize!(cp1d.ion[1].element, 1)
        cp1d.ion[1].element[1].a = 2.0
        cp1d.ion[1].element[1].z_n = 1.0
        cp1d.ion[1].multiple_states_flag = 0
        if haskey(p.profiles, "ni")
            _, v = get_prof("ni", 1e20)
            cp1d.ion[1].density_thermal = v
        end
        if haskey(p.profiles, "ti")
            _, v = get_prof("ti", 1e3)
            cp1d.ion[1].temperature = v
        end
        if haskey(p.profiles, "vtor1")
            _, v = get_prof("vtor1", 1e3)
            cp1d.ion[1].velocity.toroidal = v
        end
        if haskey(p.profiles, "vpol1")
            _, v = get_prof("vpol1", 1e3)
            cp1d.ion[1].velocity.poloidal = v
        end

        # impurity: C12 (trace density from ni if nz1 missing — OMFIT fallback)
        cp1d.ion[2].label = "C12"
        resize!(cp1d.ion[2].element, 1)
        cp1d.ion[2].element[1].a = 12.0
        cp1d.ion[2].element[1].z_n = 6.0
        cp1d.ion[2].multiple_states_flag = 0
        if haskey(p.profiles, "nz1")
            _, v = get_prof("nz1", 1e20)
            cp1d.ion[2].density_thermal = v
        elseif haskey(p.profiles, "ni")
            _, v = get_prof("ni", 1e-6)
            cp1d.ion[2].density_thermal = v  # trace
        end
        if haskey(p.profiles, "ti")
            _, v = get_prof("ti", 1e3)
            cp1d.ion[2].temperature = v
        end
    end

    return dd
end
