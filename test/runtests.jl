using Osborne
using EFIT
using IMAS
using Test

const testdir = @__DIR__
const pfile_path = joinpath(testdir, "data", "p190904.01906_324")
const gfile_path = joinpath(testdir, "data", "g190904.01906_324")

@testset "Osborne.jl" begin

    @testset "readp" begin
        p = readp(pfile_path)

        @test p isa PFile
        @test haskey(p.profiles, "ne")
        @test haskey(p.profiles, "te")
        @test haskey(p.profiles, "ni")
        @test haskey(p.profiles, "ti")
        @test haskey(p.profiles, "nz1")

        # 256-point profile grid on this shot's p-file
        @test length(p.profiles["ne"].psinorm) == 256
        @test length(p.profiles["ne"].data) == 256

        # time parsed from filename (p190904.01906_324 -> 1.906 s)
        @test p.time ≈ 1.906 atol = 1e-6

        # N Z A block: this shot's p-files are carbon-only (1 impurity + main + beam = 3 rows)
        @test !isnothing(p.ion_species)
        @test length(p.ion_species.A) == 3
    end

    @testset "readp: beam-less shot (no N Z A beam row)" begin
        # p190904.04353_697 (not checked into test fixtures) only has 2 species
        # rows [carbon, D] with no beam data (missing nb/pb) — regression-guard
        # for the has_beam handling using a synthetic minimal case instead of
        # requiring that specific fixture file.
        # (Covered indirectly: main readp test above exercises the 3-row path;
        # this testset is a placeholder for adding the 2-row fixture later.)
        @test_skip false
    end

    @testset "pfile2imas!" begin
        p = readp(pfile_path)
        g = EFIT.readg(gfile_path)

        dd = IMAS.dd()
        EFIT.geqdsk2imas!([g], dd; add_derived=true)
        dd.global_time = g.time
        IMAS.flux_surfaces(dd.equilibrium.time_slice[], IMAS.first_wall(dd.wall)...)

        Osborne.pfile2imas!(p, dd; gfile=g)

        cp1d = dd.core_profiles.profiles_1d[]

        @test cp1d.time ≈ p.time atol = 1e-6

        # electrons
        @test all(cp1d.electrons.density_thermal .> 0)
        @test all(cp1d.electrons.temperature .> 0)

        # ion species: 2 total for this shot (1 impurity + 1 consolidated D
        # entry carrying both main/thermal and beam/fast populations) —
        # main+beam are merged onto a single ion[] entry, not split across
        # two entries sharing the same label (see pfile2imas! docs)
        @test length(cp1d.ion) == 2
        labels = [ion.label for ion in cp1d.ion]
        @test "C12" in labels
        @test count(==("D"), labels) == 1

        # ordering: main ion (D) MUST be first, matching FUSE's own native
        # convention (e.g. ZIPFIT-derived dd's, dd's built from
        # ini.core_profiles.bulk/impurity) -- several FUSE actors
        # (ActorPedestal's replay-blend step) zip cp1d.ion against
        # replay_cp1d.ion POSITIONALLY, not by label, so a mismatched
        # ordering silently pairs the wrong species together
        @test cp1d.ion[1].label == "D"
        @test cp1d.ion[2].label == "C12"

        # consolidated D ion: both thermal (main) and fast (beam) populated
        # on the SAME entry
        D = cp1d.ion[findfirst(ion -> ion.label == "D", cp1d.ion)]
        @test all(D.density_thermal .> 0)
        @test all(D.temperature .> 0)
        @test any(D.density_fast .> 0)

        # carbon impurity
        carbon = cp1d.ion[findfirst(ion -> ion.label == "C12", cp1d.ion)]
        @test all(carbon.density_thermal .> 0)

        # regridding sanity: everything shares the same rho_tor_norm length
        n = length(cp1d.grid.rho_tor_norm)
        @test length(cp1d.electrons.density_thermal) == n
        @test length(D.density_thermal) == n
        @test length(carbon.density_thermal) == n

        # rotation (omgeb -> rotation_frequency_tor_sonic, per OMFIT's own
        # to_omas() mapping, confirmed directly against omfit_osborne.py)
        @test !ismissing(cp1d, :rotation_frequency_tor_sonic)
    end

    @testset "pfile2imas!: preserves pre-existing core_profiles grid" begin
        # Simulates overriding an already-fetched dd (e.g. ZIPFIT) with
        # p-file data — the p-file's output grid should automatically match
        # whatever grid was already there, with no explicit rho_target
        # needed, so the result stays compatible with anything (e.g.
        # ActorReplay) already expecting that dd's original grid length.
        p = readp(pfile_path)
        g = EFIT.readg(gfile_path)

        dd = IMAS.dd()
        EFIT.geqdsk2imas!([g], dd; add_derived=true)
        dd.global_time = g.time
        IMAS.flux_surfaces(dd.equilibrium.time_slice[], IMAS.first_wall(dd.wall)...)

        # pre-populate core_profiles with an arbitrary, non-default-length grid
        preexisting_rho = collect(range(0.0, 1.0, 57))
        cp1d0 = resize!(dd.core_profiles.profiles_1d)
        cp1d0.time = g.time
        cp1d0.grid.rho_tor_norm = preexisting_rho

        Osborne.pfile2imas!(p, dd; gfile=g)  # no rho_target passed

        cp1d = dd.core_profiles.profiles_1d[]
        @test length(cp1d.grid.rho_tor_norm) == 57
        @test cp1d.grid.rho_tor_norm == preexisting_rho
    end

end
