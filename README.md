# Osborne

[![Build Status](https://github.com/fmunguia-wulf/Osborne.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/fmunguia-wulf/Osborne.jl/actions/workflows/CI.yml?query=branch%3Amain)

[Osborne pedestal-fit files](https://fusion.gat.com/theory/Pfile) ("p-files") are
tanh-fit pedestal profiles produced by OMFIT's Osborne fitting tool, giving
smoother edge density/temperature/rotation profiles than raw ZIPFIT data.

`Osborne.jl` reads p-files and converts them into an
[IMAS](https://github.com/ProjectTorreyPines/IMAS.jl) `dd.core_profiles`
structure, playing the same role for p-files that
[`EFIT.jl`](https://github.com/JuliaFusion/EFIT.jl) plays for g-files (EFIT
GEQDSK equilibrium reconstructions): both are native-Julia readers that parse
a DIII-D file format directly into IMAS data structures, avoiding a round
trip through OMFIT/Python.

```julia-repl
julia> using Osborne, EFIT, IMAS

julia> p = readp("p190904.01906_324");

julia> g = EFIT.readg("g190904.01906_324");

julia> dd = IMAS.dd();

julia> EFIT.geqdsk2imas!([g], dd; add_derived=true);

julia> Osborne.pfile2imas!(p, dd; gfile=g);

julia> dd.core_profiles.profiles_1d[].electrons.temperature
```
