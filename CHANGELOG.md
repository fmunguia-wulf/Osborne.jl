# Changelog

All notable changes to this project are documented here.
Format based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This project follows [Semantic Versioning](https://semver.org/).

## [1.0.0] - 2026-07-31

Initial release.

### Added
- `readp`: parse Osborne pedestal-fit p-files.
- `pfile2imas!`: convert a parsed p-file into `dd.core_profiles`, matching
  OMFIT's `OMFITpFile.to_omas()` conventions (psinorm -> rho_tor_norm via
  g-file RHOVN, unit conversions, ion species handling).
- Main ion placed first in `cp1d.ion[]`, impurities after — matches FUSE's
  own convention so replay-based actors (e.g. `ActorPedestal`) pair species
  correctly.
- Main and beam populations of the same species consolidated onto one
  `cp1d.ion[]` entry (`density_thermal` + `density_fast`).
- `rho_target` auto-detection: reuses a pre-existing `dd.core_profiles` grid
  when present, so results are ready to use as `ActorReplay.replay_dd`
  without manual grid-matching.
