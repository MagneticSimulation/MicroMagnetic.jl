# Release notes

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).



## Version [v0.5.0]

### Added

- Web based GUI

### Changed

- change relax_m_interval to relax_m_every in sim_with function
- change dynamic_m_interval to dynamic_m_every in sim_with function
- change relax_data_interval to relax_data_every in sim_with function

### Removed

- add_exch(sim::MicroSim, geo::Shape, A::Number; name="exch")
- set_Ms(sim::AbstractSim, shape::Union{CSGNode,Shape}, Ms::Number)

## Version [v0.4.0]

### Added
- voronoi
- FEM (exch, demag)
- add_stt, add_sot
- inertia llg

### Changed
- she_torque -> sahe_torque
- EnergyMinimization -> SD

### Removed
- JLD2 

## Version [v0.3.9] - 2025-02-16

### Added

* Eigen methods (not finished yet)
* Hexagonal anisotropy
* SHE torque
* Biquadratic exchange

### Fixed

* mu0 for stochastic field