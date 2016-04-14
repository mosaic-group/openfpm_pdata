# Change Log
All notable changes to this project will be documented in this file.

## [0.3.0] -

### Added
- Molacular Dynamic example

### Fixed
- Nothing to report

### Changed
- Eliminated global_v_cluster, init_global_v_cluster, delete_global_v_cluster, 
  substituted by 
  create_vcluster, openfpm_init, openfpm_delete
- CartDecomposition parameter for the distributed structures is now optional
- template getPos<0>(), substituted by getPos()

## [0.2.1] - 

### Changed
- GoogleChart name function changed: AddPointGraph to AddLinesGraph and AddColumsGraph to AddHistGraph

## [0.2.0] - 2016-03-25
### Added
- Added Load Balancing and Dynamic Load Balancing on Beta
- PSE 1D example with multiple precision
- Plot example for GoogleChart plotting
- Distributed data structure now support 128bit floating point precision (on Beta)

### Fixed
- Detection 32 bit system and report as an error
- Bug in rounding off for periodic boundary condition

### Changed
- Nothing to report

## [0.1.0] - 2016-02-05
### Added
- PSE 1D example
- Cell list example
- Verlet list example
- Kickstart for OpenFPM_numeric
- Automated dependency installation for SUITESPRASE EIGEN OPENBLAS(LAPACK)


### Fixed
- CRITICAL BUG in periodic bondary condition
- BOOST auto updated to 1.60
- Compilation with multiple .cpp files

### Changed
- Nothing to report

