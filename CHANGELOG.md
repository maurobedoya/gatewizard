# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- 

### Changed
-

### Fixed
-

### Removed

## [1.0.19] - 2026-02-12

### Fixed
- Fixed equilibration bug where bilayer PDB with CRYST1 record was not being copied to output directory
- Fixed NAMDEquilibrationManager initialization to use correct directory (output_dir where files are copied)
- Box dimensions now correctly read from bilayer_*_lipid.pdb file for NAMD equilibration configurations

## [1.0.18] - 2026-02-12

### Fixed
- Fixed missing equilibration templates in package distribution by correcting MANIFEST.in
- Equilibration NAMD template files (.inp) are now properly included when installing the package

## [1.0.17] - 2026-01-21

### Fixed
- Corrected water model availability: replaced non-existent `tip4p` with `tip4pd` (TIP4P-D)
- Fixed leaprc file mappings for all water models (tip4pd, tip4pew, opc3, spceb, fb3)
- Each water model now correctly maps to its specific leaprc file

## [1.0.16] - 2026-01-21

### Changed
- Replaced Unicode icons with basic ASCII symbols for better terminal compatibility
- Warning icon (⚠️) replaced with "WARNING:" text
- Bullet points (•) replaced with dashes (-)
- Check marks (✓) replaced with "[OK]"
- Cross marks (✗) replaced with "[ERROR]"

## [1.0.15] - 2026-01-21

### Changed
- Force field validation now shows warnings instead of blocking preparation for unvalidated combinations
- Updated force field compatibility matrix based on literature validation references
- Added comprehensive validation references ([1-6]) to force field module with DOI links
- Updated water model compatibility: TIP3P validated with ff14SB+lipid17/lipid21, OPC with ff19SB+lipid21
- ff15ipq now correctly marked as protein-only (not validated with lipid force fields)
- Updated recommendations: membrane (TIP3P+ff14SB+lipid21), latest (OPC+ff19SB+lipid21)

### Fixed
- Force field combinations now properly categorized as validated, recommended, or unvalidated
- Users can now test experimental force field combinations at their own risk with clear warnings

## [1.0.14] - 2026-01-15
### Changed
- Documentation update

## [1.0.13] - 2026-01-15

### Added
- Zenodo link creation

### Fixed
- PyPI badged was fixed

## [1.0.12] - 2026-01-14

### Fixed
- Fixed PyPI publishing workflow by removing duplicate publish job
- Corrected workflow configuration to properly trigger PyPI upload on release

## [1.0.11] - 2026

### Fixed
- Automatic release creation to GitHub and PyPI from tags

## [1.0.10] - 2026

### Fixed
- Automatic release creation from tags testing


## [1.0.9] - 2026

### Fixed
- Automatic release creation from tags testing


## [1.0.8] - 2026

### Fixed
- Automatic release creation from tags

## [1.0.7] - 2026

### Added
- Automatic release creation from tags
- CHANGELOG.md file to track project changes
- It was switched to PyPI Trusted Publishing

### Changed
- Updated test documentation to clarify pytest installation requirement

### Fixed
- NPAT equilibration protocol was updated to match the NPgT 

### Removed
- GateWizard version from the main GUI

## [1.0.6] - 2025

### Added
- CHANGELOG.md file to track project changes
- It was switched to PyPI Trusted Publishing

### Changed
- Updated test documentation to clarify pytest installation requirement

### Fixed
- NPAT equilibration protocol was updated to match the NPgT 

### Removed
- GateWizard version from the main GUI

## [1.0.5] - 2025

### Added
- CHANGELOG.md file to track project changes

### Changed
- Updated test documentation to clarify pytest installation requirement

### Fixed
- NPAT equilibration protocol was updated to match the NPgT 

### Removed
- GateWizard version from the main GUI

## [1.0.4] - 2025

### Added
- Initial release of GateWizard
- Modern GUI built with CustomTkinter
- Protein structure preparation and cleaning
- Propka integration for pKa calculations
- Protein capping with ACE/NME termini
- Amber force field support (ff14SB, ff19SB)
- Membrane system building
- NAMD equilibration protocol generation
- NAMD log file analysis tools
- Molecular viewer integration
- Comprehensive API documentation
- User guide and troubleshooting documentation

[Unreleased]: https://github.com/maurobedoya/gatewizard/compare/v1.0.19...HEAD
[1.0.19]: https://github.com/maurobedoya/gatewizard/compare/v1.0.18...v1.0.19
[1.0.18]: https://github.com/maurobedoya/gatewizard/compare/v1.0.17...v1.0.18
[1.0.17]: https://github.com/maurobedoya/gatewizard/compare/v1.0.16...v1.0.17
[1.0.16]: https://github.com/maurobedoya/gatewizard/compare/v1.0.15...v1.0.16
[1.0.15]: https://github.com/maurobedoya/gatewizard/compare/v1.0.14...v1.0.15
[1.0.14]: https://github.com/maurobedoya/gatewizard/compare/v1.0.13...v1.0.14
[1.0.13]: https://github.com/maurobedoya/gatewizard/compare/v1.0.12...v1.0.13
[1.0.12]: https://github.com/maurobedoya/gatewizard/compare/v1.0.11...v1.0.12
[1.0.11]: https://github.com/maurobedoya/gatewizard/compare/v1.0.10...v1.0.11
[1.0.10]: https://github.com/maurobedoya/gatewizard/compare/v1.0.9...v1.0.10
[1.0.9]: https://github.com/maurobedoya/gatewizard/compare/v1.0.8...v1.0.9
[1.0.8]: https://github.com/maurobedoya/gatewizard/compare/v1.0.7...v1.0.8
[1.0.7]: https://github.com/maurobedoya/gatewizard/compare/v1.0.6...v1.0.7
[1.0.6]: https://github.com/maurobedoya/gatewizard/compare/v1.0.5...v1.0.6
[1.0.0]: https://github.com/maurobedoya/gatewizard/releases/tag/v1.0.0