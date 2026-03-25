# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.25] - 2026-03-25

### Added
- **Transform Structure** dialog in the Visualize frame with three tabs: Rotate, Translate, and Align
  - Rotate: rotate selected or all atoms around X/Y/Z axis by arbitrary angle, with selectable pivot (selection centroid or origin)
  - Translate: translate atoms by displacement vector (X, Y, Z in Å), plus Center at Origin button
  - Align: align a selection's principal direction (SVD) to a target axis, with optional secondary axis alignment
  - MDAnalysis selection expressions supported for all operations
  - Preview shows transformed atom positions (yellow glow at destination) before applying
  - Non-modal dialog allows rotating the 3D view while the dialog is open
- `rotate_atoms()`, `translate_atoms()`, `center_atoms()`, `align_to_axis()` methods in `MolecularViewer` API
- `_axis_rotation_matrix()` and `_rotation_matrix_from_vectors()` helper functions in `gatewizard.core.viewer`
- `ProteinStructure.refresh_residue_coords()` method to sync residue CA/O coords after transforms
- `_reassign_ss()` / `_reassign_ss_gui()` helpers to recalculate secondary structure after coordinate changes
- Four new viewer examples: 12 (rotate), 13 (translate/center), 14 (align to axis), 15 (primary+secondary alignment)
- Coordinate Transformations section in API docs (`docs/api/viewer.md`) with Examples 12–15

### Fixed
- Secondary structure now updates correctly after all coordinate transformations (was stale because `Residue.ca_coord` / `o_coord` became detached after atom coord reassignment)
- `build_bonds()` now calls `refresh_residue_coords()` automatically so tube_ss, cartoon, and backbone representations reflect coordinate changes
- Auto-detect first selection (e.g. Protein) now shows delete (X) button — changed guard from `idx > 0` to `sel.name != "All"`
- Unicode arrow `→` replaced with `->` in align messagebox to avoid rendering issues in tkinter

## [1.0.24] - 2026-03-21

### Fixed
- `psique` executable now auto-sets executable permission on first use (pip strips execute bits from package data files)
- Improved psique error messages: distinguish "not found" from "no SS results for small structures"
- Updated viewer example 11 with real protein coordinates (2MVJ residues 1-20) for proper psique testing
- Added Viewer Module to mkdocs documentation navigation

## [1.0.23] - 2026-03-21

### Fixed
- Version bump to fix tag/release mismatch from v1.0.21/v1.0.22

## [1.0.20] - 2026-03-21

### Added
- VTK-based 3D molecular viewer replacing the old matplotlib Visualize frame
- New `gatewizard.core.viewer` module with `MolecularViewer` API for programmatic structure loading, inspection, selection, editing, and saving
- MDAnalysis-based PDB parsing (replaces BioPython dependency)
- Seven molecular representations: VDW, Ball & Stick, Sticks, Cartoon, Tube SS, Backbone, Surface
- Multiple named selections with independent representation, color scheme, quality, material, and SS color settings
- Selection criteria: All, Protein, Backbone, Sidechain, Water, Ligand, Chain, Residue range, Around selection
- Auto-detect molecules (protein, water, individual ligands) with sensible defaults
- Structure editing: rename chains, rename residues, renumber residues, delete atoms
- VTK offscreen rendering widget (`VTKFrame`) with mouse rotation/pan/zoom and fog/depth cueing
- SSAO ambient occlusion and shadow map rendering passes
- High-resolution image export (PNG, JPEG, TIFF, BMP) with configurable scale and transparent background
- Viewpoint save/load as JSON (camera, selections, rendering settings)
- Drag-reorder selections in the GUI panel
- Per-selection settings dialog with live preview (quality, sizes, material presets, SS colors)
- `vtk>=9.0.0` added to project dependencies
- `MolecularViewer` exported at package level (`from gatewizard import MolecularViewer`)
- `VTKFrame` widget exported from `gatewizard.gui.widgets`
- Viewer test suite (`test_viewer.py`) with 10 example scripts covering the full API
- API documentation in `docs/api/viewer.md`

### Changed
- Moved `psique` executable from `utils/` to `tools/` (with TODO for future pip import)
- Blocking `wait=True` option in `Builder.prepare_system()` so scripts can run multiple systems sequentially
- New `Builder.wait_for_completion(job_dir)` method for waiting on an already-launched job
- Ligand parametrization module (`gatewizard.tools.ligand_parametrization`) with full AMBER/GAFF2 workflow
- Automatic detection of non-standard (ligand) residues in PDB files via HETATM scanning
- Ligand extraction, antechamber atom typing, parmchk2 missing parameter generation, and tleap .lib creation
- 2D molecular structure visualization using RDKit
- GUI widget (`LigandParamWidget`) in the Builder frame for interactive ligand detection, viewing, and parametrization
- Builder integration: `--ligand_param` and `--gaff2` flags for packmol-memgen, GAFF2/ligand loading in tleap parametrization
- New `ligand_params` configuration key in Builder for passing parametrized ligand files
- Builder example tests 19–24 covering ligand detection, extraction, parametrization, command building, and 2D imaging
- API documentation for all ligand parametrization functions and classes in `docs/api/builder.md`

### Changed
- Builder `_build_command()` now appends ligand parameter flags when ligands are present
- Builder `_create_tleap_input()` loads GAFF2 and ligand .frcmod/.lib files before `loadPDB`
- Builder bash execution script includes ligand parameter loading in tleap section
- Updated `gatewizard.tools.__init__` exports with all ligand parametrization public API

### Fixed
- Fixed Pillow deprecation warning: replaced `Image.getdata()` with `get_flattened_data()` (Pillow >=12 compat)
- Test suite now auto-cleans `./systems/` output directories after example execution

### Fixed (repo)
- Resolved Git lock file issue blocking commits in GitHub Desktop
- Added `psique` executable to version control (previously untracked)

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

[Unreleased]: https://github.com/maurobedoya/gatewizard/compare/v1.0.25...HEAD
[1.0.25]: https://github.com/maurobedoya/gatewizard/compare/v1.0.24...v1.0.25
[1.0.24]: https://github.com/maurobedoya/gatewizard/compare/v1.0.23...v1.0.24
[1.0.23]: https://github.com/maurobedoya/gatewizard/compare/v1.0.20...v1.0.23
[1.0.20]: https://github.com/maurobedoya/gatewizard/compare/v1.0.19...v1.0.20
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