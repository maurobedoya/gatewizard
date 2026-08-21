# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- **Cluster probe:** collect Slurm partitions/nodes (`sinfo`) before the slower `module avail`, and batch path env queries into one SSH round-trip, so Run-on-cluster Resources can fill sooner.
- **Equilibration (testing):** experimental membrane schedule — heat/scaffold NVT, pack NPgT (γ=0), then production in the selected ensemble. Templates are `equilibration/{engine}/eq/` plus `production/{NVT,NPT,NPAT,NPgT}/`. Details in `equilibration/PROTOCOL.md`.
- **Equilibration defaults:** `get_default_stage_params` for Amber/NAMD/GROMACS/OpenMM now return the universal schedule; generators load `eq/` for mini+Eq1–6 and `production/{ensemble}` for production. Headers stamp the **stage** ensemble.
- **Equilibration pressure / surface tension:** OpenMM (`p_ref` / `p_tens`) and GROMACS (`ref_p`) now take stage `pressure` and `surface_tension` (dyn/cm; GROMACS NPgT converts γ→bar·nm ×10). Defaults remain 1 bar / 0 dyn/cm for packing.
- **Equilibration (Amber):** `ntwx` is substituted from stage `dcd_freq` (Eq6 / production 50000); `ioutfm=1` (NetCDF).
- **PlotSpec overlay:** shared y-label when all panels use the same one (structural Pub PNG keeps “RMSD (Å)” instead of “Multiple Properties”). Markers only on short overlay series.
- **PlotSpec grid:** panels may list `series_keys` to draw multiple sets on one subplot (energetic compare-by-property / by-set Pub PNG).
- **Publication plot export:** matplotlib uses the headless Agg backend in API/GUI export so Tk/Tcl is not touched from FastAPI worker threads (fixes `main thread is not in main loop` / `Tcl_AsyncDelete` log noise on WSL).

### Added

- **Builder:** `Builder.cancel_preparation` stops a running job via `process.pid` process-group kill and marks `status.json` as `cancelled`.
- **Tools Fix PBC (GROMACS):** multi-select center/output index groups merge into temporary `GW_CENTER` / `GW_OUTPUT` compound ndx entries; optional `skip_cluster`; smarter lipid multi-group recommendations when `SOLU_MEMB` is absent.
- **Cluster submit:** optional Slurm **GPU type** (`gpu_type`) for typed GRES — `#SBATCH --gres=gpu:TYPE:N` when set; untyped jobs still use `#SBATCH --gpus=N`. Types are parsed from probed node GRES.

### Fixed

- **PlotSpec overlay:** shared y-limits span every panel (union), and panel ylim falls back to global — structural APL Pub PNG no longer clips leaflets to the Mean-only window.
- **PlotSpec / Pub PNG:** keep panel `series_keys` through `normalize_plot_spec`, and use the panel **x** limits (not y) when `sync_x` is on — fixes empty energetic “one panel per set” publication PNGs where lines were missing or crushed into an invisible speck at t≈0.
- **Update manifest:** `releases/gui-versions.json` refreshed to GUI **1.0.13** / API **1.0.53** (was stuck at 1.0.11 / 1.0.49, so in-app update banners never appeared).
- **Equilibration (NAMD COM restraint):** insert `colvars on` / `colvarsConfig` **before** the first `minimize`/`run`. Appending them at the end of the conf caused `FATAL ERROR: Setting parameter colvars from script failed!` after step1 finished.
- **Equilibration Use in form / job metadata:** OpenMM (and other engines) recover **per-stage** ensembles from inputs — packing stays **NPgT**, and later NVT/NPT/NPAT stages no longer stick as NPT after the first barostat. Sticky recovered protocols in `equilibration_job.json` are healed on read so **Use in form** shows the real schedule.
- **Equilibration (NAMD):** `firsttimestep` now skips a folded Minimization stage and attributes its `minimize_steps` to Equilibration 1 (was writing `10000` on step2 instead of `135000` when the GUI protocol still listed Minimization).
- **Equilibration (all engines):** step6.2 / Amber step2 are NVT scaffold again; barostat onset is step6.3 / Amber step3 (classic CHARMM-GUI / original NPgT), not engine-specific.
- **Equilibration (Amber):** soft first-barostat (`taup=5.0`, `ntwr=5000`) and restraint `REF` refresh moved to `step3` (first packing barostat).
- **Equilibration resources (Amber):** Equilibration / production default to CPU×1 + GPU×1 (`pmemd.cuda`); minimization stays CPU-only. The **first packing barostat** stage (typically Eq3 / NPgT) defaults to CPU×6 `pmemd` to avoid GPU “box dimensions changed too much” aborts; later stages stay on GPU.
- **Equilibration (Amber):** `ensure_prmtop_box` runs whenever any stage uses NPT/NPAT/NPgT — including NVT-final protocols that still pack under NPgT — so `ifbox == 0` no longer blocks the first barostat.

## [1.0.53] - 2026-08-06

### Added

- **Analysis:** PlotSpec rendering helpers; structural per-file `file_strides`; energetic log stride helper (`energy_stride`)
- **Tools:** trajectory PBC fixing (`trajectory_tools` / `fix_pbc_worker`) with GROMACS / cpptraj / MDAnalysis paths
- **Equilibration resources:** per-stage CPU/GPU in `equilibration_resources.json` (engine-specific defaults; Amber picks `pmemd` vs `pmemd.cuda` per stage)
- **Equilibration (NAMD):** `gpu_resident` option on `setup_namd_equilibration` — writes `GPUresident` on the production stage only; equilibration keeps `reassignFreq`/`reassignTemp`; persisted on `equilibration_job.json`
- **Cluster mid-run progress:** sync `step*.log` (and related progress files) from node-local scratch → submit directory — batch scripts rsync every 60s; API helper SSHs to the allocated node for Watching/Pull on jobs already running
- **Cluster status:** job-status records allocated CPUs, node name, and node GPU type (from `sinfo` GRES) on `execution` for Watching cards
- **Equilibration job metadata:** if `equilibration_job.json` has only an `execution` block, infer protocol/ensemble/input_dir from `protocol_summary.json` and heal the job JSON (fixes **Use in form** after cluster Watching)
- **Equilibration (remote / Slurm):** `gatewizard.utils.cluster` — module/`sinfo` parsers, workdir strategies (scratch stage-in/out), editable batch templates, Slurm adapter (+ PBS stub), SSH/rsync helpers, and `execution` metadata on `equilibration_job.json`
- **Equilibration dual run scripts:** generate/refresh write `run_equilibration.sh` (local Executable) and `run_equilibration_cluster.sh` (module-friendly `namd3`/`gmx`/`python3`/`pmemd[.cuda]`); Slurm defaults to `bash run_equilibration_cluster.sh`
- **Equilibration failure detection:** scan stage logs / Slurm outs for FATAL / CUDA stub / Error in Stage (empty `.out` no longer hides NAMD failures)
- **Cluster probe:** list compute nodes via `sinfo -N`; batch templates support `#SBATCH --nodelist=` / `--constraint`
- **Equilibration (Amber):** full `AmberEquilibrationManager` — mdin templates for NVT/NPT/NPAT/NPgT, MDA GROUP positional restraints (no dihedrals), `run_equilibration.sh` with resume/resources, executable discovery (`pmemd.cuda` → `pmemd` → MPI → `sander`), and `amber_analysis` progress/energetic parsing
- **Equilibration templates:** generated inputs stamp GateWizard API version, local generation time (with timezone), and shared templates version (`testing` for now — protocols are still evolving) for traceability

### Changed

- **README:** GateWizard splash lockup at the top of the repo README
- **Equilibration (OpenMM):** default per-stage resources are CPU×1 + GPU×1 for minimization, equilibration, and production (OpenMM uses a single host thread)

### Fixed

- **Cluster upload:** Paramiko SFTP recursive put fails clearly on zero/partial uploads; `verify_remote_files` ensures the launched run script / `.slurm` exist before `sbatch`
- **Cluster probe:** reject help-text tokens from broken `module avail` output; surface probe errors; prefer GPU partitions via `prefer_partitions`
- **Equilibration (NAMD):** production (and any stage that finishes after the last `TIMING` print) no longer shows `0.0 ns/day` — use `Wall:` / final `WallClock:` for performance and elapsed time instead of overwriting wall time with `0.0` when promoting the final output step
- **Analysis (bilayer thickness):** when the membrane straddles the periodic z boundary, thickness no longer reports the water gap (`L_z − d` ≈ 100 Å); center the bilayer in z before lipyphilic and fold long PBC paths back to the headgroup–headgroup distance (~35-40 Å for POPC)
- **Analysis (bilayer time axis):** APL/thickness x-axis follows the sum of `file_times` instead of falling back to a tiny default dt
- **Equilibration (Amber/GROMACS):** if official ns/day is missing from the log, estimate it from wall elapsed × simulated time (completed stages included)
- **Equilibration (NVT):** Amber, NAMD, and GROMACS `01_NVT` templates are true constant-volume NVT through production (aligned with OpenMM); CHARMM-GUI’s NVT packs incorrectly matched NPT after early heating
- **Equilibration (NAMD NPAT):** steps 6.3–6.6 now use `useConstantArea` like CHARMM-GUI (they incorrectly had NPgT `useConstantRatio` / `SurfaceTensionTarget`)
- **Equilibration (GROMACS):** GPU `mdrun` lines now include `-ntmpi` (one rank per GPU) with `-ntomp`, required by GROMACS 2026 when OpenMP threads and GPUs are combined
- **Equilibration (GROMACS):** step0 energy minimisation runs on CPU only (`-ntmpi 1 -ntomp … -nb cpu -pme cpu`) — PME GPU does not support `steep`/`cg`, and GPU-capable `gmx` still requires `-ntmpi` with `-ntomp` even without GPU offload flags; MD stages still use GPU flags when requested
- **Equilibration (GROMACS):** omit unused `POSRES_FC_*` macros from MDP `define` lines (fixes grompp warnings for WATER/OTHER when those restraints are not in the topology)
- **Equilibration (GROMACS):** parse GROMACS 2026 energy-minimisation start banner (`Started Steepest Descents` / `Conjugate Gradients`) so job cards show elapsed wall time for minimization
- **Equilibration (GROMACS):** Kill MD / SIGTERM often still writes `Performance` + `Finished mdrun` before `nsteps` is reached — no longer treat that as stage completion or force progress to 100%; keep actual steps/ns and mark the stage interrupted

## [1.0.52] - 2026-07-25

### Fixed

- **Equilibration (GROMACS):** `run_equilibration.sh` now applies UI CPU/GPU settings on every `mdrun` (`-ntomp`, `-nb gpu`, `-pme gpu`, `-gpu_id`)
- **Equilibration (OpenMM):** run scripts pass `--device` / `--threads`; CPU compute target sets `PLATFORM=CPU`; `openmm_run.py` honors those flags via platform properties
- **Equilibration resume:** refreshing `run_equilibration.sh` preserves CPU/GPU settings from `equilibration_resources.json`

### Notes

- **NAMD** already emitted `+p` / `+devices`; unchanged. Amber equilibration remains unimplemented.

## [1.0.51] - 2026-07-22

### Added

- **Equilibration job metadata:** `equilibration_job.json` write/infer (`equilibration_job_metadata.py`) — input directory, ensemble, and protocol for GUI job cards and **Use in form**
- **GROMACS minimization progress:** energy-minimization logs report step count and wall time (not ns/day); `converged_early` when EM finishes before `nsteps`

### Fixed

- **Equilibration resume (OpenMM):** stage completion and resume point use log-based progress (aligned with GUI); canonical stage order/labels (`Equilibration 1`…`Production`); interrupted mid-stage-2 resumes **Equilibration 2** (not stage 1)

## [1.0.50] - 2026-07-21

### Added

- **Equilibration:** stage-level continue — `get_equilibration_resume_point()` and `RESUME=1` support in generated `run_equilibration.sh` (NAMD, GROMACS, OpenMM); skip completed stages and restart the first incomplete stage from the beginning

### Fixed

- **OpenMM:** `parse_openmm_log()` uses stage-local **Progress (%)** so per-stage simulated ns is correct when stages chain via `-irst` (cumulative Step values no longer inflate ns across rows)

## [1.0.49] - 2026-07-21

### Added

- **Engines:** `list_md_engine_candidates()` includes a `variant` field (`CPU` / `CUDA` / `OpenCL` / …) for GROMACS (`gmx --version` GPU support), NAMD (install path), and OpenMM (best available platform)

### Changed

- **Docs:** clarify conda GROMACS CPU vs CUDA — prefer CPU by default; CUDA matchspec is fragile next to OpenMM/`cudatoolkit` (solver hang, not an EULA prompt); document gatewizard-gui bootstrap defaults and `GATEWIZARD_CONDA_GROMACS_CUDA`

## [1.0.48] - 2026-07-15

### Fixed

- **Engines:** detect NAMD version from `Info: NAMD x.y.z` banner (was skipped as noise) and fall back to version in the install path (e.g. `NAMD_3.0.2_...`)

## [1.0.47] - 2026-07-13

### Added

- **Conda:** recommend `gromacs` from conda-forge in `environment.yml` (GUI bootstrap installs CUDA-then-CPU on Linux)
- **Engines:** `list_md_engine_candidates()` discovers NAMD / GROMACS (+ GMXRC) / OpenMM installs for version pickers
- **GROMACS:** run scripts only `source` GMXRC when provided (no hardcoded `/usr/local/gromacs`); works with conda `gmx`

## [1.0.46] - 2026-07-11

### Added

- **Preparation:** `strip_protein_hydrogens()` / `count_protein_hydrogens()` — remove protein H only (ligands/hetero kept)
- **Builder:** optional `remove_protein_h`; warns when protein hydrogens are detected
- **GROMACS equilibration:** positional restraints for `water`, `ions`, `other`, and custom MDAnalysis selections (NAMD/OpenMM parity); MDP macros `POSRES_FC_WATER` / `ION` / `OTHER`
- **Preparation:** `detect_terminal_caps()` / `is_already_capped()` to detect ACE/NME/NMA before re-capping

### Fixed

- **Preparation:** protein capping no longer drops ligands, waters, or ions that share a chain/segment with the protein
- **Preparation:** hydrogen stripping works when PDB topology lacks an `elements` attribute (guess elements / name fallback)
- **GROMACS equilibration:** POSRES use local atom indices for multi-copy ParmEd `system1` proteins; lipid includes target `system2` (not ions)
- **Builder:** stop premature post-process after async launch (false “No bilayer PDB” / “File conversion failed” warnings)

## [1.0.45] - 2026-07-10

### Added

- **Packmol hydration:** cavity fill via PACKMOL (volume estimate, job prep, run, custom input)
- Bundled TIP3P water template, tests, examples, and API docs

### Fixed

- PACKMOL runs via stdin on AmberTools Memgen (fixes input parsing errors)
- Job-relative paths, fixed solute in input, corrected TIP3P geometry
- Volume estimate returns free grid points for GUI ghost-water preview

### Changed

- Builder default `notprotonate=True` so packmol-memgen keeps PropKa residue names (GLH/ASH/…)

## [1.0.44] - 2026-07-03

### Added

- Builder passes packmol-memgen `--dist` with GateWizard default `12 Å`

### Changed

- Builder default water layer thickness (`dist_wat`) is now `26 Å`

## [1.0.43] - 2026-07-01

### Removed

- The `gatewizard` console command — this package is a Python library only; use [gatewizard-gui](https://github.com/franciscoadasme/gatewizard-gui) for the desktop app

## [1.0.42] - 2026-06-26

### Fixed

- PSIQUE SS results are re-keyed to structure chain IDs (e.g. `PROT` segid → `A`) so assignments are not silently dropped as coil
- Fall through to PDB records / CA-angle heuristic when PSIQUE returns no helix or strand coverage
- Protein-only PSIQUE input PDB now uses normalized chain IDs from our structure writer

## [1.0.41] - 2026-06-26

### Fixed

- PSIQUE secondary structure assignment retries on a protein-only temp PDB when the full structure fails (e.g. CHARMM/NAMD EPW extra-point waters with unknown elements)
- Shared `resolve_pdb_chain_id()` for CHARMM-style long segids (`PROT`, `MEMB`, …) so SS keys use the PDB chain letter instead of the segment name

## [1.0.40] - 2026-06-19

### Added

- `compute_orientation_transform()` and `apply_orientation_transform()` to apply MemPro orientations to a full structure via rigid-body fit
- `StructureManager.apply_mempro_orientation()` for the same workflow on a loaded structure

### Fixed

- Secondary structure assignment still tries PSIQUE first; if PSIQUE fails or returns no results (e.g. on Windows), the API falls through to PDB HELIX/SHEET records and then the CA-angle heuristic instead of leaving SS unassigned
- Export `assign_secondary_structure_map()` for callers that need per-residue SS without StructureManager
- Preparation shell script exports `AMBERHOME` and `CONDA_PREFIX/bin` on `PATH` for tleap on macOS
- Embed absolute `loadPDB` path in `leap_parametrize.in` at generation time (macOS `sed -i` incompatibility left `PREPARED_PDB_PLACEHOLDER` unreplaced)
- Run tleap via `subprocess_argv_for_script()` with `stdin=DEVNULL` to avoid non-interactive `tl_getline` failures
- `resolve_conda_executable()` and improved `subprocess_argv_for_script()` run conda env scripts with the active env's Python, fixing PropKa and pdb4amber after the GUI runtime folder move left stale shebang paths

## [1.0.39] - 2026-06-19

### Fixed

- `subprocess_argv_for_script()` runs conda wrapper scripts via explicit interpreter argv, fixing `pdb4amber` ENOENT on macOS when the env path contains spaces (e.g. `Application Support`)
- `run_pdb4amber_with_cap_fix()` and builder `pdb4amber` step use the new subprocess argv helper

## [1.0.37] - 2026-06-19

### Fixed

- `run_pdb4amber_with_cap_fix()` resolves `pdb4amber` from `CONDA_PREFIX/bin` before subprocess execution, fixing Prepare failures in packaged desktop apps (e.g. GateWizard GUI on macOS) where `PATH` is minimal and the bare `pdb4amber` command was not found

## [1.0.36] - 2026-06-08

### Fixed

- PyPI publish blocked by MemPrO VCS dependency in optional extras — removed `mempro @ git+...` from `pyproject.toml` metadata (PyPI forbids direct URL dependencies in published wheels)

### Changed

- MemPrO install moved to `requirements-orientation.txt`, pinned to GitHub release `v0.1.0`
- Updated install docs and dependency messages for `gatewizard[full]` + orientation workflow

## [1.0.35] - 2026-06-07

### Changed

- Equilibration stage 4 (step 6.4) now uses a **1 fs** integration timestep across NAMD, GROMACS, and OpenMM default protocols and templates, instead of the CHARMM-GUI default 2 fs jump at this stage. Stages 5–6 and production remain at 2 fs. This improves stability for membrane systems still equilibrating box dimensions under NPT/NPAT/NPγT.

## [1.0.34] - 2026-05-01

### Fixed

- `run_structural_analysis()` RMSF `residue_number` mode now returns numeric residue IDs as x-values instead of mixed `"{resname}{resid}"` strings, so downstream plots and data consumers receive a consistent numeric series for this mode

## [1.0.33] - 2026-05-01

### Fixed

- Fix TrajectoryAnalyzer to use instance attributes (self.topology, self.trajectories) instead of temporary locals when creating MDAnalysis Universe and when loading individual trajectory files to count frames, preventing incorrect file/path usage when handling multiple trajectories.

## [1.0.32] - 2026-05-01

### Added

- `rmsf_xaxis_type` parameter in `run_structural_analysis()` — controls RMSF x-axis labeling with three modes: `"residue_number"` (default), `"residue_type_number"` (e.g. `ALA42`), and `"atom_index"`

### Fixed

- `TrajectoryAnalyzer` used stale `self.trajectories` / `self.topology` references instead of the resolved `tmp_trajs` / `tmp_top` variables when loading trajectories, causing incorrect universe construction

## [1.0.31] - 2026-04-30

### Fixed

- `NameError: name 'Any' is not defined` in `namd_analysis.py` — added `Any` to `typing` imports

## [1.0.30] - 2026-04-30

### Added

- API-oriented analysis helpers in `gatewizard.utils.namd_analysis`:
  - `run_structural_analysis()` for RMSD, RMSF, distance, and radius of gyration
  - `list_namd_energy_properties()` for NAMD log property discovery
  - `run_energetic_analysis()` returning JSON-serializable multi-series data


## [1.0.29] - 2026-04-27

### Added

- `parametrize_ligand_from_system_pdb()` to centralize extraction plus parametrization

## [1.0.28] - 2026-04-27

### Changed

- Removed `mempro` from PyPI dependencies (PyPI does not allow direct git URL dependencies); must be installed manually with `pip install git+https://github.com/pstansfeld/MemPrO.git`
- Added manual MemPrO installation instructions to README under External Requirements

## [1.0.27] - 2026-04-27

### Added

- Beta warning message displayed in the title bar next to the app name, indicating the app is not ready for production use

### Added

- **MemPrO integration** for membrane protein orientation in the Visualize tab
  - New collapsible "MemPrO Orientation" section with Run/Load Orient Folder buttons
  - GUI controls for key options: dual membrane, peripheral, B-factor weighting, flip, CPUs, membrane thickness, grid size, iterations
  - Ranked results list with potential/hits scores; click any rank to load the oriented PDB
  - Background execution with orange status bar indicator
  - New `gatewizard.core.mempro` API module (`MemPrO`, `OrientationResult`, `MemProError`)
  - API documentation, test suite, and 12 example scripts

### Changed

- Use official PSIQUE Python module instead of bundled binary and wrappers.

### Fixed

- **Propka crash** on proteins with OXT atoms (`ValueError: list.remove(x): x not in list`): OXT atoms are now stripped from a temporary PDB before running propka
- Keep HET atoms when adding caps (#10)

## [1.0.26] - 2026-03-29

### Added

- **Initial/Final 2D structure toggle** on each ligand card: view the PDB-derived structure (Initial) or the parametrized mol2 structure with correct bond orders (Final)
- Auto-switch to Final view after successful parametrization or when cached results are detected
- GAFF→SYBYL atom type conversion for loading antechamber mol2 files into RDKit (`_GAFF_TO_SYBYL` mapping, `_load_gaff_mol2()`, `_mol2_has_gaff_types()`)
- SS Assignment segmented button (PSIQUE / PDB / Heuristic) in Visualize controls to manually select the secondary structure method
- **MDAnalysis-based atom selections for equilibration restraints:**
  - `DEFAULT_SELECTIONS` class attribute with MDAnalysis selection strings for all 7 standard restraint categories
  - `count_selection_atoms()` — count atoms matching any MDAnalysis selection expression
  - `get_default_selections()` — build selection dict with auto-detected `ligand_<RESNAME>` entries for non-standard residues
  - `count_all_selections()` — count atoms for all selections in one call
  - `generate_restraints_file_mda()` — generate restraint PDB using MDAnalysis selections instead of heuristic
  - `generate_restraints_file()` now accepts optional `selections` parameter to use MDAnalysis mode
  - GUI: atom count labels, gear button (⚙) to edit selections, Add Selection (+) button per stage
  - GUI: auto-detect ligands when input folder is selected
  - New equilibration example 08 demonstrating MDAnalysis selection features

### Changed

- Default SS assignment order is now PSIQUE → PDB records → heuristic (was PDB → PSIQUE → heuristic)
- SearchableComboBox dropdown is now a floating overlay instead of an inline child (no longer expands the parent section)
- SearchableComboBox arrows changed from filled triangles to V-shape chevrons matching native CTk style
- Removed complex `_infer_bond_orders()` logic; 2D ligand images now use simple PDB loading for Initial view

### Fixed

- Removed DrawEngine monkey-patch; `circle_shapes` drawing method already renders V-shape arrows natively
- Replaced Unicode sort arrows (▲/▼) in Preparation frame with PIL-drawn images
- GAFF atom types (`ca`, `c3`, `os`, etc.) in antechamber mol2 files no longer cause RDKit `Element not found` errors

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
- Box dimensions now correctly read from bilayer\_\*\_lipid.pdb file for NAMD equilibration configurations

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
