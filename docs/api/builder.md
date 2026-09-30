# Builder Module

Module for preparing membrane protein systems with automated lipid bilayer construction, solvation, and parametrization. This module includes:

- Membrane system building with packmol-memgen
- Custom lipid composition (asymmetric membranes supported)
- Force field selection and validation
- Salt addition and ionic strength control
- Two-stage Propka workflow support
- Automated parametrization with pdb4amber and tleap

## Import

```python
from gatewizard.core.builder import Builder
from gatewizard.tools.force_fields import ForceFieldManager
```

## Class: Builder

Main class for building membrane protein systems with complete control over lipid composition, force fields, and system parameters.

### Constructor

```python
Builder()
```

**Parameters:** None

**Returns:** `Builder` instance

**Default Configuration:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `water_model` | `"tip3p"` | Water model for solvation |
| `protein_ff` | `"ff19SB"` | Protein force field |
| `lipid_ff` | `"lipid21"` | Lipid force field |
| `preoriented` | `True` | Protein already oriented for membrane insertion |
| `parametrize` | `True` | Run parametrization with tleap |
| `salt_concentration` | `0.15` | Salt concentration in M (molar) |
| `cation` | `"K+"` | Cation type for salt |
| `anion` | `"Cl-"` | Anion type for salt |
| `dist` | `12` | Minimum solute-to-box-boundary distance in Å |
| `dist_wat` | `26` | Water layer thickness in Å |
| `notprotonate` | `True` | Skip protonation (preserve PropKa residue names) |

### Example 1: Basic Configuration
```python
--8<-- "tests/builder_examples/builder_example_01.py"
```

---

## Configuration Methods

### set_configuration

Update configuration parameters for system building.

```python
set_configuration(**kwargs)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `water_model` | `str` | `"tip3p"` | Water model (tip3p, tip4p, spce, opc) |
| `protein_ff` | `str` | `"ff19SB"` | Protein force field (ff14SB, ff19SB) |
| `lipid_ff` | `str` | `"lipid21"` | Lipid force field (lipid17, lipid21) |
| `preoriented` | `bool` | `True` | Whether protein is pre-oriented |
| `parametrize` | `bool` | `True` | Run parametrization after packing |
| `salt_concentration` | `float` | `0.15` | Salt concentration in M |
| `cation` | `str` | `"K+"` | Cation type (Na+, K+, etc.) |
| `anion` | `str` | `"Cl-"` | Anion type (Cl-, Br-, etc.) |
| `dist` | `float` | `12` | Minimum solute-to-box-boundary distance in Å |
| `dist_wat` | `float` | `26` | Water layer thickness in Å |
| `notprotonate` | `bool` | `True` | Skip protonation during parametrization |
| `add_salt` | `bool` | `True` | Whether to add salt to system |

**Returns:** None

### Example 2: Custom Configuration
```python
--8<-- "tests/builder_examples/builder_example_02.py"
```

### Example 3: Available Water Models
```python
--8<-- "tests/builder_examples/builder_example_03.py"
```

Similarly can be done for lipids:

### Example 4: Available Lipid Models
```python
--8<-- "tests/builder_examples/builder_example_04.py"
```


### Example 5: Available Protein Force Fields
```python
--8<-- "tests/builder_examples/builder_example_05.py"
```

### Example 6: Available Lipid Force Fields
```python
--8<-- "tests/builder_examples/builder_example_06.py"
```

---

## Validation Methods

Before preparing systems, you can validate your inputs and check force field compatibility.

### validate_system_inputs

Validate all inputs before system preparation.

```python
validate_system_inputs(
    pdb_file: str,
    upper_lipids: List[str],
    lower_lipids: List[str],
    lipid_ratios: str = "",
    **kwargs
) -> Tuple[bool, str]
```

**Parameters:**

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `pdb_file` | `str` | No | Path to input PDB file, or empty/None for a bilayer-only build |
| `upper_lipids` | `List[str]` | Yes | List of lipids for upper leaflet |
| `lower_lipids` | `List[str]` | Yes | List of lipids for lower leaflet |
| `lipid_ratios` | `str` | No | Lipid ratios string |
| `**kwargs` | `Any` | No | Additional parameters to validate |

**Returns:** `Tuple[bool, str]`

- `bool`: Validation status
- `str`: Error message (empty if valid)

**Validation Checks:**

1. PDB file exists and is readable
2. All lipid types are valid/available
3. Lipid ratios match number of lipids
4. Ratios are positive numbers
5. Force field combinations are compatible
6. Ion types are valid

### Example 7: Input Validation
```python
--8<-- "tests/builder_examples/builder_example_07.py"
```

### Example 8: Force Field Validation
```python
--8<-- "tests/builder_examples/builder_example_08.py"
```

---

## Core Preparation Methods

### prepare_system

Prepare a complete membrane protein system with lipid bilayer, water, and ions.

```python
prepare_system(
    pdb_file: str,
    working_dir: str,
    upper_lipids: List[str],
    lower_lipids: List[str],
    lipid_ratios: str = "",
    **kwargs
) -> Tuple[bool, str, Optional[Path]]
```

**Parameters:**

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `pdb_file` | `str` | No | Path to input PDB file (oriented protein). Omit or pass `None` for a bilayer-only build |
| `working_dir` | `str` | Yes | Working directory for output files |
| `upper_lipids` | `List[str]` | Yes | List of lipid types for upper leaflet |
| `lower_lipids` | `List[str]` | Yes | List of lipid types for lower leaflet |
| `lipid_ratios` | `str` | No | Lipid molar ratios (format: "ratio1:ratio2//ratio3:ratio4") |
| `output_folder_name` | `str` | No | Custom output folder name (kwarg) |
| `distxy_fix` | `float` | No | Membrane XY size in Å (`--distxy_fix`). **Required** when `pdb_file` is omitted |
| `solutes` | `list` | No | Free molecules to pack: `[{"pdb": "TEA.pdb", "concentration": "4", "in_membrane": False}]` |
| `solute_inmem` | `bool` | No | Place free molecules in the membrane instead of water (`--solute_inmem`) |
| `solute_prot_dist` | `float` | No | Keep free molecules this far from the protein, in Å (`--solute_prot_dist`) |
| `wait` | `bool` | No | Block until the job completes or errors (default `False`) |
| `wait_timeout` | `float` | No | Maximum seconds to wait when `wait=True` (`None` = unlimited) |
| `wait_poll_interval` | `float` | No | Seconds between status checks (default `5.0`) |
| `wait_verbose` | `bool` | No | Print elapsed-time progress while waiting (default `True`) |
| `**kwargs` | `Any` | No | Additional configuration parameters (override defaults) |

**Returns:** `Tuple[bool, str, Optional[Path]]`

- `bool`: Success status
- `str`: Status message
- `Optional[Path]`: Job directory path (None if failed)

**Lipid Ratios Format:**

The `lipid_ratios` string specifies molar ratios for each leaflet:

- Format: `"upper_ratio1:upper_ratio2//lower_ratio1:lower_ratio2"`
- Ratios must match the order of lipids in the lipid lists
- Ratios are normalized automatically (don't need to sum to 1.0)
- Use `//` to separate upper and lower leaflet ratios

**Available Lipids:**

| Category | Lipid Types |
|----------|------------|
| **Phosphatidylcholine (PC)** | DHPC, DLPC, DMPC, DPPC, DSPC, DOPC, POPC, PLPC, SOPC |
| **Phosphatidylethanolamine (PE)** | DHPE, DLPE, DMPE, DPPE, DSPE, DOPE, POPE, PLPE, SOPE |
| **Phosphatidylserine (PS)** | DHPS, DLPS, DMPS, DPPS, DSPS, DOPS, POPS, PLPS, SOPS |
| **Phosphatidylglycerol (PG)** | DHPG, DLPG, DMPG, DPPG, DSPG, DOPG, POPG, PLPG, SOPG |
| **Phosphatidic Acid (PA)** | DHPA, DLPA, DMPA, DPPA, DSPA, DOPA, POPA, PLPA, SOPA |
| **Sphingomyelins** | PSM, ASM, LSM, MSM, HSM, OSM |
| **Cholesterol** | CHL1 |
| **Other** | CARDIOLIPIN, PIP, PIP2 |

**Output Files:**

In `{output_folder_name}/` directory:

- `bilayer_*.pdb` - Packed system (protein in membrane)
- `bilayer_*_lipid.pdb` - Lipid bilayer with CRYST1 info (for equilibration)
- `system.prmtop` - AMBER topology file (if parametrize=True)
- `system.inpcrd` - AMBER coordinate file (if parametrize=True)
- `system_solv.pdb` - Solvated system in PDB format
- `system_4amber.pdb` - PDB prepared by pdb4amber
- `tleap.in` - tleap input script
- `logs/preparation.log` - Main execution log
- `logs/parametrization.log` - Parametrization log
- `status.json` - Job status tracking
- `run_preparation.sh` - Execution script

**Raises:**

- `FileNotFoundError` - If input PDB file doesn't exist
- `BuilderError` - If validation fails or system preparation fails

**Note:** By default, system preparation runs in the background and the method returns immediately.
Pass `wait=True` to block until the job finishes — useful for scripting sequential preparations.
For asynchronous monitoring, use the log files or the `JobMonitor` class (see examples below).

### Bilayer-only and free molecules

packmol-memgen can pack a membrane **without a protein**. Omit `pdb_file` and set `distxy_fix` (Å) for the XY size. Free copies of a molecule go in `solutes` (`--solute` / `--solute_con`); parametrize them first and pass `ligand_params` plus GAFF2. The membrane-protein path (`pdb_file` required, `--pdb` + `--preoriented`) is unchanged.

```python
from gatewizard.core.builder import Builder

builder = Builder()

# Membrane only (no protein): XY size is required
success, msg, job_dir = builder.generate_preparation_inputs(
    pdb_file=None,
    working_dir="./systems",
    upper_lipids=["DOPE", "DOPG"],
    lower_lipids=["DOPE", "DOPG"],
    lipid_ratios="3:1//3:1",
    distxy_fix=100,
)

# Free ligand in water around a protein
success, msg, job_dir = builder.generate_preparation_inputs(
    pdb_file="1BL8.pdb",
    working_dir="./systems",
    upper_lipids=["DOPE", "DOPG"],
    lower_lipids=["DOPE", "DOPG"],
    lipid_ratios="3:1//3:1",
    solutes=[{"pdb": "TEA.pdb", "concentration": "4"}],
    solute_prot_dist=10,
    ligand_params={"TEA": {"frcmod": "TEA.frcmod", "lib": "TEA.lib"}},
)

# Ligand + bilayer, no protein (place copies in the membrane)
success, msg, job_dir = builder.generate_preparation_inputs(
    pdb_file=None,
    working_dir="./systems",
    upper_lipids=["POPC"],
    lower_lipids=["POPC"],
    lipid_ratios="1//1",
    distxy_fix=75,
    solutes=[{"pdb": "TEA.pdb", "concentration": "0.1M", "in_membrane": True}],
    ligand_params={"TEA": {"frcmod": "TEA.frcmod", "lib": "TEA.lib"}},
)
```

`concentration` is a molecule count (`"4"`) or a packmol-memgen concentration (`"0.1M"`, `"2%"`). `--solute_inmem` is set when any solute has `in_membrane=True` or when `solute_inmem=True`. `--solute_prot_dist` is only passed when a protein PDB is present.

### Example 9: Simple Symmetric Membrane
```python
--8<-- "tests/builder_examples/builder_example_09.py"
```

### Example 10: Monitoring Job Progress
```python
--8<-- "tests/builder_examples/builder_example_10.py"
```

### Example 11: Asymmetric Membrane with Multiple Lipids
```python
--8<-- "tests/builder_examples/builder_example_11.py"
```

### Example 12: Complex Composition (Plasma Membrane Mimic)
```python
--8<-- "tests/builder_examples/builder_example_12.py"
```

### Example 13: Packing Only (No Parametrization)
```python
--8<-- "tests/builder_examples/builder_example_13.py"
```

### Example 14: Custom Salt Concentration
```python
--8<-- "tests/builder_examples/builder_example_14.py"
```

### Example 15: No Salt (Charge Neutralization Only)
```python
--8<-- "tests/builder_examples/builder_example_15.py"
```

---

### run_preparation

```python
success, message = builder.run_preparation(job_dir)
```

Launch a previously generated preparation job (`run_preparation.sh` in `job_dir`). Returns `(False, ...)` if the script is missing or a job is already running.

**Parameters:**
- `job_dir` (str or Path): Job directory created by `prepare_system` / `generate_preparation_inputs`.

### cancel_preparation

```python
result = builder.cancel_preparation(job_dir) -> dict
```

Cancel a running preparation job. Kills the detached `run_preparation.sh` process group (PID in `process.pid`) and marks `status.json` as `cancelled`. If the job already finished, returns `stopped=False` and the existing status.

**Returns:** `{success, job_dir, stopped, status, message}`.

### prepare_system_stage1_for_propka

```python
success, message, job_dir = builder.prepare_system_stage1_for_propka(
    pdb_file, working_dir, upper_lipids, lower_lipids, lipid_ratios="", **kwargs
)
```

Stage 1 of the two-stage PropKa workflow: pack only (`parametrize=False`). After packing, apply PropKa residue names, then call stage 2.

### prepare_system_stage2_for_propka

```python
success, message, job_dir = builder.prepare_system_stage2_for_propka(
    packed_pdb_file, working_dir, **kwargs
)
```

Stage 2: parametrize a packed PDB that already has PropKa residue names (`notprotonate=True`).

### modify_residue_names_for_propka

```python
ok, message = builder.modify_residue_names_for_propka(
    packed_pdb_file, propka_results, output_file=None
)
```

Rewrite residue names in a packed PDB from a `{resid: new_name}` (or similar) PropKa map before stage 2.

See [Example 25](#example-25-blocking-synchronous-preparation-with-waittrue) for `wait=True` on `prepare_system`.

### wait_for_completion

Block until a preparation job completes or fails by polling `status.json`.

```python
wait_for_completion(
    job_dir,
    poll_interval: float = 5.0,
    timeout: Optional[float] = None,
    verbose: bool = True,
) -> Tuple[bool, str]
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `job_dir` | `Path \| str` | — | Path to the job directory |
| `poll_interval` | `float` | `5.0` | Seconds between status checks |
| `timeout` | `float \| None` | `None` | Maximum seconds to wait (`None` = no limit) |
| `verbose` | `bool` | `True` | Print elapsed-time progress to stdout |

**Returns:** `Tuple[bool, str]` — (success, message)

### Example 25: Blocking (synchronous) Preparation with wait=True
```python
--8<-- "tests/builder_examples/builder_example_25.py"
```

---

## Class: ForceFieldManager

```python
from gatewizard.tools.force_fields import ForceFieldManager

ff = ForceFieldManager()
```

### get_water_models

```python
ff.get_water_models() -> list[str]
```

See [Example 3](#example-3-available-water-models).

### get_protein_force_fields

See [Example 5](#example-5-available-protein-force-fields).

### get_lipid_force_fields

See [Example 6](#example-6-available-lipid-force-fields).

### get_available_lipids

Lipid residue names accepted by packmol-memgen / Lipid21.

### get_available_cations

### get_available_anions

### validate_combination

```python
valid, message, is_warning = ff.validate_combination(water, protein_ff, lipid_ff)
```

See [Example 8](#example-8-force-field-validation).

### get_recommendations

```python
ff.get_recommendations(system_type="membrane") -> dict
```

### get_force_field_info

Return metadata for one force-field name.

### validate_lipid

```python
ff.validate_lipid(lipid_name) -> bool
```

### validate_ion

```python
ok, charge = ff.validate_ion(ion_name)
```

## Job Monitoring

Since system preparation runs in the background, GateWizard provides the `JobMonitor` class to track progress, check status, and manage running jobs.

### Class: JobMonitor

Monitor and track system preparation jobs.

### Example 16: JobMonitor class
```python
--8<-- "tests/builder_examples/builder_example_16.py"
```

**Constructor Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `working_directory` | `Path` or `str` | `Path.cwd()` | Directory to monitor for jobs |

**Returns:** `JobMonitor` instance

---

### scan_for_jobs

Scan the working directory for preparation jobs.

```python
monitor.scan_for_jobs(force=False)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `force` | `bool` | `False` | Force scan even if within scan interval (3 seconds) |

**Returns:** None (updates internal job list)

**Description:** Scans the working directory for subdirectories containing `status.json` files, which indicate active or completed preparation jobs.

---

### get_active_jobs

Get all currently running jobs.

```python
active_jobs = monitor.get_active_jobs()
```

**Returns:** `Dict[str, JobInfo]`

- Dictionary mapping job IDs to `JobInfo` objects
- Only includes jobs with status `RUNNING`

---

### get_completed_jobs

Get all completed or errored jobs.

```python
completed_jobs = monitor.get_completed_jobs()
```

**Returns:** `Dict[str, JobInfo]`

- Dictionary mapping job IDs to `JobInfo` objects
- Includes jobs with status `COMPLETED` or `ERROR`

---

### get_job

Get specific job by ID.

```python
job = monitor.get_job(job_id)
```

**Parameters:**

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `job_id` | `str` | Yes | Job identifier (usually the full job directory path) |

**Returns:** `Optional[JobInfo]`

- `JobInfo` object if job exists
- `None` if job not found

---

### remove_job

```python
removed = monitor.remove_job(job_id) -> bool
```

Drop a job from the in-memory monitor. Does **not** delete files on disk.

### cleanup_stale_jobs

```python
monitor.cleanup_stale_jobs(max_age_seconds: float = 3600)
```

Remove jobs whose last update is older than `max_age_seconds`.

### get_job_statistics

```python
stats = monitor.get_job_statistics() -> dict
```

Counts of `total`, `running`, `completed`, `error`, and `unknown` jobs.

See [Example 16](#example-16-jobmonitor-class) and [Example 18](#example-18-monitoring-batch-job-management).

### refresh_job

Force refresh of a specific job's status.

```python
success = monitor.refresh_job(job_id)
```

**Parameters:**

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `job_id` | `str` | Yes | Job identifier to refresh |

**Returns:** `bool`

- `True` if job was refreshed successfully
- `False` if job not found or refresh failed

---

### JobInfo Object

Each job is represented by a `JobInfo` object with the following attributes:

| Attribute | Type | Description |
|-----------|------|-------------|
| `job_id` | `str` | Unique job identifier |
| `job_dir` | `Path` | Job directory path |
| `status` | `JobStatus` | Current status (RUNNING, COMPLETED, ERROR, UNKNOWN) |
| `progress` | `float` | Progress percentage (0.0 - 100.0) |
| `current_step` | `str` | Name of current processing step |
| `steps_completed` | `List[str]` | List of completed step names |
| `elapsed_time` | `float` | Elapsed time in seconds |
| `start_time` | `datetime` | Job start timestamp |
| `end_time` | `Optional[datetime]` | Job end timestamp (None if running) |

**JobStatus Values:**

- `RUNNING` - Job is currently executing
- `COMPLETED` - Job finished successfully
- `ERROR` - Job encountered an error
- `UNKNOWN` - Status cannot be determined

---

### Example 17: Real-time Progress Tracking

```python
--8<-- "tests/builder_examples/builder_example_17.py"
```

---

### Example 18: Monitoring Batch Job Management

```python
--8<-- "tests/builder_examples/builder_example_18.py"
```

---

### Monitoring Tips

**Real-time monitoring:**

- Use `scan_for_jobs(force=True)` to bypass the 3-second scan interval
- Check status every 1-2 seconds for responsive progress updates
- Use `flush=True` in print statements for immediate display

**Long-running jobs:**

- Check status periodically (every 30-60 seconds)
- Log progress to file for later review
- Use `elapsed_time` to estimate remaining time

**Multiple jobs:**

- Monitor entire working directory for batch processing
- Filter jobs by status for targeted monitoring
- Use job_dir.name for user-friendly display

**Debugging:**

- Check `job_dir / 'logs/preparation.log'` for detailed output
- Check `job_dir / 'logs/parametrization.log'` if parametrization fails
- Examine `status.json` for detailed job state

**Manual log inspection:**
```bash
# Follow preparation log in real-time
tail -f systems/popc_membrane/logs/preparation.log

# Check for errors
grep -i error systems/popc_membrane/logs/*.log

# View job status
cat systems/popc_membrane/status.json
```

---

## Tips and Best Practices

### Lipid Selection

**For simple systems:**

- Use 100% POPC for both leaflets (symmetric membrane)
- Well-characterized, stable, and widely used
- Good starting point for method development

**For realistic membranes:**

- Include cholesterol (10-30%) for membrane stability
- Use PE lipids (POPE) for inner leaflet (more physiological)
- Add PS lipids (POPS) for negative charge (cytoplasmic side)

**For specialized studies:**

- Use asymmetric compositions (different upper/lower)
- Match experimental lipid compositions when available
- Consider protein's native membrane environment

### Force Field Selection

**Default recommendation (most stable):**

- Water: tip3p
- Protein: ff14SB
- Lipid: lipid21

**For latest parameters:**

- Water: opc
- Protein: ff19SB
- Lipid: lipid21

### Salt Concentration

**Physiological (default):**

- 0.15 M (150 mM) NaCl or KCl
- Mimics intracellular or extracellular conditions

**High salt studies:**

- 0.5-1.0 M for ionic strength effects
- May require longer equilibration

**Charge neutralization only:**

- Set `salt_concentration=0.0`
- System will still be neutralized with minimum ions

### Water Layer Thickness

**Default (26 Å):**

- Sufficient for most membrane proteins
- ~3-4 water layers above/below membrane

**Boundary distance (12 Å default):**

- Set `dist` to control the minimum distance from the solute extents to box boundaries
- packmol-memgen may choose a larger real distance because it uses the worst-case scenario

**Large proteins or complexes:**

- Increase to 20-25 Å
- Ensures adequate solvation

**Minimize system size:**

- Decrease to 15 Å (minimum recommended)
- Faster simulations but ensure full solvation

### Protein Orientation

**Pre-oriented proteins:**

- Use [OPM database](https://opm.phar.umich.edu/)
- Protein already positioned in membrane plane
- Set `preoriented=True` (default)

**Non-oriented proteins:**

- Let packmol-memgen use MEMEMBED for orientation
- Set `preoriented=False`
- Adds extra time to preparation

### Troubleshooting

**"PDB file not found":**

- Check absolute/relative paths
- Verify file permissions
- Ensure PDB file is valid format

**"Packing failed":**

- Protein may be too large for default settings
- Try increasing `dist_wat` parameter
- Check protein orientation (OPM database)
- Verify protein is pre-oriented if using `preoriented=True`

**"Parametrization failed":**

- Check logs/parametrization.log for details
- May have non-standard residues
- Verify protonation states (use Propka workflow)
- Check for missing atoms or unusual residue names

**"Invalid lipid ratios":**

- Ensure ratios match number of lipids
- Format: "ratio1:ratio2//ratio3:ratio4"
- Ratios don't need to sum to 1.0 (auto-normalized)
- Example: "7:3//5:5" for 2 lipids per leaflet

**"Force field compatibility error":**

- Use ForceFieldManager to check combinations
- Not all force fields are compatible
- Use recommended combinations for stability

### Performance Tips

**Faster preparation:**

- Use symmetric membranes (same upper/lower)
- Fewer lipid types (1-2 per leaflet)
- Pre-oriented proteins
- Smaller water layers (15-20 Å)

**More realistic systems:**

- Asymmetric membranes (different leaflets)
- Multiple lipid types (3-4 per leaflet)
- Include cholesterol
- Larger water layers (20-25 Å)
- Use Propka workflow for correct protonation

---

## Ligand Parametrization

Module for detecting, extracting, and parametrizing non-standard (ligand) residues found in PDB files.
Uses the AMBER/GAFF workflow (antechamber → parmchk2 → tleap) to generate force field parameters
that integrate with the Builder's packmol-memgen and final tleap steps.
Supports both GAFF and GAFF2 atom type sets, with recommended pairings
per the AMBER manual: **gaff/bcc** and **gaff2/abcg2**.

### Import

```python
from gatewizard.tools.ligand_parametrization import (
    detect_ligands,
    extract_ligand_pdb,
    parametrize_ligand,
    parametrize_all_ligands,
    get_ligand_2d_image,
    get_ligand_2d_image_from_pdb_lines,
    build_ligand_param_args,
    build_tleap_ligand_lines,
    LigandInfo,
    LigandParametrizationError,
    STANDARD_RESIDUES,
    CHARGE_METHODS,
    ATOM_TYPES,
    DEFAULT_CHARGE_METHOD,
    DEFAULT_ATOM_TYPE,
    RECOMMENDED_COMBOS,
    NON_RECOMMENDED_COMBOS,
)
```

### Class: LigandInfo

Information container for a detected ligand residue.

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | `str` | 3-letter residue name (e.g., `"AAA"`) |
| `chain` | `str` | Chain identifier |
| `res_id` | `int` | Residue sequence number |
| `num_atoms` | `int` | Number of atoms in the ligand |
| `elements` | `Dict[str, int]` | Element counts (e.g., `{'C': 9, 'H': 10, 'O': 1}`) |
| `pdb_lines` | `List[str]` | Raw HETATM lines from the PDB file |
| `formula` | `str` (property) | Molecular formula string (e.g., `"C9H10O"`) |

**Methods:**

| Method | Returns | Description |
|--------|---------|-------------|
| `to_dict()` | `Dict[str, Any]` | Convert to serializable dictionary |

### Constants

| Constant | Type | Description |
|----------|------|-------------|
| `STANDARD_RESIDUES` | `set` | Residue names excluded from detection (amino acids, water, ions, lipids, capping groups) |
| `CHARGE_METHODS` | `dict` | Antechamber charge methods: `'bcc'` (AM1-BCC), `'abcg2'` (ABCG2), `'gas'` (Gasteiger), `'mul'` (Mulliken), `'cm2'` (CM2), `'rc'` (Read-in), `'resp'` (RESP — requires Gaussian), `'esp'` (ESP — requires Gaussian) |
| `ATOM_TYPES` | `dict` | Antechamber atom type sets: `'gaff2'` → `'GAFF2'`, `'gaff'` → `'GAFF'` |
| `DEFAULT_CHARGE_METHOD` | `str` | `'bcc'` |
| `DEFAULT_ATOM_TYPE` | `str` | `'gaff2'` |
| `RECOMMENDED_COMBOS` | `set` | Recommended (atom_type, charge_method) pairings per AMBER manual: `{('gaff', 'bcc'), ('gaff2', 'abcg2')}` |
| `NON_RECOMMENDED_COMBOS` | `set` | Non-recommended pairings (warning shown in GUI): `{('gaff2', 'bcc'), ('gaff', 'abcg2')}` |

!!! warning "Recommended Atom Type / Charge Method Pairings"
    Per the AMBER manual, the efficient charge models `bcc` (AM1-BCC) and `abcg2` (ABCG2)
    should be paired with specific atom type sets:

    | Atom Type | Charge Method | Status |
    |-----------|---------------|--------|
    | `gaff` | `bcc` | **Recommended** ✓ |
    | `gaff2` | `abcg2` | **Recommended** ✓ |
    | `gaff2` | `bcc` | Not recommended ✗ |
    | `gaff` | `abcg2` | Not recommended ✗ |

    The GUI shows an amber warning label when a non-recommended combination is selected,
    and displays a confirmation dialog before proceeding with parametrization.

!!! info "External QM Software Requirements"
    The `resp` and `esp` charge methods require **Gaussian** (external quantum-mechanics
    software, not included in AmberTools). All other methods (`bcc`, `abcg2`, `gas`,
    `mul`, `cm2`) use **sqm**, which is bundled with AmberTools.

### detect_ligands

Detect non-standard (ligand) residues in a PDB file by scanning HETATM records.

```python
detect_ligands(pdb_file: str) -> List[LigandInfo]
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `pdb_file` | `str` | Path to the PDB file to analyze |

**Returns:** `List[LigandInfo]` — one entry per unique ligand residue name

**Raises:** `LigandParametrizationError` if the file cannot be read

### Example 19: Detect Ligands

```python
--8<-- "tests/builder_examples/builder_example_19.py"
```

---

### extract_ligand_pdb

Extract a single ligand from a PDB file into its own file.

```python
extract_ligand_pdb(pdb_file: str, ligand_name: str, output_dir: str) -> str
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `pdb_file` | `str` | Path to the source PDB file |
| `ligand_name` | `str` | 3-letter residue name of the ligand |
| `output_dir` | `str` | Directory to write the extracted PDB file |

**Returns:** `str` — path to the extracted PDB file (`output_dir/LIGAND.pdb`)

**Raises:** `LigandParametrizationError` if the ligand is not found or extraction fails

### Example 20: Extract a Ligand

```python
--8<-- "tests/builder_examples/builder_example_20.py"
```

---

### parametrize_ligand

Parametrize a single ligand using the AMBER/GAFF workflow (antechamber → parmchk2 → tleap).

```python
parametrize_ligand(
    ligand_pdb: str,
    ligand_name: str,
    output_dir: str,
    charge: int = 0,
    charge_method: str = 'bcc',
    multiplicity: int = 1,
    atom_type: str = 'gaff2',
) -> Dict[str, str]
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `ligand_pdb` | `str` | — | Path to the ligand PDB file |
| `ligand_name` | `str` | — | 3-letter residue name (must match PDB) |
| `output_dir` | `str` | — | Directory for output files |
| `charge` | `int` | `0` | Net charge of the ligand |
| `charge_method` | `str` | `'bcc'` | Charge method (see `CHARGE_METHODS`) |
| `multiplicity` | `int` | `1` | Spin multiplicity |
| `atom_type` | `str` | `'gaff2'` | Atom type set — `'gaff2'` or `'gaff'` (see `ATOM_TYPES`) |

**Returns:** Dictionary with paths to generated files:

| Key | Description |
|-----|-------------|
| `'mol2'` | Typed MOL2 file (GAFF/GAFF2 atom types + charges) |
| `'frcmod'` | Force field modification file (missing parameters) |
| `'lib'` | Residue library file for tleap |
| `'prmtop'` | AMBER topology file |
| `'inpcrd'` | AMBER coordinate file |

**Raises:** `LigandParametrizationError` if any step fails

!!! warning "Critical: Variable Name Must Match Residue Name"
    The tleap variable name **must** match the residue name in the PDB.
    For example, `AAA = loadmol2 AAA.mol2` — using a generic name like
    `mol = loadmol2 AAA.mol2` will cause `saveoff` to store the unit under the
    wrong name and tleap will fail to recognize the residue when loading the full system PDB.

---

### parametrize_all_ligands

Detect and parametrize all ligands in a PDB file in one call.

```python
parametrize_all_ligands(
    pdb_file: str,
    output_dir: str,
    charges: Optional[Dict[str, int]] = None,
    charge_method: str = 'bcc',
    atom_type: str = 'gaff2',
) -> Dict[str, Dict[str, str]]
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pdb_file` | `str` | — | PDB file containing ligands |
| `output_dir` | `str` | — | Base directory (each ligand gets a subdirectory) |
| `charges` | `Dict[str, int]` or `None` | `None` | Map of ligand name → net charge (default: 0 for all) |
| `charge_method` | `str` | `'bcc'` | Charge method for antechamber |
| `atom_type` | `str` | `'gaff2'` | Atom type set — `'gaff2'` or `'gaff'` (see `ATOM_TYPES`) |

**Returns:** `Dict[str, Dict[str, str]]` — maps ligand names to their file paths (same structure as `parametrize_ligand()`)

### Example 22: Parametrize All Ligands

```python
--8<-- "tests/builder_examples/builder_example_22.py"
```

---

### get_ligand_2d_image

Generate a publication-quality 2D molecular structure image using RDKit.

```python
get_ligand_2d_image(
    ligand_pdb_or_mol2: str,
    output_image: str,
    width: int = 400,
    height: int = 300,
    *,
    remove_nonpolar_h: bool = True,
    remove_all_h: bool = False,
    dpi: int = 150,
    bond_line_width: float = 2.5,
    atom_label_font_size: int = 0,
    background_color: tuple = (0.11, 0.11, 0.11, 1.0),
    padding: float = 0.15,
    kekulize: bool = True,
    wedge_bonds: bool = True,
    atom_palette: dict | None = None,
    highlight_atoms: list[int] | None = None,
    highlight_color: tuple = (1.0, 0.8, 0.0, 0.3),
    transparent_background: bool = False,
) -> Optional[str]
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `ligand_pdb_or_mol2` | `str` | — | Path to PDB or MOL2 file |
| `output_image` | `str` | — | Path for the output PNG |
| `width` | `int` | `400` | Image width in pixels (before DPI scaling) |
| `height` | `int` | `300` | Image height in pixels (before DPI scaling) |
| `remove_nonpolar_h` | `bool` | `True` | Remove non-polar hydrogens (C-H) for cleaner figures while keeping polar H (O-H, N-H, etc.) |
| `remove_all_h` | `bool` | `False` | Remove **all** explicit hydrogens. Overrides `remove_nonpolar_h` |
| `dpi` | `int` | `150` | DPI multiplier. Use `300` for print-quality images |
| `bond_line_width` | `float` | `2.5` | Thickness of bond lines |
| `atom_label_font_size` | `int` | `0` | Font size for atom labels (0 = auto) |
| `background_color` | `tuple` | `(0.11, 0.11, 0.11, 1.0)` | RGBA background colour (0–1 range) |
| `padding` | `float` | `0.15` | Fractional padding around the molecule (0–1) |
| `kekulize` | `bool` | `True` | Draw aromatic bonds as alternating single/double |
| `wedge_bonds` | `bool` | `True` | Draw stereo wedge/dash bonds |
| `atom_palette` | `dict` | `None` | Custom `{atomic_number: (r, g, b)}` colour map. Uses built-in dark palette when `None` |
| `highlight_atoms` | `list[int]` | `None` | 0-based atom indices to highlight |
| `highlight_color` | `tuple` | `(1.0, 0.8, 0.0, 0.3)` | RGBA colour for highlights |
| `transparent_background` | `bool` | `False` | Produce a transparent PNG |

**Returns:** Path to the generated PNG image, or `None` if generation failed

**Built-in palettes:**

| Palette | Import | Best for |
|---------|--------|----------|
| `_DEFAULT_DARK_PALETTE` | (used automatically) | Dark backgrounds |
| `LIGHT_PALETTE` | `from gatewizard.tools.ligand_parametrization import LIGHT_PALETTE` | White / light backgrounds |

### get_ligand_2d_image_from_pdb_lines

Generate a 2D image directly from PDB HETATM lines (no file required).
All keyword arguments are forwarded to `get_ligand_2d_image`.

```python
get_ligand_2d_image_from_pdb_lines(
    pdb_lines: List[str],
    output_image: str,
    width: int = 400,
    height: int = 300,
    **kwargs,
) -> Optional[str]
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pdb_lines` | `List[str]` | — | HETATM lines for the ligand |
| `output_image` | `str` | — | Path for the output PNG |
| `width` | `int` | `400` | Image width |
| `height` | `int` | `300` | Image height |
| `**kwargs` | | | All keyword options from `get_ligand_2d_image` |

**Returns:** Path to the generated PNG image, or `None` if failed

### Example 24: Generate 2D Structure Images

```python
--8<-- "tests/builder_examples/builder_example_24.py"
```

---

### build_ligand_param_args

Build `--ligand_param` command-line arguments for packmol-memgen.

```python
build_ligand_param_args(
    ligand_files: Dict[str, Dict[str, str]]
) -> List[str]
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `ligand_files` | `Dict[str, Dict[str, str]]` | Ligand name → file paths (as returned by `parametrize_ligand`) |

**Returns:** List of CLI arguments — each ligand produces `['--ligand_param', 'path.frcmod:path.lib']`

!!! info "One Flag Per Ligand"
    packmol-memgen requires a **separate** `--ligand_param` flag per ligand.
    Combining multiple ligands into a single flag will not work.

### build_tleap_ligand_lines

Build tleap input lines to load GAFF/GAFF2 and ligand parameters. Insert these **before** `loadPDB`.

```python
build_tleap_ligand_lines(
    ligand_files: Dict[str, Dict[str, str]],
    atom_type: str = 'gaff2',
) -> str
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `ligand_files` | `Dict[str, Dict[str, str]]` | Ligand name → file paths |
| `atom_type` | `str` | Atom type set — `'gaff2'` or `'gaff'` (default: `'gaff2'`). Determines the `leaprc` source line |

**Returns:** Multi-line string with tleap commands (`source leaprc.gaff2` or `source leaprc.gaff`, `loadamberparams`, `loadoff`)

### Example 21: Build packmol-memgen and tleap Arguments

```python
--8<-- "tests/builder_examples/builder_example_21.py"
```

---

### Builder Integration: `ligand_params` Config Key

When ligands are parametrized (either through the GUI or programmatically), the Builder uses the `ligand_params` configuration key to integrate them into the build process.

**Config Structure:**

```python
config['ligand_params'] = {
    'AAA': {
        'frcmod': '/path/to/AAA/AAA.frcmod',
        'lib': '/path/to/AAA/AAA.lib',
    },
    'BBB': {
        'frcmod': '/path/to/BBB/BBB.frcmod',
        'lib': '/path/to/BBB/BBB.lib',
    },
}
```

**What happens during build:**

1. **packmol-memgen** receives `--ligand_param frcmod:lib` for each ligand, plus `--gaff2` (or `--gaff` depending on the atom type used)
2. **tleap parametrization** loads `source leaprc.gaff2` (or `source leaprc.gaff`), then `loadamberparams` and `loadoff` for each ligand before `loadPDB`
3. The atom type is automatically extracted from the parametrized ligand results, so matching `leaprc` is always used

### Example 23: Full Membrane System with Ligand Parametrization

```python
--8<-- "tests/builder_examples/builder_example_23.py"
```

### Example 26: Atom Type Selection and Recommended Pairings

```python
--8<-- "tests/builder_examples/builder_example_26.py"
```

### Ligand Parametrization Troubleshooting

**"Antechamber failed":**

- Check that AMBER/AmberTools is installed and in `$PATH`
- Verify ligand PDB has proper HETATM records
- Try a different charge method (e.g., `gas` is faster and may work when `bcc` fails)
- Check `logs/antechamber.log` for detailed error messages

**"No ligands detected":**

- Ligands must use HETATM records, not ATOM
- Residue names must NOT be in `STANDARD_RESIDUES` set
- Water (HOH, WAT), ions (NA, CL), and lipids (POPC, etc.) are excluded by design

**"tleap: unknown residue name":**

- The tleap variable name must match the 3-letter residue name exactly
- Ensure `.lib` file was generated with the correct residue name in `saveoff`
- Check that `source leaprc.gaff2` (or `source leaprc.gaff`) is loaded before `loadamberparams`

**"Non-recommended combination" warning (gaff2/bcc or gaff/abcg2):**

- The AMBER manual recommends `gaff/bcc` and `gaff2/abcg2` pairings
- Using `gaff2/bcc` or `gaff/abcg2` may produce suboptimal parameters
- The GUI shows an amber warning label and a confirmation dialog
- To silence the warning, switch to a recommended pairing
- Check `RECOMMENDED_COMBOS` and `NON_RECOMMENDED_COMBOS` constants for the full list

**"resp/esp: Gaussian not found":**

- The `resp` and `esp` charge methods require Gaussian (external QM software)
- Gaussian is **not** included in AmberTools
- Use `bcc` or `abcg2` instead (both use the built-in `sqm` engine)

**"RDKit not available":**

- Install with `conda install -c conda-forge rdkit`
- RDKit is a required dependency for ligand 2D visualization

---

## See Also

- [Preparation Module](preparation.md) - Protonation state analysis
- [Equilibration Module](equilibration.md) - MD equilibration protocols
- [User Guide - Builder](../user-guide.md#builder) - GUI workflow
- [Builder examples](https://github.com/maurobedoya/gatewizard/tree/main/tests/builder_examples) - Numbered scripts included on this page
