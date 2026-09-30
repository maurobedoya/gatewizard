# Analysis Module API

Module for analyzing MD trajectories and NAMD energy data with easy-to-use wrapper classes supporting all GUI options.

## Overview

The analysis module provides three main trajectory/energy classes:

- **`EnergyAnalyzer`** - Parse and plot NAMD energy data with full customization
- **`TrajectoryAnalyzer`** - Calculate and plot RMSD, RMSF, distances, and radius of gyration with complete control
- **`BilayerTrajectoryAnalyzer`** - Calculate area per lipid (**FATSLiM** CLI default; EVAPL experimental; official lipyphilic `AreaPerLipid` via `apl_method='lipyphilic'`) and membrane thickness (lipyphilic)

**Key Features:**

- Simple 2-line API for quick analysis
- Full customization with all GUI options
- Multi-file support with proper time scaling
- Unit conversions (Å/nm, ps/ns/µs, kcal/kJ)
- Publication-quality plots (300 DPI)
- **PlotSpec** — shared JSON plot contract for API matplotlib export and GUI publication PNG

## PlotSpec (energetic plots)

All engines share a **PlotSpec** description (`gatewizard.utils.plot_spec`) consumed by
`matplotlib_renderer.render_energetic()` for publication PNGs:

```python
from gatewizard.utils.plot_spec import build_plot_spec_from_series
from gatewizard.utils.matplotlib_renderer import render_energetic_to_bytes

data = run_energetic_analysis(log_files, properties=["TEMP", "PRESSURE"])
spec = build_plot_spec_from_series(data["series"], layout="grid", cols=2)
png_bytes = render_energetic_to_bytes(data, spec)
```

| Field | Values | Notes |
|-------|--------|-------|
| `layout` | `overlay` / `grid` | Overlay = one axis; grid = one property per panel |
| `sync_x` | bool | Multi-panel shared time window |
| `global` | colors, units, `xlim`, `dpi`, … | Dark-theme defaults match GUI |
| `panels[]` | `key`, `line_color`, `xlim`, `ylim`, … | Per-series overrides |

GUI **Analysis → Pub PNG** calls the same renderer via `POST /analysis-render-plot`.
Amber and GROMACS: `plot_amber_properties()` / `plot_gromacs_properties()`.


```python
from gatewizard.utils.namd_analysis import EnergyAnalyzer, get_equilibration_progress
from gatewizard.utils.trajectory_analysis import TrajectoryAnalyzer
from gatewizard.utils.lipid_bilayer_analysis import (
    BilayerTrajectoryAnalyzer,
    run_bilayer_analysis,
)
```

---

## Using Custom Time Input in the GUI

### For Trajectory Files (.dcd, .xtc, .trr)

When analyzing multiple trajectory files in the GUI **Structural Analysis** tab:

1. **Add Files**: Click "Add Files" or "Add Folder" to load trajectory files
2. **Input Time**: Each file will appear in the list with a time input field showing "Time (ns)"
3. **Enter Duration**: Enter the simulation duration for each file in **nanoseconds**
   - Example: For a 50 ps simulation, enter `0.05` (50 ps = 0.05 ns)
4. **Auto-Detect**: Click "Auto-Detect Time" to automatically extract times from DCD headers
5. **Reorder**: Drag files using the `:::` handle to reorder if needed

**Example File List:**
```
:::  step1_equilibration.dcd    [0.05] ns
:::  step2_equilibration.dcd    [0.05] ns  
:::  step3_equilibration.dcd    [0.05] ns
```
Total time: 0.15 ns (150 ps) continuous

### For Energy Log Files (.log)

When analyzing NAMD log files in the GUI **Energetic Analysis** tab:

1. **Add Files**: Click "Add Files" or "Add Folder" to load NAMD log files
2. **Input Time**: Each file will appear with a time input field
3. **Enter Duration**: Enter the simulation duration for each file in **nanoseconds**
   - Example: For step1_equilibration.log (50 ps), enter `0.05`
4. **Auto-Detect**: Click "Auto-Detect Time" to extract from log file timestamps
5. **Plot**: The time axis will scale correctly across all files

**Example File List:**
```
:::  step1_equilibration.log    [0.05] ns
:::  step2_equilibration.log    [0.05] ns
:::  step7_production.log       [0.10] ns
```

**Important Notes:**

- Time input is the **duration** of each file, not cumulative
- Files are processed in order - time is sequential (file1: 0-50ps, file2: 50-100ps, etc.)
- Leave empty or 0 to use default timestep-based calculation
- Time unit is always **nanoseconds** in the input fields
- You can change display units (ps/ns/µs) in Plot Settings

---

## Class: EnergyAnalyzer

Easy-to-use wrapper for NAMD energy analysis with built-in plotting and customization.

### Constructor

```python
EnergyAnalyzer(
    log_file: Union[Path, str, List[Union[Path, str]]],
    file_times: Optional[Dict[str, float]] = None
)
```

Initialize energy analyzer with NAMD log file(s).

**Parameters:**

`log_file` (Path, str, or List): Path to NAMD log file, or list of paths for multi-file analysis

  - Single file: `"step1_equilibration.log"`
  - Multiple files: `["step1.log", "step2.log", "step3.log"]`

`file_times` (Dict[str, float], optional): Mapping of log filenames to their simulation durations in **nanoseconds**

  - Example: `{"step1.log": 0.05, "step2.log": 0.05, "step3.log": 0.05}`
  - **Important:** Use just the filename (not full path) as the key
  - Time is the **duration** of each file, not cumulative
  - **Note:** In the GUI, you can input times for each file directly in the file list

**Examples:**
```python
# Single file
analyzer = EnergyAnalyzer("step1_equilibration.log")

# Multiple files with custom time
analyzer = EnergyAnalyzer(
    ["step1.log", "step2.log", "step3.log"],
    file_times={
        "step1.log": 0.05,  # 50 ps = 0.05 ns
        "step2.log": 0.05,
        "step3.log": 0.05
    }
)
```

---

### plot_energy

Create a 4-panel energy analysis plot with full customization.

```python
plot_energy(
    energy_units: str = "kcal/mol",
    time_units: str = "ns",
    target_temperature: Optional[float] = None,
    target_pressure: Optional[float] = None,
    bg_color: str = "#2b2b2b",
    fig_bg_color: str = "#212121",
    text_color: str = "Auto",
    show_grid: bool = True,
    title: Optional[str] = None,
    save: Optional[str] = None,
    show: bool = False,
    figsize: tuple = (12, 10),
    dpi: int = 300
)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `energy_units` | str | `"kcal/mol"` | Energy units: `"kcal/mol"` or `"kJ/mol"` |
| `time_units` | str | `"ns"` | Time units: `"ps"`, `"ns"`, or `"µs"` |
| `target_temperature` | float or None | `None` | Target temperature for reference line (K). Defaults to 300.0 K if not provided. |
| `target_pressure` | float or None | `None` | Target pressure for reference line (atm). Defaults to 1.0 atm if not provided. |
| `bg_color` | str | `"#2b2b2b"` | Plot area background color (hex or `"none"`) |
| `fig_bg_color` | str | `"#212121"` | Figure border background color (hex or `"none"`) |
| `text_color` | str | `"Auto"` | Text/axes color (`"Auto"`, color name, or hex) |
| `show_grid` | bool | `True` | Show grid lines on plots |
| `title` | str or None | `None` | Custom main title (auto-generated if None) |
| `save` | str or None | `None` | Filename to save plot |
| `show` | bool | `False` | Display plot interactively |
| `figsize` | tuple | `(12, 10)` | Figure size (width, height) in inches |
| `dpi` | int | `300` | Resolution for saved figure |

**Plots Generated:**

1. **Total Energy** - Convergence over time
2. **Potential & Kinetic** - Energy components
3. **Temperature** - Stability with target line
4. **Pressure** - Fluctuations (or energy components if no pressure)

---

### plot_properties

Create custom plots of selected energy/thermodynamic properties with full control over visualization.

This method provides complete API access to all GUI plotting capabilities, allowing users to:

- Analyze multiple log files with custom time assignments
- Select specific properties to plot
- Combine or separate data from multiple files
- Apply full customization (colors, units, grid, limits, etc.)

```python
plot_properties(
    properties: Union[str, List[str]],
    separate_plots: bool = False,
    energy_units: str = "kcal/mol",
    time_units: str = "ns",
    line_colors: Optional[Union[str, List[str]]] = None,
    bg_color: str = "#2b2b2b",
    fig_bg_color: str = "#212121",
    text_color: str = "Auto",
    grid_color: str = "#444444",
    show_grid: bool = True,
    xlim: Optional[tuple] = None,
    ylim: Optional[tuple] = None,
    title: Optional[str] = None,
    save_prefix: Optional[str] = None,
    show: bool = False,
    figsize: tuple = (10, 6),
    dpi: int = 300
)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `properties` | str or List[str] | - | Property name(s) to plot. See available properties below. |
| `separate_plots` | bool | `False` | If True, create individual plots per file; if False, combine all files on same plot |
| `energy_units` | str | `"kcal/mol"` | Energy units: `"kcal/mol"` or `"kJ/mol"` |
| `time_units` | str | `"ns"` | Time units: `"ps"`, `"ns"`, or `"µs"` |
| `line_colors` | str or List[str] | `None` | Line color(s). Single color or list for multiple files. Auto-assigned if None. |
| `bg_color` | str | `"#2b2b2b"` | Plot area background color (hex, color name, or `"none"`) |
| `fig_bg_color` | str | `"#212121"` | Figure border background color (hex, color name, or `"none"`) |
| `text_color` | str | `"Auto"` | Text/axes color. `"Auto"` auto-detects from bg luminance. |
| `grid_color` | str | `"#444444"` | Grid line color (hex or color name) |
| `show_grid` | bool | `True` | Show grid lines |
| `xlim` | tuple or None | `None` | X-axis limits: `(min, max)` |
| `ylim` | tuple or None | `None` | Y-axis limits: `(min, max)` |
| `title` | str or None | `None` | Custom plot title (auto-generated if None) |
| `save_prefix` | str or None | `None` | Prefix for saved filenames. Files saved as `{prefix}_{property}.png` |
| `show` | bool | `False` | Display plot interactively |
| `figsize` | tuple | `(10, 6)` | Figure size (width, height) in inches |
| `dpi` | int | `300` | Resolution for saved figures |

**Available Properties:**

Use `get_available_properties()` method to see all properties in your log file(s).

**Note:** Property names are **case-insensitive** and support multiple formats. For example, to plot temperature, you can use:

- `"TEMP"`, `"temp"`, `"Temp"` (short names)
- `"Temperature"`, `"TEMPERATURE"`, `"temperature"` (full names)
- Any mixed case variation like `"TeMpErAtUrE"` will also work!

The same applies to all properties - use whatever format is most convenient for you.

| Property | Alternative Names | Description | Units |
|----------|------------------|-------------|-------|
| `"TOTAL"` | `"Total"`, `"total"`, `"Total Energy"`, etc. | Total energy | kcal/mol or kJ/mol |
| `"KINETIC"` | `"Kinetic"`, `"kinetic"`, `"kin"`, `"Kinetic Energy"` | Kinetic energy | kcal/mol or kJ/mol |
| `"POTENTIAL"` | `"Potential"`, `"potential"`, `"pot"`, `"Potential Energy"` | Potential energy | kcal/mol or kJ/mol |
| `"ELECT"` | `"Elect"`, `"elect"`, `"elec"`, `"Electrostatic Energy"` | Electrostatic energy | kcal/mol or kJ/mol |
| `"VDW"` | `"vdw"`, `"Vdw"`, `"Van der Waals Energy"` | van der Waals energy | kcal/mol or kJ/mol |
| `"TEMP"` | `"temp"`, `"Temp"`, `"Temperature"` | Temperature | K |
| `"PRESSURE"` | `"pressure"`, `"Pressure"`, `"press"`, `"Press"` | Pressure | bar |
| `"VOLUME"` | `"volume"`, `"Volume"`, `"vol"`, `"Vol"` | System volume | Å² |
| `"BOND"` | `"bond"`, `"Bond"`, `"Bond Energy"` | Bond energy | kcal/mol or kJ/mol |
| `"ANGLE"` | `"angle"`, `"Angle"`, `"Angle Energy"` | Angle energy | kcal/mol or kJ/mol |
| `"DIHEDRAL"` | `"dihedral"`, `"Dihedral"`, `"Dihedral Energy"` | Dihedral energy | kcal/mol or kJ/mol |
| `"IMPROPER"` | `"improper"`, `"Improper"`, `"Improper Energy"` | Improper energy | kcal/mol or kJ/mol |

### Example 1: Basic Energy Analysis

Simplest energy analysis using `plot_energy()` - 4-panel plot with default settings.

```python
--8<-- "tests/analysis_examples/analysis_example_01.py"
```

**Output figure:**

![Energy Plot](../images/api/energy_analysis_example_01.png)

*Figure: Default energy plot using the specified NAMD .log files.*


---

**Multi-File Color Assignment:**

When plotting multiple files on the same plot, colors are automatically assigned:

- If `line_colors` is a single color: all files use that color
- If `line_colors` is a list: each file gets corresponding color
- If `line_colors=None`: auto-assigned from matplotlib color cycle

**Time Scaling with `file_times`:**

The `file_times` parameter in the constructor allows mapping data points to correct time ranges:
```python
# Example: 3 files with different simulation times
# File 1: 0-50 ns (5000 frames)
# File 2: 50-125 ns (7500 frames)  
# File 3: 125-200 ns (7500 frames)

analyzer = EnergyAnalyzer(
    log_file=["file1.log", "file2.log", "file3.log"],
    file_times=[(0, 50), (50, 125), (125, 200)],
    # or
    #file_times={
    #    "file1.log": 0.05,  # 50 ps = 0.05 ns
    #    "file2.log": 0.05,
    #    "file3.log": 0.05
    #}
)

# Now time axis will correctly show 0-200 ns with proper frame distribution
analyzer.plot_properties(properties="TEMP", time_units="ns")
```

---

### get_available_properties

Get list of all plottable properties available in the loaded log file(s).

```python
get_available_properties() -> List[str]
```

**Returns:** List of property names that can be passed to `plot_properties()`

### Example 2: Get available properties

```python
--8<-- "tests/analysis_examples/analysis_example_02.py"
```

---


### Example 3: Individual energy plots
```python
--8<-- "tests/analysis_examples/analysis_example_03.py"
```

**Example Output:**

![Multiple Energy Plot](../images/api/energy_analysis_example_03.png)

*Figure: Multiple energy properties plot.*

---

### get_statistics

Get statistical summary of all energy components.

```python
get_statistics() -> Dict[str, Dict[str, float]]
```

**Returns:** Dictionary with statistics for each energy term

**Statistics Provided:**

- `mean` - Average value
- `std` - Standard deviation
- `min` - Minimum value
- `max` - Maximum value
- `initial` - First value
- `final` - Last value

**Available Energy Terms:**

- `timestep`, `total`, `kinetic`, `potential`
- `elect`, `vdw`, `boundary`, `misc`
- `temp`, `pressure`, `volume`

### Example 4: Individual energy plots

```python
--8<-- "tests/analysis_examples/analysis_example_04.py"
```

---

## Class: TrajectoryAnalyzer

Easy-to-use wrapper for MD trajectory analysis with built-in plotting and full customization.

### Constructor

```python
TrajectoryAnalyzer(
    topology: Path,
    trajectory: Path | List[Path],
    file_times: Optional[Dict[str, float]] = None
)
```

Initialize trajectory analyzer with topology and trajectory file(s).

**Parameters:**

- `topology` (Path or str): Path to topology file (PSF, PDB, IMPCRD, etc.)
- `trajectory` (Path, str, or List): Path(s) to trajectory file(s) (DCD, XTC, TRR, etc.)
- Single file: `"trajectory.dcd"`
- Multiple files: `["eq1.dcd", "eq2.dcd", "prod.dcd"]`
- `file_times` (Dict[str, float], optional): Mapping of trajectory filenames to their simulation durations in **nanoseconds**

Example: `{"eq1.dcd": 0.05, "eq2.dcd": 0.05, "prod.dcd": 0.05}`

  - **Important:** Use just the filename (not full path) as the key
  - Time is the **duration** of that file, not cumulative
  - **Note:** In the GUI, you can input times for each file directly in the file list

**Examples:**

**Single Trajectory:**
```python
analyzer = TrajectoryAnalyzer("system.pdb", "trajectory.dcd")
```

**Multiple Trajectories with Time Scaling:**
```python
analyzer = TrajectoryAnalyzer(
    "system.pdb",
    ["step1_equilibration.dcd", "step2_equilibration.dcd", "step7_production.dcd"],
    file_times={
        "step1_equilibration.dcd": 0.05,  # 50 ps
        "step2_equilibration.dcd": 0.05,  # 50 ps
        "step7_production.dcd": 0.05      # 50 ps
    }
)
# Time axis will be continuous: [0...0.05] [0.05...0.10] [0.10...0.15] ns
```

---

### calculate_rmsd

Calculate RMSD for selected atoms.

```python
calculate_rmsd(
    selection: str = "protein and backbone",
    reference_frame: int = 0,
    align: bool = True
) -> Dict[str, np.ndarray]
```

**Parameters:**

- `selection` (str): MDTraj selection string
- `reference_frame` (int): Frame to use as reference (0 = first frame)
- `align` (bool): If True, perform structural alignment (rotation + translation) before RMSD.
                 If False, calculate raw coordinate RMSD without alignment.

**Returns:** Dictionary with:

- `'time'`: Time array in **nanoseconds**
- `'rmsd'`: RMSD values in **Angstroms**

**Example:**
```python
data = analyzer.calculate_rmsd("protein and backbone", align=True)
print(f"Time range: {data['time'][0]:.3f} - {data['time'][-1]:.3f} ns")
print(f"RMSD range: {data['rmsd'].min():.2f} - {data['rmsd'].max():.2f} Å")
```

---

### plot_rmsd

Plot RMSD with full customization options including line styling, threshold highlighting, and convergence line control.

```python
plot_rmsd(
    selection: str = "protein and backbone",
    reference_frame: int = 0,
    align: bool = True,
    distance_units: str = "Å",
    time_units: str = "ns",
    line_color: str = "blue",
    line_width: float = 1.2,
    line_style: str = "-",
    bg_color: str = "#2b2b2b",
    fig_bg_color: str = "#212121",
    text_color: str = "Auto",
    show_grid: bool = True,
    xlim: Optional[tuple] = None,
    ylim: Optional[tuple] = None,
    title: Optional[str] = None,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    highlight_threshold: Optional[float] = None,
    highlight_color: str = "orange",
    highlight_alpha: float = 0.2,
    show_convergence: bool = True,
    convergence_color: str = "red",
    convergence_style: str = "--",
    convergence_width: float = 1.5,
    hlines: Optional[List[float]] = None,
    hline_colors: Optional[List[str]] = None,
    hline_styles: Optional[List[str]] = None,
    hline_widths: Optional[List[float]] = None,
    vlines: Optional[List[float]] = None,
    vline_colors: Optional[List[str]] = None,
    vline_styles: Optional[List[str]] = None,
    vline_widths: Optional[List[float]] = None,
    save: Optional[str] = None,
    show: bool = False,
    figsize: tuple = (10, 6),
    dpi: int = 300
)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `selection` | str | `"protein and backbone"` | MDTraj selection string |
| `reference_frame` | int | `0` | Reference frame index |
| `align` | bool | `True` | Perform structural alignment before RMSD |
| `distance_units` | str | `"Å"` | Distance units: `"Å"` or `"nm"` |
| `time_units` | str | `"ns"` | Time units: `"ps"`, `"ns"`, or `"µs"` |
| `line_color` | str | `"blue"` | Plot line color (matplotlib color or hex) |
| `line_width` | float | `1.2` | Line thickness/width |
| `line_style` | str | `"-"` | Line style: `"-"` (solid), `"--"` (dashed), `"-."` (dash-dot), `":"` (dotted) |
| `bg_color` | str | `"#2b2b2b"` | Plot area background (hex or `"none"`) |
| `fig_bg_color` | str | `"#212121"` | Figure border background (hex or `"none"`) |
| `text_color` | str | `"Auto"` | Text/axes color (`"Auto"`, color, or hex) |
| `show_grid` | bool | `True` | Show grid lines |
| `xlim` | tuple or None | `None` | X-axis limits (min, max) |
| `ylim` | tuple or None | `None` | Y-axis limits (min, max) |
| `title` | str or None | `None` | Custom title (auto-generated if None) |
| `xlabel` | str or None | `None` | Custom X-label (auto with units if None) |
| `ylabel` | str or None | `None` | Custom Y-label (auto with units if None) |
| `highlight_threshold` | float or None | `None` | Highlight regions above this RMSD value |
| `highlight_color` | str | `"orange"` | Color for highlight region and line |
| `highlight_alpha` | float | `0.2` | Alpha transparency for highlight fill |
| `show_convergence` | bool | `True` | Show convergence line (mean of last 20% of trajectory) |
| `convergence_color` | str | `"red"` | Color for convergence line |
| `convergence_style` | str | `"--"` | Line style for convergence line |
| `convergence_width` | float | `1.5` | Width of convergence line |
| `hlines` | list or None | `None` | List of Y values for horizontal reference lines |
| `hline_colors` | list or None | `None` | List of colors for horizontal lines |
| `hline_styles` | list or None | `None` | List of line styles for horizontal lines (default: `"--"`) |
| `hline_widths` | list or None | `None` | List of line widths for horizontal lines (default: `1.0`) |
| `vlines` | list or None | `None` | List of X values for vertical reference lines |
| `vline_colors` | list or None | `None` | List of colors for vertical lines |
| `vline_styles` | list or None | `None` | List of line styles for vertical lines (default: `"--"`) |
| `vline_widths` | list or None | `None` | List of line widths for vertical lines (default: `1.0`) |
| `save` | str or None | `None` | Filename to save plot |
| `show` | bool | `False` | Display plot interactively |
| `figsize` | tuple | `(10, 6)` | Figure size (width, height) |
| `dpi` | int | `300` | Resolution for saved figure |

### Example 5: RMSD Plot

```python
--8<-- "tests/analysis_examples/analysis_example_05.py"
```

---

### calculate_rmsf

Calculate RMSF (Root Mean Square Fluctuation) for selected atoms.

```python
calculate_rmsf(
    selection: str = "protein and name CA"
) -> Dict[str, np.ndarray]
```

**Parameters:**

- `selection` (str): MDTraj selection string

**Returns:** Dictionary with:

- `'resids'`: Residue numbers
- `'rmsf'`: RMSF values in **Angstroms**
- `'resnames'`: Residue names (e.g., ALA, GLY, VAL)
- `'atom_indices'`: Atom indices

### Example 6: RMSF Calculation

```python
--8<-- "tests/analysis_examples/analysis_example_06.py"
```

---

### plot_rmsf

Plot RMSF with full customization including line styling, residue labeling, and flexibility highlighting.

```python
plot_rmsf(
    selection: str = "protein and name CA",
    xaxis_type: str = "residue_number",
    show_residue_labels: bool = True,
    residue_name_format: str = "single",
    label_frequency: str = "auto",
    distance_units: str = "Å",
    line_color: str = "blue",
    line_width: float = 1.2,
    line_style: str = "-",
    bg_color: str = "#2b2b2b",
    fig_bg_color: str = "#212121",
    text_color: str = "Auto",
    show_grid: bool = True,
    xlim: Optional[tuple] = None,
    ylim: Optional[tuple] = None,
    title: Optional[str] = None,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    highlight_threshold: Optional[float] = None,
    highlight_color: str = "orange",
    highlight_alpha: float = 0.2,
    hlines: Optional[List[float]] = None,
    hline_colors: Optional[List[str]] = None,
    hline_styles: Optional[List[str]] = None,
    hline_widths: Optional[List[float]] = None,
    vlines: Optional[List[float]] = None,
    vline_colors: Optional[List[str]] = None,
    vline_styles: Optional[List[str]] = None,
    vline_widths: Optional[List[float]] = None,
    save: Optional[str] = None,
    show: bool = False,
    figsize: tuple = (12, 6),
    dpi: int = 300
)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `selection` | str | `"protein and name CA"` | MDTraj selection string |
| `xaxis_type` | str | `"residue_number"` | X-axis type: `"residue_number"`, `"residue_type_number"`, `"atom_index"` |
| `show_residue_labels` | bool | `True` | Show residue labels on X-axis |
| `residue_name_format` | str | `"single"` | Residue format: `"single"` (A, G, V) or `"triple"` (ALA, GLY, VAL) |
| `label_frequency` | str | `"auto"` | Label frequency: `"all"`, `"auto"`, `"every_2"`, `"every_5"`, `"every_10"`, `"every_20"` |
| `distance_units` | str | `"Å"` | Distance units: `"Å"` or `"nm"` |
| `line_color` | str | `"blue"` | Plot line color (matplotlib color or hex) |
| `line_width` | float | `1.2` | Line thickness/width |
| `line_style` | str | `"-"` | Line style: `"-"` (solid), `"--"` (dashed), `"-."` (dash-dot), `":"` (dotted) |
| `highlight_threshold` | float or None | `None` | If set, highlight residues above this RMSF value |
| `highlight_color` | str | `"orange"` | Color for highlight region and threshold line |
| `highlight_alpha` | float | `0.2` | Alpha transparency for highlight fill region |
| `hlines` | list or None | `None` | List of Y values for horizontal reference lines |
| `hline_colors` | list or None | `None` | List of colors for horizontal lines (cycles if fewer than lines) |
| `hline_styles` | list or None | `None` | List of line styles for horizontal lines (default: `"--"`) |
| `hline_widths` | list or None | `None` | List of line widths for horizontal lines (default: 1.0) |
| `vlines` | list or None | `None` | List of X values (residue numbers) for vertical reference lines |
| `vline_colors` | list or None | `None` | List of colors for vertical lines (cycles if fewer than lines) |
| `vline_styles` | list or None | `None` | List of line styles for vertical lines (default: `"--"`) |
| `vline_widths` | list or None | `None` | List of line widths for vertical lines (default: 1.0) |
| *...other params* | | | Same as `plot_rmsd()` for colors, grid, limits, etc. |

**X-Axis Types:**

- `"residue_number"` - Simple numbers: 1, 2, 3, ...
- `"residue_type_number"` - Type + number: ALA1, GLY2, VAL3, ...
- `"atom_index"` - Atom indices: 0, 15, 30, ...

**Label Frequencies:**

- `"all"` - Every residue labeled
- `"auto"` - Smart frequency based on count (<20: all, <50: every_2, <100: every_5, <200: every_10, else: every_20)
- `"every_2"`, `"every_5"`, `"every_10"`, `"every_20"` - Specific intervals

### Example 7: RMSF Plotting

```python
--8<-- "tests/analysis_examples/analysis_example_07.py"
```

**Output figure:**

![RMSF plot](../images/api/trajectory_analysis_example_07.png)

*Figure: RMSF plot for multiple with threshold highlight.*

---

### calculate_distances

Calculate distances between atom selections over time.

```python
calculate_distances(
    selections: Dict[str, tuple]
) -> Dict[str, Dict[str, np.ndarray]]
```

**Parameters:**

- `selections` (dict): Dictionary of `{name: (selection1, selection2)}`

**Returns:** Dictionary with distance data for each named selection:

- `'time'`: Time array in **nanoseconds**
- `'distance'`: Distance array in **Angstroms**

**Example:**
```python
results = analyzer.calculate_distances({
    "gate": ("resid 1-3 and name CA", "resid 18-20 and name CA"),
    "binding_site": ("resid 100-110", "resname LIG")
})

for name, data in results.items():
    print(f"{name}:")
    print(f"  Time: {data['time'][0]:.3f} - {data['time'][-1]:.3f} ns")
    print(f"  Distance: {data['distance'].mean():.2f} ± {data['distance'].std():.2f} Å")
```

---

### plot_distances

Plot distances with full customization including line styling.

```python
plot_distances(
    selections: Dict[str, tuple],
    distance_units: str = "Å",
    time_units: str = "ns",
    line_colors: Optional[List[str]] = None,
    line_width: float = 1.2,
    line_style: str = "-",
    bg_color: str = "#2b2b2b",
    fig_bg_color: str = "#212121",
    text_color: str = "Auto",
    show_grid: bool = True,
    xlim: Optional[tuple] = None,
    ylim: Optional[tuple] = None,
    title: Optional[str] = None,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    show_mean_lines: bool = True,
    hlines: Optional[List[float]] = None,
    hline_colors: Optional[List[str]] = None,
    hline_styles: Optional[List[str]] = None,
    hline_widths: Optional[List[float]] = None,
    vlines: Optional[List[float]] = None,
    vline_colors: Optional[List[str]] = None,
    vline_styles: Optional[List[str]] = None,
    vline_widths: Optional[List[float]] = None,
    save: Optional[str] = None,
    show: bool = False,
    figsize: tuple = (10, 6),
    dpi: int = 300
)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `selections` | dict | - | Dictionary of `{name: (sel1, sel2)}` |
| `distance_units` | str | `"Å"` | Distance units: `"Å"` or `"nm"` |
| `time_units` | str | `"ns"` | Time units: `"ps"`, `"ns"`, or `"µs"` |
| `line_colors` | list or None | `None` | List of colors for each distance pair |
| `line_width` | float | `1.2` | Line thickness/width |
| `line_style` | str | `"-"` | Line style: `"-"` (solid), `"--"` (dashed), `"-."` (dash-dot), `":"` (dotted) |
| `show_mean_lines` | bool | `True` | Show dashed mean distance lines |
| `hlines` | list or None | `None` | List of Y values for horizontal reference lines |
| `hline_colors` | list or None | `None` | List of colors for horizontal lines (cycles if fewer than lines) |
| `hline_styles` | list or None | `None` | List of line styles for horizontal lines (default: `"--"`) |
| `hline_widths` | list or None | `None` | List of line widths for horizontal lines (default: 1.0) |
| `vlines` | list or None | `None` | List of X values (time points) for vertical reference lines |
| `vline_colors` | list or None | `None` | List of colors for vertical lines (cycles if fewer than lines) |
| `vline_styles` | list or None | `None` | List of line styles for vertical lines (default: `"--"`) |
| `vline_widths` | list or None | `None` | List of line widths for vertical lines (default: 1.0) |
| *...other params* | | | Same as `plot_rmsd()` |


### Example 8: Distances

```python
--8<-- "tests/analysis_examples/analysis_example_08.py"
```

**Output figure:**

![Distances plot](../images/api/trajectory_analysis_example_08_distances.png)

*Figure: Distances plot for multiple defined regions on the protein.*


---

### calculate_radius_of_gyration

Calculate radius of gyration over trajectory.

```python
calculate_radius_of_gyration(
    selection: str = "protein"
) -> Dict[str, np.ndarray]
```

**Parameters:**

- `selection` (str): MDTraj selection string

**Returns:** Dictionary with:

- `'time'`: Time array in **nanoseconds**
- `'rg'`: Radius of gyration in **Angstroms**

**Example:**
```python
data = analyzer.calculate_radius_of_gyration("protein")
print(f"Rg range: {data['rg'].min():.2f} - {data['rg'].max():.2f} Å")
print(f"Mean Rg: {data['rg'].mean():.2f} Å")
```

---

### plot_radius_of_gyration

Plot radius of gyration with full customization including line styling and convergence control.

```python
plot_radius_of_gyration(
    selection: str = "protein",
    distance_units: str = "Å",
    time_units: str = "ns",
    line_color: str = "purple",
    line_width: float = 1.2,
    line_style: str = "-",
    bg_color: str = "#2b2b2b",
    fig_bg_color: str = "#212121",
    text_color: str = "Auto",
    show_grid: bool = True,
    xlim: Optional[tuple] = None,
    ylim: Optional[tuple] = None,
    title: Optional[str] = None,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    show_convergence: bool = True,
    convergence_color: str = "red",
    convergence_style: str = "--",
    convergence_width: float = 1.5,
    hlines: Optional[List[float]] = None,
    hline_colors: Optional[List[str]] = None,
    hline_styles: Optional[List[str]] = None,
    hline_widths: Optional[List[float]] = None,
    vlines: Optional[List[float]] = None,
    vline_colors: Optional[List[str]] = None,
    vline_styles: Optional[List[str]] = None,
    vline_widths: Optional[List[float]] = None,
    save: Optional[str] = None,
    show: bool = False,
    figsize: tuple = (10, 6),
    dpi: int = 300
)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `selection` | str | `"protein"` | MDTraj selection string |
| `distance_units` | str | `"Å"` | Distance units: `"Å"` or `"nm"` |
| `time_units` | str | `"ns"` | Time units: `"ps"`, `"ns"`, or `"µs"` |
| `line_color` | str | `"purple"` | Plot line color (matplotlib color or hex) |
| `line_width` | float | `1.2` | Line thickness/width |
| `line_style` | str | `"-"` | Line style: `"-"` (solid), `"--"` (dashed), `"-."` (dash-dot), `":"` (dotted) |
| `show_convergence` | bool | `True` | Show convergence line (mean of last 20% of trajectory) |
| `convergence_color` | str | `"red"` | Color for convergence line |
| `convergence_style` | str | `"--"` | Line style for convergence line |
| `convergence_width` | float | `1.5` | Width of convergence line |
| `hlines` | list or None | `None` | List of Y values for horizontal reference lines |
| `hline_colors` | list or None | `None` | List of colors for horizontal lines (cycles if fewer than lines) |
| `hline_styles` | list or None | `None` | List of line styles for horizontal lines (default: `"--"`) |
| `hline_widths` | list or None | `None` | List of line widths for horizontal lines (default: 1.0) |
| `vlines` | list or None | `None` | List of X values (time points) for vertical reference lines |
| `vline_colors` | list or None | `None` | List of colors for vertical lines (cycles if fewer than lines) |
| `vline_styles` | list or None | `None` | List of line styles for vertical lines (default: `"--"`) |
| `vline_widths` | list or None | `None` | List of line widths for vertical lines (default: 1.0) |
| *...other params* | | | Same as `plot_rmsd()` for colors, grid, limits, etc. |

### Example 9: Radius of gyration

```python
--8<-- "tests/analysis_examples/analysis_example_09.py"
```

---

### plot_summary

Create a 4-panel analysis plot with RMSD, RMSF, Rg, and statistics.

```python
plot_summary(
    save: Optional[str] = None,
    show: bool = False,
    figsize: tuple = (14, 10),
    dpi: int = 300
)
```

**Parameters:**

- `save` (str, optional): Filename to save plot
- `show` (bool): Display plot interactively
- `figsize` (tuple): Figure size
- `dpi` (int): Resolution

**Example:**
```python
analyzer.plot_summary(save="comprehensive_analysis.png")
```

---

## Tips & Best Practices

### 1. Always Provide file_times for Multiple Files

```python
# ✅ GOOD: Proper time scaling
file_times = {
    "step1_equilibration.dcd": 0.05,  # 50 ps
    "step2_equilibration.dcd": 0.05,
    "step3_equilibration.dcd": 0.05
}

# ❌ BAD: Frame-based time (incorrect)
file_times = None  # Will use frame indices
```

### 2. Match Units to Your Analysis

```python
# Short simulations (< 1 ns)
time_units="ps", distance_units="Å"

# Medium simulations (1-100 ns)
time_units="ns", distance_units="Å"

# Long simulations (> 100 ns)
time_units="ns", distance_units="nm"
```

### 3. RMSF Labeling for Different Protein Sizes

```python
# Small protein (< 100 residues)
label_frequency="every_2"

# Medium protein (100-300 residues)
label_frequency="every_5"

# Large protein (> 300 residues)
label_frequency="every_10" or "every_20"
```

### 4. Color Schemes

```python
# Dark theme (presentations)
bg_color="#2b2b2b", text_color="Auto"

# Light theme (publications)
bg_color="#ffffff", text_color="black"

# Transparent (overlays)
bg_color="none", fig_bg_color="none"
```

---

## Function: get_equilibration_progress()

Get progress information from NAMD equilibration logs.

```python
get_equilibration_progress(
    log_directory: Path,
    output_file: Optional[Path] = None
) -> Dict[str, NAMDProgress]
```

**Parameters:**

- `log_directory` (Path): Directory containing NAMD log files
- `output_file` (Path, optional): Save progress summary to JSON

**Returns:** Dictionary of progress information for each stage

**Example:**
```python
from gatewizard.utils.namd_analysis import get_equilibration_progress

progress = get_equilibration_progress("equilibration/")
for stage, info in progress.items():
    print(f"{stage}: {info.progress_percent:.1f}% complete")
```

---

## Line Styling and Customization

All trajectory analysis plotting methods (`plot_rmsd()`, `plot_rmsf()`, `plot_distances()`, `plot_radius_of_gyration()`) support line styling options for full control over plot appearance.

### Common Line Parameters

Available in all plotting methods:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `line_color` | str | varies | Color name (e.g., `"blue"`) or hex code (e.g., `"#1f77b4"`) |
| `line_width` | float | `1.2` | Line thickness (e.g., `0.8`, `1.5`, `2.5`) |
| `line_style` | str | `"-"` | Line style (see options below) |

### Line Style Options

| Style | Value | Visual | Use Case |
|-------|-------|--------|----------|
| Solid | `"-"` | ──────── | Default, standard plots |
| Dashed | `"--"` | ── ── ── | Highlighting, predictions |
| Dash-dot | `"-."` | ── · ── · | Alternative highlighting |
| Dotted | `":"` | · · · · · | Subtle lines, backgrounds |

### RMSD & Radius of Gyration Specific

Additional convergence line control for `plot_rmsd()` and `plot_radius_of_gyration()`:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `show_convergence` | bool | `True` | Toggle convergence line on/off |
| `convergence_color` | str | `"red"` | Color for convergence line |
| `convergence_style` | str | `"--"` | Line style for convergence |
| `convergence_width` | float | `1.5` | Width of convergence line |

### RMSD Specific

Threshold highlighting for `plot_rmsd()`:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `highlight_threshold` | float or None | `None` | Highlight regions above this RMSD value (in Å or nm) |

### RMSF Specific

Flexibility highlighting for `plot_rmsf()`:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `highlight_threshold` | float or None | `None` | Highlight residues above this RMSF value (in Å or nm) |

### Quick Reference Examples

**Thick Dashed Line:**
```python
analyzer.plot_rmsd(line_width=2.5, line_style="--")
```

**Thin Dotted Line:**
```python
analyzer.plot_rmsf(line_width=0.8, line_style=":")
```

**Custom Convergence (Green Dash-dot):**
```python
analyzer.plot_rmsd(
    convergence_color="green",
    convergence_style="-.",
    convergence_width=2.0
)
```

**No Convergence Line:**
```python
analyzer.plot_radius_of_gyration(show_convergence=False)
```

**RMSD with Threshold:**
```python
analyzer.plot_rmsd(highlight_threshold=3.0)  # Highlight > 3.0 Å
```

**RMSF with Flexibility Threshold:**
```python
analyzer.plot_rmsf(highlight_threshold=2.0)  # Highlight > 2.0 Å
```

---

## MDTraj Selection Syntax

### Common Selections

```python
"protein"                          # All protein atoms
"protein and backbone"             # Backbone atoms (N, CA, C, O)
"protein and name CA"              # C-alpha atoms only
"resid 50-150"                     # Residues 50 to 150
"resid 50-70 and name CA"          # C-alpha of residues 50-70
"resname LIG"                      # Ligand with residue name LIG
"protein and not name H*"          # Protein without hydrogens
"within 5 of resname LIG"          # Atoms within 5Å of ligand
```

### Combining Selections

```python
"(resid 1-50) or (resid 100-150)"    # Two regions
"protein and (name CA or name N)"    # Multiple atom names
"resid 25 and name NH1 NH2"          # Multiple atoms in residue
```

---

## Class: BilayerTrajectoryAnalyzer

Lipid bilayer analysis built on [lipyphilic](https://lipyphilic.readthedocs.io/) and MDAnalysis.

### Constructor

```python
BilayerTrajectoryAnalyzer(
    topology: Path,
    trajectory: Path | List[Path],
    file_times: Optional[Dict[str, float]] = None
)
```

Same topology/trajectory interface as `TrajectoryAnalyzer`. Leaflet assignment is performed automatically before each analysis.

### calculate_area_per_lipid

Calculate the area per lipid. Leaflet assignment uses [lipyphilic `AssignLeaflets`](https://lipyphilic.readthedocs.io/).

- **`fatslim` / `auto` (default):** original archived [FATSLiM](https://github.com/FATSLiM/fatslim) CLI via subprocess. Requires a companion conda env — run [`scripts/install_fatslim_env.sh`](https://github.com/maurobedoya/gatewizard/blob/main/scripts/install_fatslim_env.sh) from an API **source checkout** (not from a GUI `/usr/bin` install or a plain `pip install gatewizard`), then set `GATEWIZARD_FATSLIM`. Full GUI-vs-API notes: [analysis.md](../analysis.md#where-the-install-scripts-live-read-this-first). `exclude_sel` → FATSLiM interacting group. GateWizard remains MIT (no FATSLiM import).
- **`evapl`:** GateWizard **EVAPL** via [freud](https://freud.readthedocs.io/) — periodic XY Voronoi, then `exclude_sel` atoms in the leaflet headgroup Z-range shrink lipid cells with successive half-plane clips. **Experimental / not yet validated.** Default `exclude_sel` is `"protein"`.
- **`lipyphilic`:** official [`lipyphilic.AreaPerLipid`](https://lipyphilic.readthedocs.io/). Empty `exclude_sel` works on PyPI 0.12.x. Non-empty `exclude_sel` requires upstream exclude support (PR #164); install `pip install -r requirements-lipyphilic-git.txt` or use `fatslim` / `evapl`.
- **`gridmat`:** GateWizard **GW experimental** GridMAT-style grid (closer to GridMAT-MD.pl, not bit-identical).
- **`gridmat_md`:** original [GridMAT-MD.pl](https://github.com/jalemkul/gridmat-md) via Perl (external).
- **`vtmc`:** GateWizard **GW experimental** Voronoi + Monte Carlo (Mori et al.).

```python
calculate_area_per_lipid(
    lipid_sel: str = "name PO4",
    leaflet_lipid_sel: Optional[str] = None,
    exclude_sel: Optional[str] = "protein",
    exclude_cutoff: float = 30.0,
    exclude_dim: int = 3,
    apl_method: Optional[str] = "fatslim",
    fatslim_nthreads: int = 1,
    fatslim_jobs: int = 1,
    start: Optional[int] = None,
    stop: Optional[int] = None,
    step: Optional[int] = None,
    verbose: bool = False,
) -> Dict[str, np.ndarray]
```

**Parameters:**

| Parameter | Description |
|-----------|-------------|
| `lipid_sel` | Atoms used for Voronoi tessellation / FATSLiM headgroups. MARTINI: `"name GL1 GL2 ROH"`. All-atom: `"name PO4"` or phosphate selections. |
| `leaflet_lipid_sel` | Selection for leaflet assignment. Defaults to `lipid_sel`. |
| `exclude_sel` | Non-lipid occupant atoms (e.g. `"protein"`, a peptide, DNA). Empty/`None` disables exclusion. |
| `exclude_cutoff` | Å cutoff for exclude atoms near the leaflet (`0` = all exclude atoms). Default `30` Å (3.0 nm). Used by lipyphilic AreaPerLipid (when supported), GridMAT, VTMC. **Ignored by FATSLiM and EVAPL.** |
| `exclude_dim` | `3` = 3D distance to leaflet atoms; `1` = \|z − leaflet midplane\|. |
| `apl_method` | `fatslim` / `auto` (default; external FATSLiM CLI), `evapl` (experimental), `lipyphilic` (official `AreaPerLipid`), `gridmat` (GW experimental), `gridmat_md` (external `.pl`), or `vtmc` (GW experimental). |
| `fatslim_nthreads` | FATSLiM `--nthreads` per process (`-1` = all CPUs). Default `1`. Prefer 1–4 when `fatslim_jobs` > 1. |
| `fatslim_jobs` | Parallel frame chunks (`--begin-frame` / `--end-frame`, half-open). Default `1`. Prefer several jobs with small nthreads over one job with all CPUs. |
| `start`, `stop`, `step` | Trajectory frame range. |

**Returns:** `time` (ns), `areas` (n_lipids × n_frames, Å²), `mean_area_per_lipid`, `mean_upper_leaflet`, `mean_lower_leaflet`, `resids`, `resnames`.

See [Example 14: Area per Lipid](#example-14-area-per-lipid).

### calculate_membrane_thickness

Calculate bilayer thickness from interleaflet headgroup distances ([lipyphilic `MembThickness`](https://lipyphilic.readthedocs.io/en/latest/reference/analysis/memb_thickness.html)).

```python
calculate_membrane_thickness(
    lipid_sel: str = "name PO4",
    leaflet_lipid_sel: Optional[str] = None,
    leaflet_filter_sel: Optional[str] = None,
    n_bins: int = 1,
    interpolate: bool = False,
    start: Optional[int] = None,
    stop: Optional[int] = None,
    step: Optional[int] = None,
    verbose: bool = False,
) -> Dict[str, np.ndarray]
```

**Parameters:**

| Parameter | Description |
|-----------|-------------|
| `lipid_sel` | Headgroup atoms for the thickness calculation. |
| `leaflet_lipid_sel` | Selection for leaflet assignment. Defaults to `lipid_sel`. |
| `leaflet_filter_sel` | Optional filter passed to `AssignLeaflets.filter_leaflets()` to exclude species (e.g. cholesterol). |
| `n_bins` | Grid resolution for intrinsic surface construction (`1` = global mean height). |
| `interpolate` | Fill missing grid values when `n_bins` is large (slower). |

**Returns:** `time` (ns), `thickness` (Å).

See [Example 15: Membrane Thickness](#example-15-membrane-thickness).

### plot_area_per_lipid / `plot_membrane_thickness()`

Generate time-series plots with the same styling options as `TrajectoryAnalyzer` plots. See Examples 14 and 15.

### run_bilayer_analysis

JSON-serializable API for GUI and downstream consumers.

```python
run_bilayer_analysis(
    topology_file: Union[str, Path],
    trajectory_files: List[Union[str, Path]],
    analysis_type: str,
    lipid_sel: str = "name PO4",
    leaflet_lipid_sel: Optional[str] = None,
    leaflet_filter_sel: Optional[str] = None,
    n_bins: int = 1,
    interpolate: bool = False,
    file_times: Optional[Dict[str, float]] = None,
    start: Optional[int] = None,
    stop: Optional[int] = None,
    step: Optional[int] = None,
    verbose: bool = False,
) -> Dict[str, Any]
```

**Supported `analysis_type` values:** `area_per_lipid`, `membrane_thickness` (aliases: `apl`, `thickness`, `memb_thickness`).

**Area per lipid return shape:**

```python
{
    "analysis_type": "area_per_lipid",
    "x": [...],           # time (ns)
    "y": [...],           # mean area per lipid (Å²)
    "mean_upper_leaflet": [...],
    "mean_lower_leaflet": [...],
    "lipid_resids": [...],
    "lipid_resnames": [...],
    "per_lipid_areas": [[...], ...],
    "stats": {"mean", "std", "min", "max"},
}
```

**Membrane thickness return shape:**

```python
{
    "analysis_type": "membrane_thickness",
    "x": [...],
    "y": [...],           # thickness (Å)
    "n_bins": 1,
    "stats": {"mean", "std", "min", "max"},
}
```

---

## Examples

Complete working examples are available in [Analysis examples](https://github.com/maurobedoya/gatewizard/tree/main/tests/analysis_examples). Each example demonstrates specific analysis capabilities and can be run directly.

### Example 10: Dark Theme RMSD Analysis - Full Customization


Complete RMSD customization with dark themes showing **all available plot options**.

```python
--8<-- "tests/analysis_examples/analysis_example_10.py"
```

**Output figures** Three dark-themed RMSD plots:

![Dark RMSD 1](../images/api/dark_theme_rmsd_example_10.png)

*Figure: Dark theme RMSD plot*

![Dark RMSD 2](../images/api/dark_theme_rmsd_with_lines_example_10.png)

*Figure: Dark theme RMSD with lines*

![Dark RMSD 3](../images/api/dark_theme_rmsd_minimal_example_10.png)

*Figure: Dark theme RMSD minimal*

---

## Class: OpenMMLogAnalyzer

Parse and analyze **OpenMM StateDataReporter** log files. The interface mirrors
[EnergyAnalyzer](#class-energyanalyzer) so the two can be used interchangeably
in analysis workflows.

```python
from gatewizard.utils.openmm_analysis import OpenMMLogAnalyzer
```

**Constructor**

```python
OpenMMLogAnalyzer(
    log_file: Union[Path, str, List[Union[Path, str]]],
    file_times: Optional[Dict[str, float]] = None,
)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `log_file` | `str \| Path \| list` | Path to a single OpenMM log file, or a list of paths for multi-stage analysis (times are concatenated). |
| `file_times` | `dict` | Optional `{filename: duration_ns}` mapping that overrides the time axis from the log's `"Time (ps)"` column. Useful when logs lack the time column or have restarted counters. |

**Key differences from EnergyAnalyzer:**

| Feature | `EnergyAnalyzer` (NAMD) | `OpenMMLogAnalyzer` (OpenMM) |
|---------|------------------------|------------------------------|
| Native energy units | kcal/mol | **kJ/mol** |
| Color kwarg in `plot_properties` | `line_colors` | **`colors`** |
| Log format | space-delimited `.log` | tab-delimited (StateDataReporter) |

**Available property keys** (depend on what StateDataReporter reports):

| Key | Description |
|-----|-------------|
| `potential` | Potential energy (kJ/mol) |
| `kinetic` | Kinetic energy (kJ/mol) |
| `total` | Total energy (kJ/mol) |
| `temp` | Temperature (K) |
| `volume` | Box volume (nm³) |
| `density` | Density (g/mL) |

---

### get_statistics

```python
analyzer.get_statistics() -> Dict[str, Dict[str, float]]
```

Returns mean / std / min / max / initial / final for every numeric column
that contains data.

**Returns** `{key: {"mean", "std", "min", "max", "initial", "final"}}`

---

### plot_energy

```python
analyzer.plot_energy(
    energy_units: str = "kJ/mol",
    time_units: str = "ns",
    bg_color: str = "#2b2b2b",
    fig_bg_color: str = "#212121",
    text_color: str = "Auto",
    show_grid: bool = True,
    title: Optional[str] = None,
    target_temperature: Optional[float] = None,
    target_pressure: Optional[float] = None,
    save: Optional[str] = None,
    show: bool = False,
    figsize: tuple = (12, 10),
    dpi: int = 300,
)
```

2×2 energy summary plot — total energy, potential + kinetic, temperature,
volume / density.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `energy_units` | `'kJ/mol'` | `'kJ/mol'` or `'kcal/mol'` |
| `time_units` | `'ns'` | `'ps'`, `'ns'`, or `'µs'` |
| `target_temperature` | `None` | Reference temperature (K) drawn as a horizontal dashed line |
| `save` | `None` | Filename for saved figure |
| `show` | `False` | Display interactively |

---

### plot_properties

```python
analyzer.plot_properties(
    properties: Optional[List[str]] = None,
    energy_units: str = "kJ/mol",
    time_units: str = "ns",
    bg_color: str = "#2b2b2b",
    fig_bg_color: str = "#212121",
    text_color: str = "Auto",
    show_grid: bool = True,
    separate_plots: bool = False,
    save: Optional[str] = None,
    show: bool = False,
    figsize: Optional[tuple] = None,
    dpi: int = 300,
    xlim: Optional[tuple] = None,
    ylim: Optional[tuple] = None,
    colors: Optional[List[str]] = None,
)
```

Plot selected properties vs time. If `properties` is not provided, all
non-empty numeric columns are plotted automatically.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `properties` | `None` (all) | List of keys to plot, e.g. `["potential", "temp", "volume"]` |
| `energy_units` | `'kJ/mol'` | `'kJ/mol'` or `'kcal/mol'` |
| `separate_plots` | `False` | Save/show each property as its own figure |
| `colors` | `None` | Line color list (one per property). Unlike `EnergyAnalyzer`, the kwarg is `colors`, not `line_colors`. |
| `save` | `None` | Filename (or prefix when `separate_plots=True`) |
| `show` | `False` | Display interactively |

---

### Example 11: Basic Energy Analysis

```python
--8<-- "tests/analysis_examples/analysis_example_11.py"
```

---

### Example 12: Multi-Stage Analysis

```python
--8<-- "tests/analysis_examples/analysis_example_12.py"
```

---

### Example 13: Custom Property Plots

```python
--8<-- "tests/analysis_examples/analysis_example_13.py"
```

---

### Example 14: Area per Lipid

Calculate and plot the area per lipid using the default **FATSLiM** CLI (`exclude_sel="protein"` → interacting group). Requires companion `fatslim` — see `scripts/install_fatslim_env.sh`. For the experimental freud path use `apl_method="evapl"`.

```python
--8<-- "tests/analysis_examples/analysis_example_14.py"
```

For other force fields, adjust `lipid_sel` to match your headgroup atoms, e.g. `"name GL1 GL2 ROH"` (MARTINI) or `"name PO4"` (all-atom CHARMM).

---

### Example 15: Membrane Thickness

Calculate and plot bilayer thickness as the mean interleaflet headgroup distance ([lipyphilic `MembThickness`](https://lipyphilic.readthedocs.io/en/latest/reference/analysis/memb_thickness.html)) on the same equilibration trajectories.

```python
--8<-- "tests/analysis_examples/analysis_example_15.py"
```

To exclude cholesterol from a mixed bilayer thickness calculation, pass `leaflet_filter_sel="resname DPPC DOPC"` and set `lipid_sel` to the headgroup atoms of the lipids included in the calculation.

---

## Class: GROMACSLogEnergyAnalyzer

```python
from gatewizard.utils.gromacs_analysis import GROMACSLogEnergyAnalyzer

analyzer = GROMACSLogEnergyAnalyzer(log_files, file_times=None, file_strides=None)
```

Parse GROMACS `.log` energy blocks. Same plotting surface as `EnergyAnalyzer`: `get_statistics()`, `plot_energy()`, `plot_properties()`.

## Class: AmberLogEnergyAnalyzer

```python
from gatewizard.utils.amber_analysis import AmberLogEnergyAnalyzer

analyzer = AmberLogEnergyAnalyzer(mdout_files, file_times=None, file_strides=None)
```

Parse Amber `mdout` energy frames. Use `run_amber_energetic_analysis` for the one-shot helper.

### run_energetic_analysis

```python
from gatewizard.utils.namd_analysis import run_energetic_analysis

data = run_energetic_analysis(log_files, properties=["TEMP", "PRESSURE"])
```

One-shot NAMD energetic analysis used by the GUI and scripts. Returns series data for `plot_spec` rendering.

### run_openmm_energetic_analysis

```python
from gatewizard.utils.openmm_analysis import run_openmm_energetic_analysis
```

Same role as `run_energetic_analysis` for OpenMM StateDataReporter logs.

### run_amber_energetic_analysis

```python
from gatewizard.utils.amber_analysis import run_amber_energetic_analysis
```

One-shot Amber mdout energetic analysis.

### parse_namd_log

```python
from gatewizard.utils.namd_analysis import parse_namd_log
```

Low-level NAMD log parser used by `EnergyAnalyzer`.

### prepare_structural_inputs

```python
from gatewizard.utils.trajectory_analysis import prepare_structural_inputs
```

Normalize topology/trajectory lists (split PDB-as-traj, fill missing box dimensions) before `TrajectoryAnalyzer`.

### run_structural_analysis

```python
from gatewizard.utils.trajectory_analysis import run_structural_analysis
```

One-shot RMSD/RMSF/distance/Rg runner used by the Analysis page.

### clear_analysis_cache

```python
analyzer.clear_analysis_cache()
```

Drop cached MDAnalysis universes / prepared frames on a `TrajectoryAnalyzer` after changing stride or inputs.

### run_fatslim_apl

```python
from gatewizard.utils.fatslim_apl import run_fatslim_apl
```

Call the archived FATSLiM CLI (GPLv3 companion env). Discovery: `GATEWIZARD_FATSLIM` → `PATH` → conda env `fatslim-py38`. See [Analysis features](../analysis.md#where-the-install-scripts-live-read-this-first).

### run_gridmat_md_apl

```python
from gatewizard.utils.gridmat_md_apl import run_gridmat_md_apl
```

Run original GridMAT-MD.pl (Perl). Requires `GATEWIZARD_GRIDMAT` or a clone on `PATH`.

## See Also

- [User Guide](../user-guide.md) - Desktop Analysis page
- [Examples](https://github.com/maurobedoya/gatewizard/tree/main/tests/analysis_examples) - Working code examples
- [Troubleshooting](../troubleshooting.md) - Common issues
