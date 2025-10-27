# Analysis Examples# Analysis Examples# Analysis Examples



This directory contains working examples demonstrating the NAMD analysis API for energy and trajectory analysis.



## Available ExamplesThis directory contains working examples demonstrating the NAMD analysis API for energy and trajectory analysis.This directory contains example scripts demonstrating the usage of GateWizard's analysis utilities for MD trajectory and NAMD log analysis.



All examples are self-contained and demonstrate specific analysis capabilities. Each follows a simple pattern similar to equilibration examples.



### Energy Analysis## Available Examples## 📚 Documentation Examples (Matches Official API Docs)



1. **analysis_example_01.py** - Basic Energy Analysis ⭐ RECOMMENDED

   - Simplest example using `plot_energy()`
   - 4-panel plot with default settings
   - Shows Total Energy, Potential & Kinetic, Temperature, and Pressure
   - **Best starting point** for energy analysis



2. **analysis_example_02.py** - Energy Properties Selection### Energy Analysis1. **analysis_example_01_multi_file_trajectory.py** - Example 1: Multi-File Trajectory Analysis

   - Using `plot_properties()` to select specific components

   - Demonstrates plotting TEMP, TOTAL, and PRESSURE together   - Load multiple trajectory files with time scaling

   - Shows focused analysis of selected properties

1. **analysis_example_01.py** - Basic Energy Analysis ⭐ RECOMMENDED   - RMSD analysis with proper alignment

3. **analysis_example_03.py** - Multi-File Energy Analysis

   - Analyzing multiple log files with custom time scaling   - Simplest example using `plot_energy()`   - RMSF analysis with residue labeling and highlighting

   - Essential for multi-stage equilibration protocols
   - Demonstrates continuous energy plot across files
   - 4-panel plot with default settings

### Trajectory Analysis



4. **analysis_example_04.py** - Basic Trajectory Analysis ⭐ RECOMMENDED   - **Best starting point** for energy analysis   - Complete energy analysis from NAMD log files

   - RMSD and RMSF calculations from trajectory files

   - Fundamental structural analysis   - Statistical summaries and convergence assessment

   - Shows backbone RMSD and per-residue RMSF

   - **Best starting point** for trajectory analysis2. **analysis_example_02.py** - Energy Properties Selection   - Custom plot styling



5. **analysis_example_05.py** - Multi-File Trajectory Analysis   - Using `plot_properties()` to select specific components

   - Analyzing multiple trajectory files with time scaling

   - Continuous RMSD across equilibration stages   - Demonstrates plotting TEMP, TOTAL, and PRESSURE together3. **analysis_example_03_publication_figures.py** - Example 3: Publication-Quality Figures

   - Proper time handling for multi-stage protocols

   - Shows focused analysis of selected properties   - White backgrounds and high resolution (300 DPI)

6. **analysis_example_06.py** - Distance and Radius of Gyration

   - Distance measurements between selections   - Clean formatting suitable for journals

   - Protein compactness analysis (Rg)

   - Advanced structural analysis3. **analysis_example_03.py** - Multi-File Energy Analysis   - Suggested figure captions



### Publication-Ready Figures   - Analyzing multiple log files with custom time scaling



7. **analysis_example_07.py** - Publication-Quality Figures   - Essential for multi-stage equilibration protocols### Running Documentation Examples

   - White backgrounds and high resolution (300 DPI)

   - Creates publication-ready energy, RMSD, and RMSF plots   - Demonstrates continuous energy plot across files

   - Clean styling suitable for journal submission

```bash

### Dark Theme Examples (Full Customization) ⭐ NEW

### Trajectory Analysiscd tests/analysis_examples

8. **analysis_example_08.py** - Dark Theme Energy Analysis

   - Complete dark theme customization for energy plots
   - Multiple color schemes (GitHub dark, blue-black, pure black)
   - Shows **all available customization options** for energy analysis

   - Demonstrates separate plots per property

   - RMSD and RMSF calculations from trajectory filespython analysis_example_01_multi_file_trajectory.py

9. **analysis_example_09.py** - Dark Theme RMSD Analysis

   - Full RMSD customization with dark themes   - Fundamental structural analysis

   - Threshold highlighting and reference lines

   - Convergence line customization   - Shows backbone RMSD and per-residue RMSF# Example 2: Energy analysis

   - Multiple dark color schemes demonstrated

   - Shows **all available RMSD plot options**   - **Best starting point** for trajectory analysispython analysis_example_02_energy_analysis.py



10. **analysis_example_10.py** - Dark Theme RMSF Analysis

    - Complete RMSF customization with dark backgrounds

    - Flexibility highlighting for high RMSF regions5. **analysis_example_05.py** - Multi-File Trajectory Analysis# Example 3: Publication-quality figures

    - Different x-axis representations (residue number, type+number, atom index)

    - Vertical lines for marking specific residues   - Analyzing multiple trajectory files with time scalingpython analysis_example_03_publication_figures.py

    - Shows **all available RMSF plot options**

   - Continuous RMSD across equilibration stages```

## Quick Start

   - Proper time handling for multi-stage protocols

### Energy Analysis (2 lines of code)

---

```python

from gatewizard.utils.namd_analysis import EnergyAnalyzer6. **analysis_example_06.py** - Distance and Radius of Gyration



analyzer = EnergyAnalyzer("equilibration_folder/step1_equilibration.log")   - Distance measurements between selections## 🎯 How to Input Custom Time in the GUI

analyzer.plot_energy(save="energy.png")

```   - Protein compactness analysis (Rg)



### Trajectory Analysis (3 lines of code)   - Advanced structural analysis### For Trajectory Files (.dcd, .xtc, .trr)



```python

from gatewizard.utils.namd_analysis import TrajectoryAnalyzer

### Publication-Ready Figures**Structural Analysis Tab:**

analyzer = TrajectoryAnalyzer("equilibration_folder/system.pdb", "equilibration_folder/trajectory.dcd")

analyzer.plot_rmsd(selection="protein and backbone", save="rmsd.png")1. Click "Add Files" to load trajectory files

```

7. **analysis_example_07.py** - Publication-Quality Figures2. Each file shows a time input box: `[___] ns`

## Requirements

   - White backgrounds and high resolution (300 DPI)3. Enter the duration in nanoseconds (e.g., `0.05` for 50 ps)

### Essential

- matplotlib and numpy (for plotting)   - Creates publication-ready energy, RMSD, and RMSF plots4. Click "Auto-Detect Time" to extract from DCD headers automatically

- NAMD log files (.log) for energy analysis

- Topology (PDB, PSF, PRMTOP) and trajectory files (DCD, XTC, TRR) for trajectory analysis   - Clean styling suitable for journal submission5. Files process sequentially (file1: 0-50ps, file2: 50-100ps, etc.)



### For Trajectory Analysis

- MDAnalysis: `conda install -c conda-forge mdanalysis`

## Quick Start**Example:**

## Running Examples

```

All examples use relative paths and expect to be run from this directory:

### Energy Analysis (2 lines of code):::  step1_equilibration.dcd    [0.05] ns  ← Enter 0.05 here

```bash

cd tests/analysis_examples:::  step2_equilibration.dcd    [0.05] ns  ← Enter 0.05 here

python analysis_example_01.py  # Basic energy analysis

python analysis_example_04.py  # Basic trajectory analysis```python:::  step3_equilibration.dcd    [0.05] ns  ← Enter 0.05 here

python analysis_example_08.py  # Dark theme energy (full options)

python analysis_example_09.py  # Dark theme RMSD (full options)from gatewizard.utils.namd_analysis import EnergyAnalyzer```

python analysis_example_10.py  # Dark theme RMSF (full options)

```Result: Continuous time axis from 0 to 150 ps



## Test Data Locationanalyzer = EnergyAnalyzer("step1_equilibration.log")



All examples now read data from the `equilibration_folder/` subdirectory:analyzer.plot_energy(save="energy.png")### For Energy Log Files (.log)

- `equilibration_folder/system.pdb` - Topology file

- `equilibration_folder/step*_equilibration.dcd` - Trajectory files```

- `equilibration_folder/step*_equilibration.log` - NAMD log files

**Energetic Analysis Tab:**

This ensures consistent data access and makes examples self-contained.

### Trajectory Analysis (3 lines of code)1. Click "Add Files" to load NAMD log files

## Features Demonstrated

2. Each file shows a time input box: `[___] ns`

### Energy Analysis

- ✅ Single and multiple log files```python3. Enter the duration in nanoseconds

- ✅ Custom time scaling for multi-file analysis

- ✅ 4-panel summary plotsfrom gatewizard.utils.namd_analysis import TrajectoryAnalyzer4. Click "Auto-Detect Time" to extract from log timestamps

- ✅ Selective property plotting

- ✅ Unit conversions (kcal/mol ↔ kJ/mol, ps ↔ ns ↔ µs)5. Plot will show continuous time across all files

- ✅ Statistical summaries

- ✅ Publication-quality styling (white background)analyzer = TrajectoryAnalyzer("system.pdb", "trajectory.dcd")

- ✅ Dark theme customization (GitHub dark, blue-black, pure black)

- ✅ Separate plots per propertyanalyzer.plot_rmsd(selection="protein and backbone", save="rmsd.png")**Example:**

- ✅ **All plot customization options** (Example 08)

``````

### Trajectory Analysis

- ✅ RMSD calculations with alignment:::  step1_equilibration.log    [0.05] ns  ← 50 ps simulation

- ✅ RMSF per-residue analysis

- ✅ Distance measurements## Requirements:::  step2_equilibration.log    [0.05] ns  ← 50 ps simulation

- ✅ Radius of gyration

- ✅ Multi-file trajectory support:::  step7_production.log       [0.10] ns  ← 100 ps simulation

- ✅ Custom time scaling

- ✅ MDAnalysis selection syntax### Essential```

- ✅ Publication-quality figures (white background)

- ✅ Dark theme figures (multiple color schemes)- matplotlib and numpy (for plotting)Result: Time axis shows 0-50ps (file1), 50-100ps (file2), 100-200ps (file3)

- ✅ Threshold highlighting

- ✅ Reference lines (horizontal and vertical)- NAMD log files (.log) for energy analysis

- ✅ Convergence line customization

- ✅ **All plot customization options** (Examples 09-10)- Topology (PDB, PSF, PRMTOP) and trajectory files (DCD, XTC, TRR) for trajectory analysis**Important Notes:**



## Dark Theme Options Demonstrated- ⏱️ Time input is the **duration** of each file, not cumulative



Examples 08-10 showcase **complete customization** including:### For Trajectory Analysis- 📊 Input unit is always **nanoseconds** (0.05 = 50 ps, 1.0 = 1000 ps)



### Color Schemes- MDAnalysis: `conda install -c conda-forge mdanalysis`- 🔄 Files are processed in order - drag `:::` to reorder

- GitHub Dark (`#0d1117`, `#010409`)

- Blue-Black (`#1a1a2e`, `#16213e`)- 📈 Display units (ps/ns/µs) can be changed in Plot Settings

- Dark Slate (`#2b2d42`, `#1a1b2e`)

- Pure Black (`#000000`)## Running Examples- 🤖 Leave empty or 0 to use default timestep calculation



### Line Colors

- Cyan (`#00d9ff`)

- Magenta (`#ff006e`)All examples use relative paths and expect to be run from this directory:---

- Chartreuse (`#7fff00`)

- Gold (`#ffd700`)

- White (`#ffffff`)

- Custom RGB values```bash## ⭐ Quick Start Examples (Legacy - Still Useful!)



### Customization Optionscd tests/analysis_examples

- Line width and style

- Background and figure colorspython analysis_example_01.py  # Energy analysis### Simple Examples

- Text and grid colors

- Threshold highlightingpython analysis_example_04.py  # Trajectory analysis

- Reference lines (horizontal and vertical)

- Convergence lines```1. **simple_example_01.py** - Energy analysis in 3 lines of code!

- Axis limits

- Figure size and DPI   - Initialize `EnergyAnalyzer` with log file



## See Also## Test Data   - Create energy plot with `plot_energy()`



- [Analysis Module Documentation](../../docs/api/analysis.md) - Complete API reference   - Get statistics with `get_statistics()`

- [Equilibration Examples](../equilibration_examples/) - Setup examples for creating simulation files

- [User Guide](../../docs/user-guide.md) - Complete usage guideThe examples use test data from the `equilibration_folder/` subdirectory:



## Example Pattern- `system.pdb` - Topology file2. **simple_example_02.py** - Trajectory analysis made easy (conceptual)



All examples follow a consistent, simple pattern:- `step*_equilibration.dcd` - Trajectory files3. **simple_example_03.py** - Complete workflow (conceptual)



```python- `step*_equilibration.log` - NAMD log files

from pathlib import Path

from gatewizard.utils.namd_analysis import EnergyAnalyzer  # or TrajectoryAnalyzer### Running Simple Examples



# Get directory## Features Demonstrated

script_dir = Path(__file__).parent

data_dir = script_dir / "equilibration_folder"```bash



# Initialize analyzer### Energy Analysis# Energy analysis (works with demo data)

analyzer = EnergyAnalyzer(data_dir / "file.log")

- ✅ Single and multiple log filespython simple_example_01.py

# Create plot

analyzer.plot_energy(save="output.png")- ✅ Custom time scaling for multi-file analysis```



# Print confirmation- ✅ 4-panel summary plots

print(f"Analysis complete! Plot saved: output.png")

```- ✅ Selective property plotting---



This pattern ensures examples are:- ✅ Unit conversions (kcal/mol ↔ kJ/mol, ps ↔ ns ↔ µs)

- **Simple** - Minimal code, easy to understand

- **Self-contained** - Uses relative paths to equilibration_folder, no configuration needed- ✅ Statistical summaries### For Energy Analysis (EnergyAnalyzer)

- **Reusable** - Easy to adapt for your own data

- **Clear** - One example = one concept- ✅ Publication-quality styling```bash

- **Complete** - Dark theme examples show ALL customization options

# No additional dependencies! Uses built-in Python libraries

### Trajectory Analysis```

- ✅ RMSD calculations with alignment

- ✅ RMSF per-residue analysis### For Trajectory Analysis (TrajectoryAnalyzer)

- ✅ Distance measurements```bash

- ✅ Radius of gyration# Install MDAnalysis

- ✅ Multi-file trajectory supportconda install -c conda-forge mdanalysis

- ✅ Custom time scaling

- ✅ MDAnalysis selection syntax# Or with pip

- ✅ Publication-quality figurespip install MDAnalysis

```

## See Also

### For Plotting

- [Analysis Module Documentation](../../docs/api/analysis.md) - Complete API reference```bash

- [Equilibration Examples](../equilibration_examples/) - Setup examples for creating simulation files# Usually already installed, but if needed:

- [User Guide](../../docs/user-guide.md) - Complete usage guideconda install matplotlib numpy

```

## Example Pattern

## Quick Reference

All examples follow a consistent, simple pattern:

### Energy Analysis API

```python

from pathlib import Path```python

from gatewizard.utils.namd_analysis import EnergyAnalyzer  # or TrajectoryAnalyzerfrom gatewizard.utils.namd_analysis import EnergyAnalyzer



# Get directory# Initialize

script_dir = Path(__file__).parentanalyzer = EnergyAnalyzer("step1_equilibration.log")



# Initialize analyzer# Plot

analyzer = EnergyAnalyzer(script_dir / "file.log")analyzer.plot_energy(save="energy.png")



# Create plot# Statistics

analyzer.plot_energy(save="output.png")stats = analyzer.get_statistics()

print(f"Temperature: {stats['temp']['mean']:.1f} K")

# Print confirmation```

print(f"Analysis complete! Plot saved: output.png")

```### Trajectory Analysis API



This pattern ensures examples are:```python

- **Simple** - Minimal code, easy to understandfrom gatewizard.utils.namd_analysis import TrajectoryAnalyzer

- **Self-contained** - Uses relative paths, no configuration needed

- **Reusable** - Easy to adapt for your own data# Initialize

- **Clear** - One example = one conceptanalyzer = TrajectoryAnalyzer("system.psf", "trajectory.dcd")


# Individual plots
analyzer.plot_rmsd(save="rmsd.png")
analyzer.plot_rmsf(save="rmsf.png")

# Complete analysis in one plot
analyzer.plot_summary(save="full.png")

# Custom distance analysis
analyzer.plot_distances(
    selections={"gate": ("resid 50-70", "resid 150-170")},
    save="distances.png"
)
```

### Progress Monitoring API

```python
from gatewizard.utils.namd_analysis import get_equilibration_progress
from pathlib import Path

progress = get_equilibration_progress(Path("equilibration/namd"))

for stage, data in progress.items():
    print(f"{stage}: {data.status}")
    if data.timing:
        print(f"  Performance: {data.timing.ns_per_day:.4f} ns/day")
```

## Documentation
For complete API documentation, see:
- `docs/api/analysis.md` - Complete API reference with all methods
- Inline documentation in `gatewizard/utils/namd_analysis.py`

## Notes

- Simple examples use GateWizard's high-level wrapper classes
- Original examples demonstrate underlying concepts with full matplotlib code
- All examples include detailed comments and explanations
- Simple API handles edge cases and provides better error messages
- Publication-quality plots (300 DPI) by default

## Tips

1. **Start with simple examples** - They're much easier to understand
2. **Use simple API for production** - Less code = fewer bugs
3. **Reference original examples** - For understanding implementation details
4. **Customize as needed** - All methods accept optional parameters

## Getting Help

If you encounter issues:

1. Check the documentation: `docs/api/analysis.md`
2. Look at simple examples first
3. Verify MDAnalysis is installed for trajectory analysis
4. Check that log files exist and contain ENERGY lines
from gatewizard.tools.equilibration import NAMDEquilibrationManager
from gatewizard.utils.namd_analysis import get_equilibration_progress

# Setup equilibration
eq_manager = NAMDEquilibrationManager(working_dir=Path("equilibration"))

# Monitor progress
progress = get_equilibration_progress(Path("equilibration/namd"))
print(f"Progress: {progress['equilibration_1'].progress_percent:.1f}%")
```
