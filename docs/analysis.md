# Analysis Features

Detailed documentation of GateWizard's trajectory and energy analysis capabilities.

## Overview

GateWizard provides two main types of analysis:

1. **Structural Analysis**: For trajectory files (.dcd, .xtc, .trr)
2. **Energetic Analysis**: For NAMD log files

Both types support multiple files, custom time assignment, and flexible unit selection.

## Structural Analysis

### Supported Formats

- **Topology**: PSF, PDB, GRO
- **Trajectories**: DCD, XTC, TRR, NetCDF

### Analysis Types

#### RMSD (Root Mean Square Deviation)

**Purpose**: Measure structural deviation from a reference structure.

**Usage**:
```
1. Select RMSD from Analysis Type
2. Set reference frame (default: first frame)
3. Choose atom selection (e.g., "protein and name CA")
4. Enable/disable alignment
5. Run analysis
```

**Output**:
- Plot: Time vs RMSD
- Data: Time points and RMSD values
- Units: Ångstroms (Å) or nanometers (nm)

**Interpretation**:
- Low RMSD (~1-3 Å): Structure remains close to reference
- High RMSD (>5 Å): Significant structural changes
- Plateau: System has equilibrated
- Increasing trend: System still evolving

#### RMSF (Root Mean Square Fluctuation)

**Purpose**: Identify flexible and rigid regions of the protein.

**Usage**:
```
1. Select RMSF from Analysis Type
2. Choose atom selection (typically C-alpha atoms)
3. Run analysis
```

**Output**:
- Plot: Residue number vs RMSF
- Data: Per-residue fluctuation values
- Units: Ångstroms (Å) or nanometers (nm)

**Interpretation**:
- Low RMSF (<1 Å): Rigid regions (core, secondary structures)
- High RMSF (>2 Å): Flexible regions (loops, termini)
- Peaks indicate highly mobile residues

#### Distance Analysis

**Purpose**: Track distances between specific atoms or groups.

**Usage**:
```
1. Select Distances from Analysis Type
2. Define atom pairs:
   - Selection 1: "resid 10 and name CA"
   - Selection 2: "resid 50 and name CA"
3. Run analysis
```

**Output**:
- Plot: Time vs Distance
- Data: Time points and distance values
- Units: Ångstroms (Å) or nanometers (nm)

**Use Cases**:
- Salt bridge formation
- Domain movements
- Ligand-protein interactions
- Conformational changes

#### Radius of Gyration

**Purpose**: Measure overall compactness of the structure.

**Usage**:
```
1. Select Radius of Gyration
2. Choose atom selection (typically all protein atoms)
3. Run analysis
```

**Output**:
- Plot: Time vs Rg
- Data: Time points and Rg values
- Units: Ångstroms (Å) or nanometers (nm)

**Interpretation**:
- Constant Rg: Stable, compact structure
- Increasing Rg: Unfolding or expansion
- Decreasing Rg: Compaction or folding

### Time Management for Trajectories

#### Without Time Assignment

Uses frame numbers × time step:
- Automatically calculated from trajectory metadata
- Suitable for single, continuous trajectories

#### With Time Assignment

Distributes frames evenly across assigned time:

**Example**:
```
File1.dcd: 1000 frames, assigned 10 ns
  → Frame 0: 0.0 ns
  → Frame 500: 5.0 ns
  → Frame 999: 10.0 ns

File2.dcd: 2000 frames, assigned 50 ns
  → Frame 0: 10.0 ns (continues from File1)
  → Frame 1000: 35.0 ns
  → Frame 1999: 60.0 ns
```

**Benefits**:
- Accurate time representation for concatenated runs
- Proper handling of different sampling rates
- Sequential placement of multiple simulations

### Atom Selections

GateWizard uses MDAnalysis selection syntax:

**Common Selections**:
```
protein                      # All protein atoms
protein and name CA          # C-alpha atoms only
protein and backbone         # Backbone atoms (N, CA, C, O)
resid 1:50                   # Residues 1 to 50
resname ALA                  # All alanine residues
name CA and resid 10:20      # C-alpha of residues 10-20
segid PROA                   # Atoms in segment PROA
```

**Advanced Selections**:
```
protein and not name H*              # Protein without hydrogens
(resid 10:20) or (resid 50:60)      # Multiple regions
around 5 resname LIG                 # Within 5Å of ligand
```

### Units and Conversions

**Distance Units**:
- Ångstroms (Å): Default in structural biology
- Nanometers (nm): Common in MD simulations
- Conversion: 1 nm = 10 Å

**Time Units**:
- Picoseconds (ps): Fine-grained time scale
- Nanoseconds (ns): Common for simulation lengths
- Microseconds (µs): Long simulations
- Conversions: 1 ns = 1000 ps, 1 µs = 1000 ns

## Energetic Analysis

### NAMD Log File Analysis

**Supported Data**:
- **ETITLE**: Column headers from log file
- **ENERGY**: Data lines with timesteps and values

**Available Properties**:
- TEMP: Temperature (K)
- TOTAL: Total energy (kcal/mol)
- KINETIC: Kinetic energy (kcal/mol)
- POTENTIAL: Potential energy (kcal/mol)
- ELECT: Electrostatic energy (kcal/mol)
- VDW: van der Waals energy (kcal/mol)
- PRESSURE: Pressure (bar)
- VOLUME: Volume (Ų)

### Time Distribution for Log Files

#### Without Time Assignment

Uses timestep (TS) column:
```
TS = 0:     0 ps
TS = 1000:  2 ps (if timestep=2fs, output every 1000 steps)
TS = 50000: 100 ps
```

#### With Time Assignment

Distributes data points evenly:

**Example**:
```
equilibration.log: 500 points, assigned 5 ns
  → Point 0: 0.0 ns
  → Point 250: 2.5 ns
  → Point 499: 5.0 ns

production.log: 5000 points, assigned 100 ns
  → Point 0: 5.0 ns (continues from equilibration)
  → Point 2500: 55.0 ns
  → Point 4999: 105.0 ns
```

**Benefits**:
- Correct timeline for multi-stage simulations
- Proper spacing for different output frequencies
- Sequential representation of simulation phases

### Energy Units

**kcal/mol** (default in NAMD):
- Standard in US molecular modeling
- 1 kcal/mol = 4.184 kJ/mol

**kJ/mol** (SI unit):
- Standard in international publications
- Conversion applied automatically when selected

**Use Cases**:
- Plot TOTAL energy to check stability
- Monitor TEMP for thermostat performance
- Check PRESSURE for barostat performance
- Compare KINETIC and POTENTIAL for energy distribution

### Multiple Property Plotting

Plot multiple properties simultaneously:

1. Select multiple columns (Ctrl+Click)
2. Each property shown in different color
3. Legend automatically added
4. Y-axis label adjusted for multiple properties

**Example**:
```
Plot: TOTAL, TEMP, PRESSURE
Result: Three lines showing:
  - Energy evolution (left Y-axis)
  - Temperature stability
  - Pressure fluctuations
```

## Plot Customization

### Available Options

**Title and Labels**:
- Plot Title: Custom description
- X-Axis Label: Time (with units)
- Y-Axis Label: Property (with units)

**Units**:
- X-Axis: ps, ns, µs
- Y-Axis (Distance): Å, nm
- Y-Axis (Energy): kcal/mol, kJ/mol

**Appearance**:
- Background Color: Dark Gray, Light Gray, White, Black
- Line Color: Blue, Red, Green, etc.
- Grid: Automatic with transparency

**Ranges**:
- X-Min, X-Max: Limit time range
- Y-Min, Y-Max: Limit value range
- Auto-scaling: Enabled by default

### Export Options

**Plot Export**:
- PNG: Raster image (good for presentations)
- SVG: Vector image (good for publications)
- High resolution available

**Data Export**:
- **CSV**: Spreadsheet format
  ```
  Time(ns),RMSD(Å)
  0.0,0.0
  0.1,0.5
  0.2,0.8
  ```

- **JSON**: Structured format
  ```json
  {
    "time": [0.0, 0.1, 0.2],
    "rmsd": [0.0, 0.5, 0.8],
    "units": {"time": "ns", "rmsd": "Å"}
  }
  ```

- **NumPy**: Python analysis
  ```python
  import numpy as np
  data = np.load('results.npz')
  time = data['time']
  rmsd = data['rmsd']
  ```

## Advanced Features

### Draggable File List

Reorder files by dragging:
1. Click and hold on file bar handle (☰)
2. Drag to new position
3. Release to reorder
4. Time distribution updates automatically

### Real-Time Progress

During analysis:
- Progress bar shows completion
- Status messages update
- Can cancel long-running analyses

### Batch Processing

Analyze multiple trajectories:
1. Add all files at once
2. Assign times for each
3. Run analysis once
4. Data combined automatically

## Troubleshooting Analysis

### Issue: RMSD values too high

**Possible causes**:
- Misaligned structures
- Wrong reference frame
- Including flexible loops

**Solutions**:
- Enable alignment option
- Try different reference frames
- Use more specific atom selection (e.g., secondary structure only)

### Issue: Choppy or discontinuous plots

**Possible causes**:
- Files not in correct order
- Time assignment incorrect
- Missing frames

**Solutions**:
- Reorder files using drag-and-drop
- Verify time assignments
- Check trajectory files for completeness

### Issue: Energy values seem wrong

**Possible causes**:
- Wrong units selected
- Log files not from same simulation
- Parsing errors

**Solutions**:
- Check unit selection (kcal/mol vs kJ/mol)
- Verify log files are sequential
- Check status messages for warnings

### Issue: Memory errors with large trajectories

**Solutions**:
- Analyze in smaller chunks
- Use more specific atom selections
- Increase system RAM
- Close other applications

## Best Practices

### For Structural Analysis

1. **Always check alignment**: Enable for RMSD calculations
2. **Use appropriate selections**: Backbone for overall structure, C-alpha for simplicity
3. **Verify time assignments**: Especially for concatenated trajectories
4. **Export data regularly**: Don't lose analysis results
5. **Check convergence**: RMSD should plateau

### For Energetic Analysis

1. **Sequential file ordering**: Equilibration → Production phases
2. **Assign times accurately**: Match actual simulation times
3. **Plot multiple properties**: Check system stability
4. **Monitor temperature**: Should be stable (±1-2 K)
5. **Check energy drift**: Total energy should be conserved (NVE) or stable (NVT/NPT)

### For Publication

1. **Use consistent units**: Pick one set and stick to it
2. **Export high-resolution plots**: SVG for vector graphics
3. **Save raw data**: CSV or NumPy for later analysis
4. **Document parameters**: Note selections, reference frames, etc.
5. **Check error bars**: Consider statistical analysis of multiple runs

## Example Workflows

### Workflow 1: Basic RMSD Analysis

```
1. Load topology: system.psf
2. Add trajectory: production.dcd
3. Assign time: 100 ns
4. Select: RMSD analysis
5. Atom selection: protein and name CA
6. Reference frame: 0
7. Enable alignment: Yes
8. Run analysis
9. Export: plot as PNG, data as CSV
```

### Workflow 2: Multi-File Energy Analysis

```
1. Add files:
   - min.log → 0.1 ns
   - equil.log → 1.0 ns
   - prod1.log → 50 ns
   - prod2.log → 50 ns
2. Total timeline: 0-101.1 ns
3. Select properties: TOTAL, TEMP, PRESSURE
4. Units: X=ns, Y=kcal/mol
5. Run analysis
6. Export: data as JSON
```

### Workflow 3: Residue Flexibility

```
1. Load topology and trajectory
2. Select: RMSF analysis
3. Atom selection: protein and name CA
4. Run analysis
5. Identify flexible regions (RMSF > 2 Å)
6. Cross-reference with structure
7. Export: plot and data
```

---

*For more examples and detailed tutorials, see the User Guide.*
