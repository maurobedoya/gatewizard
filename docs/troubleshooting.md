# Troubleshooting

Common issues and solutions for GateWizard.

## Installation Issues

### Issue: ImportError with numpy.compat

**Error Message**:
```
ImportError: cannot import name 'OLD_NUMERIC_TYPES' from 'numpy.compat'
```

**Cause**: Version conflict between NumPy and Parmed.

**Solution**:
```bash
conda activate gatewizard
conda install -c conda-forge parmed=4.3.0 --force-reinstall
```

### Issue: pdb4amber command not found

**Error Message**:
```
bash: pdb4amber: command not found
```

**Cause**: AmberTools not installed or environment not activated.

**Solution**:
```bash
# Check if environment is activated
conda activate gatewizard

# Verify AmberTools installation
which pdb4amber

# If not found, reinstall
conda install -c conda-forge ambertools=24
```

### Issue: CustomTkinter not found

**Error Message**:
```
ModuleNotFoundError: No module named 'customtkinter'
```

**Solution**:
```bash
conda activate gatewizard
pip install customtkinter --force-reinstall
```

### Issue: Display/GUI not appearing

**On Linux**:
```bash
# Install tkinter support
sudo apt-get install python3-tk

# Check DISPLAY variable
echo $DISPLAY

# If empty, set it
export DISPLAY=:0
```

**On macOS**:
```bash
# Install python.app for GUI support
conda install -c conda-forge python.app

# Run with pythonw
pythonw -m gatewizard
```

**On Windows WSL**:
```bash
# Install X server on Windows (VcXsrv, X410)
# In WSL, set DISPLAY
export DISPLAY=$(cat /etc/resolv.conf | grep nameserver | awk '{print $2}'):0
```

## Runtime Issues

### Issue: GateWizard crashes on startup

**Check 1**: Verify Python version
```bash
python --version  # Should be 3.8+
```

**Check 2**: Check dependencies
```bash
python -c "import customtkinter"
python -c "import MDAnalysis"
python -c "import parmed"
```

**Check 3**: Run with debug mode
```bash
gatewizard --debug
```

**Solution**: Look at error messages and reinstall problematic packages.

### Issue: GUI appears but is unresponsive

**Possible causes**:
- Large file being processed
- System resources exhausted
- Display server issues

**Solutions**:
1. Wait for operation to complete (check terminal for progress)
2. Close other applications to free memory
3. Restart display server
4. Check system resources: `top` or `htop`

### Issue: Cannot load PDB file

**Error**: "Invalid PDB file" or "Cannot parse file"

**Solutions**:
1. **Verify file format**: Ensure it's a valid PDB file
   ```bash
   head -n 5 your_file.pdb  # Should show ATOM/HETATM lines
   ```

2. **Clean the PDB**: Use Propka with "Clean PDB" option enabled

3. **Check file permissions**:
   ```bash
   ls -l your_file.pdb
   chmod 644 your_file.pdb  # If needed
   ```

4. **Try different source**: If from PDB database, download again

## Analysis Issues

### Issue: RMSD calculation fails

**Error**: "Cannot compute RMSD"

**Possible causes**:
- Topology/trajectory mismatch
- Atom selection returns no atoms
- Missing frames

**Solutions**:
1. **Verify topology matches trajectory**:
   ```bash
   # Check atom counts
   grep "^ATOM" topology.pdb | wc -l
   # Should match trajectory frame size
   ```

2. **Test atom selection**:
   - Try simpler selection: just "protein"
   - Verify selection syntax

3. **Check trajectory integrity**:
   ```python
   import MDAnalysis as mda
   u = mda.Universe("topology.psf", "trajectory.dcd")
   print(f"Frames: {len(u.trajectory)}")
   ```

### Issue: Trajectory loading is very slow

**Causes**:
- Very large trajectory (>10GB)
- Network storage (slow I/O)
- Many atoms selected

**Solutions**:
1. **Use more specific selection**: "protein and name CA" instead of "protein"
2. **Copy to local disk**: Faster than network storage
3. **Split trajectory**: Analyze in smaller chunks
4. **Use frame stride**: Skip frames if possible

### Issue: Energy plot shows weird values

**Symptoms**:
- Extremely large/small values
- Discontinuities
- Unexpected units

**Solutions**:
1. **Check unit selection**:
   - Verify kcal/mol vs kJ/mol
   - Verify ps vs ns vs µs

2. **Verify log files are sequential**:
   - Check file order
   - Verify time assignments

3. **Check for parsing errors**:
   - Look at status messages
   - Verify log files are from NAMD

4. **Plot individual files**:
   - Remove all files
   - Add one at a time
   - Find problematic file

### Issue: Time distribution is wrong

**Symptoms**:
- Files overlap on time axis
- Gaps in timeline
- Wrong total duration

**Solutions**:
1. **Check time assignments**:
   ```
   File1: 10 ns (not 10,000 ps)
   File2: 50 ns
   Result: 0-10 ns, 10-60 ns ✓
   ```

2. **Verify file order**:
   - Use drag-and-drop to reorder
   - Equilibration before production

3. **Check for typos**: 10 vs 100 ns makes big difference

4. **Clear and reassign**: Remove all files, re-add in order

### Issue: Plot is empty or not updating

**Solutions**:
1. **Verify data exists**: Check status messages
2. **Check plot ranges**: Auto-scale or adjust manually
3. **Verify columns selected**: For NAMD analysis
4. **Re-run analysis**: Click "Run Analysis" again

## Memory Issues

### Issue: Out of memory error

**Error**: "MemoryError" or system freezes

**Causes**:
- Trajectory too large for available RAM
- Multiple large analyses in memory
- Memory leak (rare)

**Solutions**:
1. **Reduce atom selection**: Use "name CA" instead of "protein"
2. **Analyze in chunks**: Split trajectory temporally
3. **Close other applications**: Free up RAM
4. **Increase swap space** (Linux):
   ```bash
   sudo fallocate -l 8G /swapfile
   sudo chmod 600 /swapfile
   sudo mkswap /swapfile
   sudo swapon /swapfile
   ```
5. **Use a machine with more RAM**: For very large systems

## Performance Issues

### Issue: Analysis is very slow

**Symptoms**:
- Progress bar barely moving
- Hours for simple calculation
- CPU not fully utilized

**Solutions**:
1. **Use simpler selections**: Fewer atoms = faster
2. **Check system load**: Close background processes
3. **Use SSD storage**: Faster I/O for trajectory reading
4. **Verify not swapping**: Check memory usage with `htop`
5. **Update libraries**: Newer versions may be faster

### Issue: GUI is laggy

**Solutions**:
1. **Reduce plot updates**: Don't update every frame
2. **Close unused tabs**: Free resources
3. **Lower display resolution**: If using remote X
4. **Restart application**: Clear accumulated state

## File/Path Issues

### Issue: Cannot find output directory

**Error**: "Output directory does not exist"

**Solutions**:
1. **Create directory**:
   ```bash
   mkdir -p /path/to/output
   ```
2. **Use absolute path**: /home/user/output not ~/output
3. **Check permissions**: Must have write access
4. **Verify path exists**: No typos in directory name

### Issue: Cannot save results

**Error**: "Permission denied" or "Cannot write file"

**Solutions**:
1. **Check write permissions**:
   ```bash
   ls -ld /path/to/directory
   # Should show 'w' for write permission
   ```
2. **Change permissions**:
   ```bash
   chmod u+w /path/to/directory
   ```
3. **Save to different location**: Your home directory
4. **Run with proper user**: Not as root (usually)

## Propka Issues

### Issue: Propka calculation fails

**Error**: "Propka failed" or "Cannot calculate pKa"

**Solutions**:
1. **Check PDB format**: Must be valid PDB
2. **Clean PDB first**: Remove heteroatoms
3. **Check residue names**: Must be standard amino acids
4. **Update Propka**:
   ```bash
   pip install --upgrade propka
   ```

### Issue: Protonation states seem wrong

**Verify**:
1. **pH value**: Is it set correctly?
2. **pKa values**: Check against literature
3. **Environment effects**: Buried residues have shifted pKa

**Solutions**:
- Review pKa table in results
- Check pKa shifts
- Compare with similar structures
- Manually adjust if needed

## Platform-Specific Issues

### Linux

**Issue**: libGL error
```
libGL error: No matching fbConfigs or visuals found
```
**Solution**:
```bash
sudo apt-get install libgl1-mesa-glx
```

**Issue**: Permission denied for conda
**Solution**:
```bash
# Fix conda ownership
sudo chown -R $USER:$USER ~/miniconda3
```

### macOS

**Issue**: App not trusted
**Solution**:
```bash
# Allow running from unidentified developer
xattr -d com.apple.quarantine /path/to/gatewizard
```

**Issue**: Retina display scaling
**Solution**: Adjust display settings in System Preferences

### Windows (WSL)

**Issue**: X server not connecting
**Solution**:
```bash
# In WSL, set DISPLAY correctly
export DISPLAY=$(cat /etc/resolv.conf | grep nameserver | awk '{print $2}'):0

# In X server (VcXsrv), disable access control
```

**Issue**: File path issues (Windows vs WSL paths)
**Solution**: Use WSL paths (/mnt/c/...) inside WSL

## Getting More Help

If you're still experiencing issues:

1. **Check logs**: Look at `gatewizard_logging.log` in working directory
2. **Run in debug mode**: `gatewizard --debug`
3. **Search issues**: Check GitHub issues (if repository is public)
4. **Provide details**: When asking for help, include:
   - Operating system
   - Python version (`python --version`)
   - GateWizard version (`gatewizard --version`)
   - Error messages (full traceback)
   - Steps to reproduce

5. **Contact developers**:
   - Constanza González: constanza.gonzalez.villagra@gmail.com
   - Mauricio Bedoya: mbedoya@ucm.cl

## FAQ

**Q: How much RAM do I need?**
A: Minimum 4GB, 8GB+ recommended. Large trajectories may need 16GB+.

**Q: Can I run GateWizard remotely?**
A: Yes, with X forwarding: `ssh -X user@server` then run `gatewizard`

**Q: Does GateWizard work with GPU acceleration?**
A: Not currently. Analysis is CPU-based.

**Q: Can I automate GateWizard?**
A: Python API is available for scripting (see API documentation).

**Q: What trajectory formats are supported?**
A: DCD, XTC, TRR, NetCDF, and others via MDAnalysis.

**Q: How do I cite GateWizard?**
A: Citation information to be added.

---

*This troubleshooting guide is continuously updated. Last update: October 2025*
