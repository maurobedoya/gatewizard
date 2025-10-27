# PROPKA Examples

This directory contains all code examples from the PROPKA API documentation.

## Files

- `propka_example_01.py` - PropkaAnalyzer initialization
- `propka_example_02.py` - Run PROPKA analysis
- `propka_example_03.py` - Extract summary from PKA file
- `propka_example_04.py` - Parse summary and analyze ligands
- `propka_example_05.py` - Apply protonation states
- `propka_example_06.py` - Get default protonation states at different pH
- `propka_example_07.py` - Get available states for residue types
- `propka_example_08.py` - Detect disulfide bonds
- `propka_example_09.py` - Apply disulfide bonds
- `propka_example_10.py` - Combined protonation and disulfide workflow
- `propka_example_11.py` - Get residue statistics
- `propka_example_12.py` - Get pH titration curves
- `propka_example_13.py` - Plot pH titration curves (requires matplotlib)
- `propka_example_14.py` - Plot pKa distribution (requires matplotlib and numpy)
- `propka_example_15.py` - Protein capping with ProteinCapper
- `propka_example_16.py` - Protein capping convenience function
- `propka_example_17.py` - Complete workflow with output directory
- `propka_example_18.py` - Workflow with protein capping (includes pdb4amber and HETATM fix)
- `propka_example_19.py` - Generate structures for multiple pH values
- `propka_example_20.py` - Residue mapping after capping (includes pdb4amber and HETATM fix)
- `propka_example_21.py` - Custom protonation states
- `propka_example_22.py` - Filtering and analysis of pKa shifts

**Note:** Examples 18 and 20 demonstrate the complete workflow with protein capping, including:
  - pdb4amber processing for MD simulation preparation
  - ACE/NME HETATM record fix (critical for packmol-memgen compatibility)

## Usage

### Running Individual Examples

```bash
cd tests/propka_examples
python propka_example_01.py
```

### Running All Tests

To run all examples as tests:

```bash
cd tests
pytest test_propka_examples.py -v
pytest test_propka_complex_structures.py -v
```

## Test PDB Files

The examples use different PDB files depending on complexity:

### Simple Examples (protein.pdb)
- Basic PROPKA functionality
- Protonation state assignment
- Disulfide bond detection
- Generated automatically in tests

### Ligand Analysis (6RV3_AB.pdb)
- Multi-chain membrane protein
- Contains ligands (Y01, DMU, PC1, KKQ, K ions)
- Used for ligand atom analysis examples
- Example 4, 14, 22

### Large Structure (8I5B.pdb)
- Human Nav1.7 sodium channel
- Multiple chains (A, B, C)
- Large number of residues
- Used for multi-chain and performance tests

## Import Dependencies

Different examples require different imports:

### Examples 1-14, 17, 19, 21-22
```python
from gatewizard.core.propka import PropkaAnalyzer
```

### Examples 13-14 (plotting)
```python
from gatewizard.core.propka import PropkaAnalyzer
import matplotlib.pyplot as plt
import numpy as np  # Example 14 only
```

### Examples 15-16
```python
from gatewizard.utils.protein_capping import ProteinCapper
# or
from gatewizard.utils.protein_capping import cap_protein
```

### Examples 18, 20
```python
from gatewizard.core.propka import PropkaAnalyzer
from gatewizard.utils.protein_capping import ProteinCapper
```

**Note:** Examples 18 and 20 use the `run_pdb4amber_with_cap_fix()` method which automatically handles both pdb4amber execution and ACE/NME HETATM fixes.

## Special Workflows

### Complete Capped Protein Workflow (Examples 18, 20)

These examples demonstrate the full pipeline for preparing capped proteins:

1. **Add ACE/NME caps** - Protect termini with standard caps
2. **Run PROPKA** - Analyze pKa values on capped structure
3. **Apply protonation** - Set residue names based on pH
4. **Apply disulfide bonds** - Convert CYS to CYX
5. **Run pdb4amber with cap fix** - Add hydrogens and fix ACE/NME HETATM records

**New Unified Method:** `run_pdb4amber_with_cap_fix()`

This method in the `PropkaAnalyzer` class combines pdb4amber execution with optional ACE/NME HETATM fix:

```python
from gatewizard.core.propka import PropkaAnalyzer

analyzer = PropkaAnalyzer()

# Run pdb4amber with automatic cap fix
result = analyzer.run_pdb4amber_with_cap_fix(
    input_pdb="protein_capped_ph7_ss.pdb",
    output_pdb="protein_prepared.pdb",
    fix_caps=True  # Set to False to disable HETATM fix
)

print(f"Success: {result['success']}")
print(f"Fixed {result['hetatm_fixed']} HETATM records")
```

**Why the fix is needed:** pdb4amber converts ACE/NME caps from ATOM to HETATM records, which causes issues with downstream tools like packmol-memgen. The fix converts them back to ATOM records.

**Parameters:**
- `fix_caps=True`: Automatically fix ACE/NME HETATM → ATOM (recommended for capped proteins)
- `fix_caps=False`: Run pdb4amber only, no HETATM fix

**GUI:** The fix is automatically applied when using "Apply States & Run pdb4amber" with protein capping enabled.

## Documentation

These examples correspond to the API documentation in:
`docs/api/propka.md`

Each example demonstrates specific features of the PROPKA module and can be adapted for your own protein analysis workflows.
