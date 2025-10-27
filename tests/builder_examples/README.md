#!/usr/bin/env python3
"""
System Builder Example 06: README

This README provides an overview of the system builder examples.
"""

CONTENT = """# System Builder Examples

This directory contains example scripts demonstrating the usage of the
GateWizard System Builder module for preparing membrane protein systems.

## Overview

The System Builder module provides tools for:
- Building membrane protein systems with lipid bilayers
- Custom lipid composition (symmetric and asymmetric)
- Force field selection and validation
- Salt addition and ionic strength control
- Integration with Propka for protonation states

## Examples

All examples correspond directly to code examples in the API documentation
(`docs/api/system-builder.md`).

### Basic API (01-05)
- **Example 01**: SystemBuilder constructor - Basic initialization
- **Example 02**: set_configuration() - Custom configuration
- **Example 03**: Available water models - ForceFieldManager queries
- **Example 04**: Available protein force fields - ForceFieldManager queries
- **Example 05**: Available lipid force fields - ForceFieldManager queries

### System Preparation (06-11)
- **Example 06**: Simple symmetric membrane (100% POPC)
- **Example 07**: Asymmetric membrane with multiple lipids
- **Example 08**: Complex composition (plasma membrane mimic)
- **Example 09**: Packing only (no parametrization)
- **Example 10**: Custom salt concentration (high salt)
- **Example 11**: No salt (charge neutralization only)

### Validation and Best Practices (12-13)
- **Example 12**: Input validation using validate_system_inputs()
- **Example 13**: Force field validation and recommendations

## Requirements

- packmol-memgen (for system building)
- AmberTools (for parametrization)
- Propka (for pH-dependent protonation)

Examples 01-05, 12-13 can run without external tools (API demonstrations).
Examples 06-11 require packmol-memgen and oriented PDB files.

## Usage

Run any example script directly:

```bash
python system_builder_example_01.py
python system_builder_example_02.py
# ... etc
```

Or use pytest to run all examples as tests:

```bash
pytest ../test_system_builder.py::TestSystemBuilderExamples -v
```

## Example Descriptions

### Example 01: SystemBuilder Constructor
Demonstrates basic SystemBuilder initialization and checking default configuration.
From: Constructor section in documentation.

### Example 02: set_configuration()
Shows how to customize configuration parameters using set_configuration().
From: set_configuration() method documentation.

### Example 03-05: Force Field Queries
Demonstrate querying available water models, protein force fields, and lipid force fields.
From: ForceFieldManager class documentation.

### Example 06: Simple Symmetric Membrane
Basic membrane preparation with 100% POPC (symmetric).
From: prepare_system() Example 1 in documentation.

### Example 07: Asymmetric Membrane
Asymmetric membrane with different lipid compositions per leaflet.
From: prepare_system() Example 2 in documentation.

### Example 08: Plasma Membrane Mimic
Complex multi-lipid system mimicking plasma membrane composition.
From: prepare_system() Example 3 in documentation.

### Example 09: Packing Only
Demonstrates packing without parametrization for visual inspection.
From: prepare_system() Example 4 in documentation.

### Example 10: Custom Salt Concentration
High salt concentration (0.5 M) for ionic strength studies.
From: prepare_system() Example 5 in documentation.

### Example 11: No Salt
Charge neutralization only without extra salt.
From: prepare_system() Example 6 in documentation.

### Example 12: Input Validation
Shows proper validation workflow before system preparation.
From: validate_system_inputs() method documentation.

### Example 13: Force Field Validation
Demonstrates force field combination validation and recommendations.
From: Tips and Best Practices section in documentation.

## Notes

- Examples 06-11 show the API usage but require actual PDB files to execute
- Each example includes commented-out executable code
- All examples match the code shown in `docs/api/system-builder.md`
- Examples are designed for both learning and testing

## See Also

- [System Builder API Documentation](../../docs/api/system-builder.md)
- [User Guide](../../docs/user-guide.md#prepare-tab-system-building)
- [Complete Workflow Example](../../examples/complete_workflow.py)
- [Test Suite](../test_system_builder.py)
"""

def main():
    print(CONTENT)

if __name__ == "__main__":
    main()
