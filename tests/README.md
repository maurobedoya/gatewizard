# GateWizard Tests

This directory contains automated tests for GateWizard's core functionality.

## Test Files

### Core Functionality Tests

- **`test_builder.py`** - Builder core functionality and documentation examples
  - Core API tests (configuration, lipids, force fields)
  - Force field manager validation
  - All documentation examples (Examples 1-17)
  - Integration tests (require external tools)
  - Auto-discovers examples from `builder_examples/`

- **`test_preparation.py`** - Preparation integration and protein protonation state analysis
  - Chain information extraction
  - pKa value parsing
  - Protein capping functionality

- **`test_equilibration.py`** - NAMD equilibration protocol generation
  - Configuration file generation
  - Restraint file creation
  - Multi-stage equilibration protocols
  - firsttimestep calculation with variable timesteps
  - Run script generation

- **`test_analysis.py`** - NAMD log file analysis
  - Log file parsing
  - Performance metrics extraction
  - Energy/temperature/pressure analysis
  - Log file discovery and pattern matching

### Example Scripts

- **`builder_examples/`** - Contains executable example scripts (01-17)
  - Each example matches the API documentation
  - Automatically tested by `test_builder.py`
  - Can be run individually or as a suite

### Test Data Files

- **`test_multichain.pdb`** - Sample multi-chain protein structure
- **`test_system.pdb`** - Sample system for testing

## Running Tests

### Run All Tests

```bash
# From the main gatewizard directory
python -m pytest tests/

# With verbose output
python -m pytest tests/ -v

# With coverage report
python -m pytest tests/ --cov=gatewizard --cov-report=html
```

### Run Specific Test Files

```bash
# Run builder tests (includes all 17 examples)
python -m pytest tests/test_builder.py -v

# Run only example tests
python -m pytest tests/test_builder.py::TestSystemBuilderExamples -v

# Run specific example (e.g., Example 08)
python -m pytest tests/test_builder.py::TestSystemBuilderExamples::test_individual_examples[08] -v

# Run examples manually (outside pytest)
cd tests && python test_builder.py

# Run preparation tests
python -m pytest tests/test_preparation.py -v

# Run equilibration tests
python -m pytest tests/test_equilibration.py -v

# Run analysis tests
python -m pytest tests/test_analysis.py -v
```

### Run Specific Test Classes or Methods

```bash
# Run a specific test class
python -m pytest tests/test_equilibration.py::TestNAMDEquilibrationManager -v

# Run a specific test method
python -m pytest tests/test_equilibration.py::TestNAMDEquilibrationManager::test_firsttimestep_calculation -v
```

## Test Requirements

All testing dependencies are included in the main environment:

```bash
conda activate gatewizard
# pytest and dependencies are already installed
```

## Test Coverage

The test suite covers:

1. **Preparation Analysis**
   - Multi-chain PDB file handling
   - pKa value extraction
   - Protonation state assignment
   - Protein capping for termini

2. **System Preparation**
   - Force field selection
   - Lipid composition handling
   - System building workflow

3. **Equilibration Protocol**
   - CHARMM-GUI template-based configuration
   - Multi-stage equilibration (6 stages + production)
   - Variable timestep handling
   - Cumulative firsttimestep calculation
   - Restraint gradients across stages
   - GPU/CPU configuration
   - Run script generation with per-stage resources

4. **Trajectory Analysis**
   - NAMD log file parsing
   - Performance metrics extraction
   - Energy time series analysis
   - Temperature and pressure monitoring

## Writing New Tests

When adding new tests:

1. Follow the existing test structure with pytest fixtures
2. Use descriptive test names: `test_<feature>_<specific_behavior>`
3. Include docstrings explaining what is being tested
4. Use tmp_path fixture for file operations
5. Clean up resources after tests
6. Group related tests in classes

Example:
```python
class TestNewFeature:
    """Test the new feature functionality."""
    
    @pytest.fixture
    def sample_data(self, tmp_path):
        """Create sample test data."""
        # Setup code
        return data
    
    def test_basic_functionality(self, sample_data):
        """Test basic feature behavior."""
        # Test code
        assert result == expected
```

## Continuous Integration

Tests can be run automatically via CI/CD pipelines. Ensure all tests pass before merging pull requests.

## Troubleshooting

**Import Errors**: Make sure gatewizard is installed in development mode:
```bash
pip install -e .
```

**File Not Found**: Tests use tmp_path fixtures for temporary files. These are automatically cleaned up.

**Timeout Issues**: Some tests may take longer on slower systems. Adjust pytest timeout settings if needed.

- Sample trajectories (not included, use your own for manual testing)

## Writing New Tests

### Test File Template

```python
"""
Test module for [feature name].
"""
import pytest
from gatewizard.module import function_to_test


def test_basic_functionality():
    """Test basic functionality of [feature]."""
    result = function_to_test(input_data)
    assert result == expected_output


def test_edge_case():
    """Test edge case handling."""
    with pytest.raises(ValueError):
        function_to_test(invalid_input)


@pytest.mark.integration
def test_integration_scenario():
    """Test integration between components."""
    # Setup
    component1 = setup_component1()
    component2 = setup_component2()
    
    # Test
    result = integrate(component1, component2)
    
    # Verify
    assert result.is_valid()
```

### Test Conventions

1. **File naming**: `test_<module>_<feature>.py`
2. **Function naming**: `test_<what_is_being_tested>`
3. **Use descriptive names**: Clear intent from name
4. **One assertion per test**: When possible
5. **Use fixtures**: For common setup
6. **Mark tests**: Use pytest markers for categories

### Pytest Fixtures

```python
import pytest

@pytest.fixture
def sample_pdb():
    """Provide a sample PDB file path."""
    return "tests/test_system.pdb"

@pytest.fixture
def mock_trajectory():
    """Create a mock trajectory for testing."""
    # Create and return mock data
    pass

def test_with_fixture(sample_pdb):
    """Test using fixture."""
    # Use sample_pdb in test
    pass
```

## Continuous Integration

### GitHub Actions (Example)

```yaml
# .github/workflows/tests.yml
name: Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - uses: conda-incubator/setup-miniconda@v2
        with:
          environment-file: environment.yml
          activate-environment: gatewizard
      - name: Run tests
        shell: bash -l {0}
        run: |
          pip install -e .
          pytest tests/
```

## Test Coverage

### Generate Coverage Report

```bash
# HTML report
python -m pytest tests/ --cov=gatewizard --cov-report=html

# View report
firefox htmlcov/index.html  # Or your browser

# Terminal report
python -m pytest tests/ --cov=gatewizard --cov-report=term-missing
```

### Coverage Goals

- **Overall**: >80% coverage
- **Core modules**: >90% coverage
- **GUI modules**: >60% coverage (GUI testing is challenging)

## Debugging Tests

### Run with Debug Output

```bash
# Show print statements
python -m pytest tests/ -s

# Stop on first failure
python -m pytest tests/ -x

# Drop into debugger on failure
python -m pytest tests/ --pdb

# Show local variables on failure
python -m pytest tests/ -l
```

### Run Individual Tests

```bash
# Run just one test for debugging
python -m pytest tests/test_file.py::test_function -v -s
```

## Common Test Issues

### Issue: Display/GUI tests fail in headless environment

**Solution**: Skip GUI tests or set up virtual display
```bash
# Skip GUI tests
pytest tests/ --ignore=tests/test_*_frame.py

# Or use xvfb (Linux)
xvfb-run pytest tests/
```

### Issue: Import errors

**Solution**: Ensure gatewizard is installed in development mode
```bash
pip install -e .
```

### Issue: Test data not found

**Solution**: Run tests from main directory
```bash
cd /path/to/gatewizard
pytest tests/
```

## Test Markers

Use pytest markers to categorize tests:

```python
@pytest.mark.slow
def test_large_trajectory():
    """Test with large trajectory (slow)."""
    pass

@pytest.mark.requires_display
def test_gui_component():
    """Test requiring display (GUI)."""
    pass

# Run only fast tests
pytest tests/ -m "not slow"

# Skip tests requiring display
pytest tests/ -m "not requires_display"
```

## Contributing Tests

When contributing:

1. **Add tests for new features**: Every new feature needs tests
2. **Fix failing tests**: Don't submit with failing tests
3. **Update test data**: Add new test data if needed
4. **Document tests**: Clear docstrings
5. **Check coverage**: Maintain or improve coverage

## Resources

- [Pytest Documentation](https://docs.pytest.org/)
- [Python Testing Best Practices](https://docs.python-guide.org/writing/tests/)
- [MDAnalysis Testing](https://docs.mdanalysis.org/stable/documentation_pages/testing.html)

---

*For questions about testing, contact the developers.*
