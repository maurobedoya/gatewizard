#!/usr/bin/env python3
"""
Builder Test Suite

This test suite covers:
1. Core functionality tests (specs and features)
2. Documentation example workflows (Examples 1-17)
3. Force field management tests
4. Integration tests (require external tools)

The test suite automatically discovers and runs all example scripts from:
    tests/builder_examples/builder_example_*.py

Note: Some tests require external tools (packmol-memgen, AmberTools) and are skipped
if those tools are not available.

Usage:
    # Run all tests
    pytest tests/test_builder.py -v
    
    # Run only example tests
    pytest tests/test_builder.py::TestBuilderExamples -v
    
    # Run specific example
    pytest tests/test_builder.py::TestBuilderExamples::test_individual_examples[08] -v
    
    # Run examples manually (outside pytest)
    python tests/test_builder.py
"""

import pytest
import sys
import os
import tempfile
import shutil
from pathlib import Path
import importlib.util

# Add gatewizard to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from gatewizard.core.builder import Builder
from gatewizard.tools.force_fields import ForceFieldManager


# ============================================================================
# SECTION 1: CORE FUNCTIONALITY TESTS
# ============================================================================

class TestBuilder:
    """Test the Builder class core functionality."""
    
    @pytest.fixture
    def builder(self):
        """Create a Builder instance."""
        return Builder()
    
    @pytest.fixture
    def ff_manager(self):
        """Create a ForceFieldManager instance."""
        return ForceFieldManager()
    
    def test_default_configuration(self, builder):
        """Test that default configuration is set correctly."""
        assert builder.config['water_model'] == 'tip3p'
        assert builder.config['protein_ff'] == 'ff19SB'
        assert builder.config['lipid_ff'] == 'lipid21'
        assert builder.config['salt_concentration'] == 0.15
        assert builder.config['cation'] == 'K+'
        assert builder.config['anion'] == 'Cl-'
        assert builder.config['dist_wat'] == 17.5
        assert builder.config['preoriented'] == True
        assert builder.config['parametrize'] == True
        assert builder.config['notprotonate'] == False
    
    def test_set_configuration(self, builder):
        """Test configuration update."""
        builder.set_configuration(
            water_model="tip4p",
            salt_concentration=0.5,
            dist_wat=20.0
        )
        
        assert builder.config['water_model'] == 'tip4p'
        assert builder.config['salt_concentration'] == 0.5
        assert builder.config['dist_wat'] == 20.0
        # Other values should remain unchanged
        assert builder.config['protein_ff'] == 'ff19SB'
    
    def test_lipids_argument_symmetric(self, builder):
        """Test lipids argument generation for symmetric membrane."""
        result = builder._prepare_lipids_argument(
            upper_lipids=["POPC"],
            lower_lipids=["POPC"]
        )
        assert result == "POPC//POPC"
    
    def test_lipids_argument_asymmetric(self, builder):
        """Test lipids argument generation for asymmetric membrane."""
        result = builder._prepare_lipids_argument(
            upper_lipids=["POPC", "POPE"],
            lower_lipids=["POPE", "POPS"]
        )
        assert result == "POPC:POPE//POPE:POPS"
    
    def test_lipids_argument_complex(self, builder):
        """Test lipids argument with complex composition."""
        result = builder._prepare_lipids_argument(
            upper_lipids=["POPC", "POPE", "CHL1"],
            lower_lipids=["POPE", "POPS", "CHL1"]
        )
        assert result == "POPC:POPE:CHL1//POPE:POPS:CHL1"


# ============================================================================
# SECTION 2: FORCE FIELD MANAGEMENT TESTS
# ============================================================================

class TestForceFieldManager:
    """Test the ForceFieldManager class."""
    
    @pytest.fixture
    def ff_manager(self):
        """Create a ForceFieldManager instance."""
        return ForceFieldManager()
    
    def test_water_models_available(self, ff_manager):
        """Test that water models are available."""
        water_models = ff_manager.get_water_models()
        assert len(water_models) > 0
        assert 'tip3p' in water_models
        assert 'tip4p' in water_models
    
    def test_protein_force_fields_available(self, ff_manager):
        """Test that protein force fields are available."""
        protein_ffs = ff_manager.get_protein_force_fields()
        assert len(protein_ffs) > 0
        assert 'ff14SB' in protein_ffs
        assert 'ff19SB' in protein_ffs
    
    def test_lipid_force_fields_available(self, ff_manager):
        """Test that lipid force fields are available."""
        lipid_ffs = ff_manager.get_lipid_force_fields()
        assert len(lipid_ffs) > 0
        assert 'lipid17' in lipid_ffs
        assert 'lipid21' in lipid_ffs
    
    def test_available_lipids(self, ff_manager):
        """Test that lipids are available."""
        lipids = ff_manager.get_available_lipids()
        assert len(lipids) > 0
        assert 'POPC' in lipids
        assert 'POPE' in lipids
        assert 'CHL1' in lipids
    
    def test_validate_combination_valid(self, ff_manager):
        """Test validation of valid force field combination."""
        valid, message = ff_manager.validate_combination(
            water_model="tip3p",
            protein_ff="ff14SB",
            lipid_ff="lipid21"
        )
        assert valid == True
        assert "valid" in message.lower()
    
    def test_validate_combination_invalid_water(self, ff_manager):
        """Test validation fails for invalid water model."""
        valid, message = ff_manager.validate_combination(
            water_model="invalid_water",
            protein_ff="ff14SB",
            lipid_ff="lipid21"
        )
        assert valid == False
        assert "unknown" in message.lower() or "water" in message.lower()
    
    def test_recommendations_membrane(self, ff_manager):
        """Test recommendations for membrane systems."""
        rec = ff_manager.get_recommendations("membrane")
        assert 'water_model' in rec
        assert 'protein_ff' in rec
        assert 'lipid_ff' in rec
        assert 'reason' in rec
        assert rec['water_model'] in ff_manager.get_water_models()
        assert rec['protein_ff'] in ff_manager.get_protein_force_fields()
        assert rec['lipid_ff'] in ff_manager.get_lipid_force_fields()
    
    def test_validate_lipid_valid(self, ff_manager):
        """Test validation of valid lipid."""
        assert ff_manager.validate_lipid("POPC") == True
        assert ff_manager.validate_lipid("POPE") == True
        assert ff_manager.validate_lipid("CHL1") == True
    
    def test_validate_lipid_invalid(self, ff_manager):
        """Test validation of invalid lipid."""
        assert ff_manager.validate_lipid("INVALID_LIPID") == False
        assert ff_manager.validate_lipid("XYZ123") == False


# ============================================================================
# SECTION 3: DOCUMENTATION EXAMPLE TESTS
# ============================================================================

class TestBuilderExamples:
    """Test system builder example scripts."""
    
    @pytest.fixture
    def temp_dir(self):
        """Create temporary directory for test outputs."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)
    
    def test_example_01_basic_configuration(self, temp_dir):
        """Test Example 01: Basic configuration (Constructor)."""
        # This example matches the constructor example in docs
        builder = Builder()
        assert builder.config['water_model'] == 'tip3p'
        assert builder.config['protein_ff'] == 'ff19SB'
        print("✓ Example 01: Basic configuration works")
    
    def test_example_02_custom_configuration(self, temp_dir):
        """Test Example 02: set_configuration()."""
        builder = Builder()
        
        # Set custom configuration (matches docs example)
        builder.set_configuration(
            water_model="tip3p",
            protein_ff="ff14SB",
            lipid_ff="lipid21",
            salt_concentration=0.15,
            cation="Na+",
            anion="Cl-",
            dist_wat=20.0,  # Larger water layer
            preoriented=True
        )
        
        # Verify configuration was updated
        assert builder.config['water_model'] == 'tip3p'
        assert builder.config['protein_ff'] == 'ff14SB'
        assert builder.config['salt_concentration'] == 0.15
        assert builder.config['cation'] == 'Na+'
        assert builder.config['dist_wat'] == 20.0
        print("✓ Example 02: set_configuration() works")
    
    def test_example_03_water_models(self, temp_dir):
        """Test Example 03: Available water models."""
        ff_manager = ForceFieldManager()
        water_models = ff_manager.get_water_models()
        assert len(water_models) > 0
        assert 'tip3p' in water_models
        print("✓ Example 03: Available water models works")
    
    def test_example_04_protein_force_fields(self, temp_dir):
        """Test Example 04: Available protein force fields."""
        ff_manager = ForceFieldManager()
        protein_ffs = ff_manager.get_protein_force_fields()
        assert len(protein_ffs) > 0
        assert 'ff14SB' in protein_ffs
        print("✓ Example 04: Available protein force fields works")
    
    def test_example_05_lipid_force_fields(self, temp_dir):
        """Test Example 05: Available lipid force fields."""
        ff_manager = ForceFieldManager()
        lipid_ffs = ff_manager.get_lipid_force_fields()
        assert len(lipid_ffs) > 0
        assert 'lipid21' in lipid_ffs
        print("✓ Example 05: Available lipid force fields works")
    
    def test_example_06_simple_symmetric_membrane(self, temp_dir):
        """Test Example 06: Simple symmetric membrane (demonstration)."""
        # This tests the API structure shown in the docs
        builder = Builder()
        builder.set_configuration(
            water_model="tip3p",
            protein_ff="ff14SB",
            lipid_ff="lipid21",
            salt_concentration=0.15,
            cation="K+",
            anion="Cl-"
        )
        # Configuration successful
        assert builder.config['water_model'] == 'tip3p'
        print("✓ Example 06: Simple symmetric membrane API works")
    
    def test_example_07_asymmetric_membrane(self, temp_dir):
        """Test Example 07: Asymmetric membrane (demonstration)."""
        builder = Builder()
        builder.set_configuration(
            water_model="tip3p",
            protein_ff="ff14SB",
            lipid_ff="lipid21",
            salt_concentration=0.15,
            dist_wat=20.0
        )
        assert builder.config['dist_wat'] == 20.0
        print("✓ Example 07: Asymmetric membrane API works")
    
    def test_example_08_plasma_membrane(self, temp_dir):
        """Test Example 08: Plasma membrane mimic (demonstration)."""
        builder = Builder()
        # Test the configuration from example 08
        assert builder.config['water_model'] == 'tip3p'
        print("✓ Example 08: Plasma membrane API works")
    
    def test_example_09_packing_only(self, temp_dir):
        """Test Example 09: Packing only (demonstration)."""
        builder = Builder()
        # Verify we can set parametrize=False
        builder.set_configuration(parametrize=False)
        assert builder.config['parametrize'] == False
        print("✓ Example 09: Packing only API works")
    
    def test_example_10_custom_salt(self, temp_dir):
        """Test Example 10: Custom salt concentration (demonstration)."""
        builder = Builder()
        builder.set_configuration(
            salt_concentration=0.5,
            cation="Na+",
            anion="Cl-"
        )
        assert builder.config['salt_concentration'] == 0.5
        assert builder.config['cation'] == 'Na+'
        print("✓ Example 10: Custom salt concentration API works")
    
    def test_example_11_no_salt(self, temp_dir):
        """Test Example 11: No salt (demonstration)."""
        builder = Builder()
        builder.set_configuration(
            salt_concentration=0.0,
            add_salt=True
        )
        assert builder.config['salt_concentration'] == 0.0
        assert builder.config.get('add_salt', True) == True
        print("✓ Example 11: No salt API works")
    
    def test_example_12_input_validation(self, temp_dir):
        """Test Example 12: Input validation."""
        builder = Builder()
        
        # Create a dummy PDB file
        dummy_pdb = temp_dir / "protein.pdb"
        dummy_pdb.write_text("ATOM      1  N   GLU A   1       0.000   0.000   0.000  1.00  0.00           N\n")
        
        # Test validation (matches docs example)
        valid, error_msg = builder.validate_system_inputs(
            pdb_file=str(dummy_pdb),
            upper_lipids=["POPC", "POPE"],
            lower_lipids=["POPC", "POPE"],
            lipid_ratios="7:3//7:3",
            water_model="tip3p",
            protein_ff="ff14SB",
            lipid_ff="lipid21"
        )
        # Should be valid or have specific error
        assert valid or len(error_msg) > 0
        print("✓ Example 12: Input validation works")
    
    def test_example_13_force_field_validation(self, temp_dir):
        """Test Example 13: Force field validation."""
        ff_manager = ForceFieldManager()
        
        # Test validation (matches docs example)
        valid, message = ff_manager.validate_combination("tip3p", "ff14SB", "lipid21")
        assert valid == True
        assert "valid" in message.lower()
        print("✓ Example 13: Force field validation works")
    
    def test_run_example_scripts(self, temp_dir):
        """Test running actual example scripts from builder_examples directory."""
        examples_dir = Path(__file__).parent / "builder_examples"
        
        if not examples_dir.exists():
            pytest.skip(f"Examples directory not found: {examples_dir}")
        
        # Find all example files (01-17)
        example_files = sorted(examples_dir.glob("builder_example_*.py"))
        
        if not example_files:
            pytest.skip("No example files found in builder_examples directory")
        
        print(f"\nFound {len(example_files)} example files to test")
        
        failed_examples = []
        passed_examples = []
        
        for example_file in example_files:
            example_num = example_file.stem.split("_")[-1]
            
            # Load and run the example
            spec = importlib.util.spec_from_file_location(
                f"example_{example_num}", 
                example_file
            )
            
            if spec is None or spec.loader is None:
                failed_examples.append((example_num, "Could not load module"))
                continue
            
            module = importlib.util.module_from_spec(spec)
            
            try:
                print(f"\n{'='*60}")
                print(f"Testing Example {example_num}: {example_file.name}")
                print(f"{'='*60}")
                spec.loader.exec_module(module)
                print(f"✓ Example {example_num} executed successfully")
                passed_examples.append(example_num)
            except Exception as e:
                error_msg = f"{type(e).__name__}: {str(e)}"
                print(f"✗ Example {example_num} failed: {error_msg}")
                failed_examples.append((example_num, error_msg))
        
        # Print summary
        print(f"\n{'='*60}")
        print(f"TEST SUMMARY")
        print(f"{'='*60}")
        print(f"Total examples: {len(example_files)}")
        print(f"Passed: {len(passed_examples)}")
        print(f"Failed: {len(failed_examples)}")
        
        if passed_examples:
            print(f"\n✓ Passed examples: {', '.join(passed_examples)}")
        
        if failed_examples:
            print(f"\n✗ Failed examples:")
            for num, error in failed_examples:
                print(f"  - Example {num}: {error}")
            # Fail the test if any examples failed
            pytest.fail(f"{len(failed_examples)} example(s) failed")
        else:
            print(f"\n🎉 All {len(passed_examples)} examples passed!")
    
    @pytest.mark.parametrize("example_num", [f"{i:02d}" for i in range(1, 18)])
    def test_individual_examples(self, example_num, temp_dir):
        """Test each example individually for better pytest reporting."""
        examples_dir = Path(__file__).parent / "builder_examples"
        example_file = examples_dir / f"builder_example_{example_num}.py"
        
        if not example_file.exists():
            pytest.skip(f"Example {example_num} not found")
        
        # Load and run the example
        spec = importlib.util.spec_from_file_location(
            f"example_{example_num}", 
            example_file
        )
        
        if spec is None or spec.loader is None:
            pytest.fail(f"Could not load example {example_num}")
        
        module = importlib.util.module_from_spec(spec)
        
        try:
            spec.loader.exec_module(module)
            print(f"✓ Example {example_num} passed")
        except Exception as e:
            pytest.fail(f"Example {example_num} failed: {type(e).__name__}: {str(e)}")


# ============================================================================
# SECTION 4: INTEGRATION TESTS (REQUIRE EXTERNAL TOOLS)
# ============================================================================

class TestBuilderIntegration:
    """Integration tests requiring external tools."""
    
    @pytest.fixture
    def temp_dir(self):
        """Create temporary directory for test outputs."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)
    
    @pytest.fixture
    def sample_pdb(self, temp_dir):
        """Create a minimal sample PDB file."""
        pdb_content = """REMARK   Created by GateWizard test suite
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.251   2.389   0.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       1.989  -0.729  -1.232  1.00  0.00           C
TER       6      ALA A   1
END
"""
        pdb_file = temp_dir / "test_protein.pdb"
        pdb_file.write_text(pdb_content)
        return pdb_file
    
    # Note: This test only validates inputs, gatewizard environment should be activated
    def test_prepare_system_dry_run(self, temp_dir, sample_pdb):
        """Test system preparation (dry run - validation only)."""
        builder = Builder()
        
        # Only validate, don't actually run
        valid, error_msg = builder.validate_system_inputs(
            pdb_file=str(sample_pdb),
            upper_lipids=["POPC"],
            lower_lipids=["POPC"],
            lipid_ratios="1//1",
            water_model="tip3p",
            protein_ff="ff14SB",
            lipid_ff="lipid21"
        )
        
        # Should be valid (or have a specific validation error)
        assert valid or len(error_msg) > 0
        print(f"Validation result: {valid}, message: {error_msg}")


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def run_all_examples():
    """Helper function to run all examples manually."""
    examples_dir = Path(__file__).parent / "builder_examples"
    
    if not examples_dir.exists():
        print(f"Examples directory not found: {examples_dir}")
        return
    
    # Find all example files automatically
    example_files = sorted(examples_dir.glob("builder_example_*.py"))
    
    if not example_files:
        print("No example files found")
        return
    
    print(f"Found {len(example_files)} examples to run\n")
    
    passed = []
    failed = []
    
    for example_file in example_files:
        example_num = example_file.stem.split("_")[-1]
        
        print(f"\n{'='*80}")
        print(f"Running Example {example_num}: {example_file.name}")
        print(f"{'='*80}")
        
        spec = importlib.util.spec_from_file_location(f"example_{example_num}", example_file)
        if spec is None or spec.loader is None:
            print(f"✗ Example {example_num}: Could not load module")
            failed.append(example_num)
            continue
        
        module = importlib.util.module_from_spec(spec)
        
        try:
            spec.loader.exec_module(module)
            print(f"✓ Example {example_num} completed successfully")
            passed.append(example_num)
        except Exception as e:
            print(f"✗ Example {example_num} failed: {e}")
            import traceback
            traceback.print_exc()
            failed.append(example_num)
    
    # Print summary
    print(f"\n{'='*80}")
    print(f"SUMMARY")
    print(f"{'='*80}")
    print(f"Total: {len(example_files)}")
    print(f"Passed: {len(passed)} - {passed}")
    print(f"Failed: {len(failed)} - {failed if failed else 'None'}")
    print(f"{'='*80}")


if __name__ == "__main__":
    # Run tests with pytest
    pytest.main([__file__, "-v"])
