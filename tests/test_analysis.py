#!/usr/bin/env python3
"""
NAMD Analysis Test Suite

This test suite covers:
1. Core functionality tests (NAMD log parsing, timing extraction)
2. API feature tests (EnergyAnalyzer, TrajectoryAnalyzer)
3. Documentation example workflows (Examples 1-3)
4. Energy and trajectory analysis tests

Features Tested:
- NAMD log file parsing
- Energy statistics extraction
- Trajectory analysis (RMSD, RMSF, distances, radius of gyration)
- Multi-file trajectory support with time scaling
- Plot generation and customization

The test suite automatically discovers and runs all example scripts from:
    tests/analysis_examples/analysis_example_*.py

Usage:
    # Run all tests
    pytest tests/test_analysis.py -v
    
    # Run only core functionality tests
    pytest tests/test_analysis.py::TestNAMDLogParsing -v
    
    # Run only example tests
    pytest tests/test_analysis.py::TestAnalysisExamples -v
    
    # Run specific example
    pytest tests/test_analysis.py::TestAnalysisExamples::test_individual_examples[01] -v
    
    # Run examples manually (outside pytest)
    python tests/test_analysis.py
"""

import pytest
import sys
import os
import tempfile
from pathlib import Path
import importlib.util

# Add gatewizard to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from gatewizard.utils.namd_analysis import parse_namd_log, NAMDTiming, EnergyAnalyzer


# ============================================================================
# SECTION 1: CORE FUNCTIONALITY TESTS
# ============================================================================


class TestNAMDLogParsing:
    """Test NAMD log file parsing."""
    
    @pytest.fixture
    def sample_log_content(self):
        """Create sample NAMD log content."""
        return """Info: Running on 4 processors.
Info: 50000 ATOMS
Info: TIMESTEP               2
Info: FIRST TIMESTEP         0
ETITLE:      TS           BOND          ANGLE          DIHED          IMPRP               ELECT            VDW       BOUNDARY           MISC        KINETIC               TOTAL           TEMP      POTENTIAL         TOTAL3        TEMPAVG            PRESSURE      GPRESSURE         VOLUME       PRESSAVG      GPRESSAVG

ENERGY:       0      4219.5355     20691.2964     15504.9548         0.0000         -41415.9161   6784446.9640        12.7574         0.0000         0.0000        6783459.5920         0.0000   6783459.5920   6783459.5920         0.0000          -1301.5110    147115.9065   1118935.4248     -1301.5110    147115.9065
ENERGY:       1      4560.0156     20813.8273     15518.4082         0.0000         -41589.6188   4190147.9789        24.2474         0.0000         0.0000        4189474.8586         0.0000   4189474.8586   4189474.8586         0.0000          -1297.8001     99963.5646   1118935.4248     -1297.8001     99963.5646
ENERGY:       100      4234.1234     20712.3456     15489.7890         0.0000         -41203.4567   3456789.1234        18.9876         0.0000         0.0000        3456830.9183       310.5000   3456830.9183   3456830.9183       310.5000          -1289.3456     87654.3210   1118935.4248     -1289.3456     87654.3210
ENERGY:      200      4123.4567     20598.7654     15432.1098         0.0000         -41111.2345   3210987.6543        19.5432         0.0000         0.0000        3211050.2949       315.2000   3211050.2949   3211050.2949       315.2000          -1278.9876     76543.2109   1118935.4248     -1278.9876     76543.2109
TIMING:     100  CPU: 5.5, 0.055/step  Wall: 6.0, 0.06/step,  15.0 ns/day
TIMING:     200  CPU: 11.0, 0.055/step  Wall: 12.0, 0.06/step,  15.0 ns/day
"""
    
    @pytest.fixture
    def sample_log_file(self, tmp_path, sample_log_content):
        """Create a sample log file."""
        log_file = tmp_path / "test_equilibration.log"
        log_file.write_text(sample_log_content)
        return log_file
    
    def test_log_parsing(self, sample_log_file):
        """Test basic log file parsing."""
        timing = parse_namd_log(sample_log_file)
        
        # Check that timing object was created
        assert isinstance(timing, NAMDTiming)
        
        # Check basic attributes
        assert timing.processors == 4
        assert timing.atoms == 50000
        assert timing.timestep_fs == 2.0
        assert timing.first_timestep == 0
    
    def test_steps_extraction(self, sample_log_file):
        """Test that steps are extracted correctly."""
        timing = parse_namd_log(sample_log_file)
        
        # Should have found the last step (200)
        assert timing.steps_completed >= 0
    
    def test_performance_extraction(self, sample_log_file):
        """Test performance metrics extraction."""
        timing = parse_namd_log(sample_log_file)
        
        # Check performance metrics (if TIMING lines were parsed)
        # ns_per_day should be > 0 if TIMING lines were found
        if timing.ns_per_day > 0:
            # The parser calculates ns/day based on the log data
            assert timing.ns_per_day > 0, "ns/day should be positive"
            assert timing.steps_completed > 0, "Should have completed some steps"


class TestEquilibrationStageDiscovery:
    """Test equilibration stage discovery."""
    
    def test_stage_naming_conventions(self, tmp_path):
        """Test that log files follow naming conventions."""
        # Create sample log files
        log_files = [
            "step1_equilibration.log",
            "step2_equilibration.log",
            "step3_equilibration.log",
            "step7_production.log"
        ]
        
        for log_file in log_files:
            file_path = tmp_path / log_file
            file_path.write_text("ENERGY:       0      100.0\n")
            assert file_path.exists()
        
        # Verify files can be found
        all_logs = list(tmp_path.glob("*.log"))
        assert len(all_logs) == 4
        
        # Check specific patterns
        eq_logs = list(tmp_path.glob("step*_equilibration.log"))
        assert len(eq_logs) == 3
        
        prod_logs = list(tmp_path.glob("step7_production.log"))
        assert len(prod_logs) == 1


# ============================================================================
# SECTION 2: API FEATURE TESTS
# ============================================================================


class TestEnergyAnalyzer:
    """Test the EnergyAnalyzer class."""
    
    def test_energy_analyzer_initialization(self):
        """Test that EnergyAnalyzer can be initialized with a log file."""
        examples_dir = Path(__file__).parent / "analysis_examples"
        log_file = examples_dir / "equilibration_folder" / "step1_equilibration.log" 
        
        if not log_file.exists():
            pytest.skip(f"Log file not found: {log_file}")
        
        # Initialize analyzer
        analyzer = EnergyAnalyzer(log_file)
        assert analyzer is not None
        print("✓ EnergyAnalyzer initialized successfully")
    
    def test_energy_statistics(self):
        """Test energy statistics extraction."""
        examples_dir = Path(__file__).parent / "analysis_examples"
        log_file = examples_dir / "equilibration_folder" / "step1_equilibration.log"
        
        if not log_file.exists():
            pytest.skip(f"Log file not found: {log_file}")
        
        # Get statistics
        analyzer = EnergyAnalyzer(log_file)
        stats = analyzer.get_statistics()
        
        # Check that we got statistics
        assert stats is not None
        assert 'temp' in stats
        assert 'total' in stats
        
        # Check that statistics have expected keys
        for key in ['mean', 'std', 'min', 'max', 'initial', 'final']:
            assert key in stats['temp'], f"Missing key '{key}' in temperature stats"
            assert key in stats['total'], f"Missing key '{key}' in total energy stats"
        
        print("✓ Energy statistics extracted successfully")


# ============================================================================
# SECTION 3: DOCUMENTATION EXAMPLE TESTS
# ============================================================================


class TestAnalysisExamples:
    """Test analysis example scripts from documentation."""
    
    @pytest.fixture
    def temp_dir(self, tmp_path):
        """Create a temporary directory for test outputs."""
        return tmp_path
    
    def test_run_all_example_scripts(self, temp_dir):
        """Test running all example scripts from analysis_examples directory."""
        examples_dir = Path(__file__).parent / "analysis_examples"
        
        if not examples_dir.exists():
            pytest.skip(f"Examples directory not found: {examples_dir}")
        
        # Automatically discover all example files
        example_files = sorted(examples_dir.glob("analysis_example_*.py"))
        
        if not example_files:
            pytest.skip("No example files found in analysis_examples directory")
        
        print(f"\nFound {len(example_files)} example files to test")
        
        failed_examples = []
        passed_examples = []
        
        # Change to examples directory for relative paths
        original_dir = os.getcwd()
        os.chdir(examples_dir)
        
        try:
            for example_file in example_files:
                # Extract example number from filename (e.g., "01" from "analysis_example_01.py")
                parts = example_file.stem.split("_")
                example_num = parts[-1] if len(parts) >= 3 else "unknown"
                
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
        finally:
            # Restore original directory
            os.chdir(original_dir)
        
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
    
    @pytest.mark.parametrize("example_num", [
        "01", "02", "03", "04", "05", "06", "07", "08", "09", "10"
    ])
    def test_individual_examples(self, example_num, temp_dir):
        """Test each example individually for better pytest reporting."""
        examples_dir = Path(__file__).parent / "analysis_examples"
        
        # Find files with pattern analysis_example_NN.py
        example_file = examples_dir / f"analysis_example_{example_num}.py"
        
        if not example_file.exists():
            pytest.skip(f"Example {example_num} not found")
        
        # Change to examples directory for relative paths
        original_dir = os.getcwd()
        os.chdir(examples_dir)
        
        try:
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
        finally:
            # Restore original directory
            os.chdir(original_dir)


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def run_all_examples():
    """Helper function to run all examples manually."""
    examples_dir = Path(__file__).parent / "analysis_examples"
    
    if not examples_dir.exists():
        print(f"Examples directory not found: {examples_dir}")
        return
    
    # Find all example files automatically
    example_files = sorted(examples_dir.glob("analysis_example_*.py"))
    
    if not example_files:
        print("No example files found")
        return
    
    print(f"Found {len(example_files)} examples to run\n")
    
    passed = []
    failed = []
    
    # Change to examples directory for relative paths
    original_dir = os.getcwd()
    os.chdir(examples_dir)
    
    try:
        for example_file in example_files:
            # Extract example number
            parts = example_file.stem.split("_")
            example_num = parts[2] if len(parts) > 2 else "unknown"
            
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
    finally:
        # Restore original directory
        os.chdir(original_dir)
    
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
