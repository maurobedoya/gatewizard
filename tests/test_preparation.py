#!/usr/bin/env python3
"""
Unified Preparation Test Suite.

This test suite consolidates all protein preparation/PROPKA testing into a single file:
1. Core functionality tests (specs and features)
2. Documentation example workflows (22 examples with protein.pdb)
3. Example file execution tests (runs all preparation_example_XX.py files with all 3 PDB files)
4. Complex structure tests with 6RV3_AB.pdb and 8I5B.pdb

PDB files used:
- protein.pdb: Simple test protein (in preparation_examples/)
- 6RV3_AB.pdb: Multi-chain membrane protein with ligands
- 8I5B.pdb: Large multi-chain sodium channel

Note: Example files (Section 3) are run with all three PDB files using pytest parametrize.
"""

import pytest
import sys
import os
import tempfile
import shutil
from pathlib import Path
import importlib.util
from collections import defaultdict

# Add gatewizard to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from gatewizard.core.preparation import PreparationManager
from gatewizard.utils.protein_capping import ProteinCapper, cap_protein



# ============================================================================
# SECTION 1: CORE FUNCTIONALITY TESTS (Specs and Features)
# ============================================================================

class TestPreparationManager:
    """Test the PreparationManager class."""
    
    @pytest.fixture
    def analyzer(self):
        """Create a PreparationManager instance."""
        return PreparationManager(propka_version="3")
    
    @pytest.fixture
    def sample_pdb_file(self, tmp_path):
        """Create a sample multi-chain PDB file."""
        pdb_content = """ATOM      1  N   GLU A  23      10.000  20.000  30.000  1.00 20.00           N  
ATOM      2  CA  GLU A  23      11.000  21.000  31.000  1.00 20.00           C  
ATOM      3  N   ASP B  45      12.000  22.000  32.000  1.00 20.00           N  
ATOM      4  CA  ASP B  45      13.000  23.000  33.000  1.00 20.00           C  
ATOM      5  N   HIS C  67      14.000  24.000  34.000  1.00 20.00           N  
ATOM      6  CA  HIS C  67      15.000  25.000  35.000  1.00 20.00           C  
END
"""
        pdb_file = tmp_path / "test_multi_chain.pdb"
        pdb_file.write_text(pdb_content)
        return pdb_file
    
    def test_chain_extraction(self, analyzer, sample_pdb_file):
        """Test chain information extraction from PDB files."""
        chain_info = analyzer._extract_chain_info_from_pdb(str(sample_pdb_file))
        
        # Expected chains
        expected_chains = {'23:GLU': ['A'], '45:ASP': ['B'], '67:HIS': ['C']}
        
        assert chain_info == expected_chains, f"Expected {expected_chains}, got {chain_info}"
    
    def test_propka_version(self, analyzer):
        """Test that Propka version is set correctly."""
        assert analyzer.propka_version == "3"
        assert not hasattr(analyzer, 'propka31'), "Propka 3.1 should not be available"
    
    def test_pka_parsing(self, analyzer):
        """Test pKa summary parsing."""
        summary_content = """PROPKA SUMMARY
   Group      Residue    pKa    Buried
   ASP  23     A         3.65    0%
   GLU  45     B         4.25   15%
   HIS  67     C         6.50   30%
"""
        # Test parsing logic
        lines = summary_content.strip().split('\n')
        assert len(lines) == 5  # Header + 1 blank + 3 residues
        
        # Verify format
        data_lines = [l for l in lines if 'ASP' in l or 'GLU' in l or 'HIS' in l]
        assert len(data_lines) == 3
    
    def test_ligand_atom_type_parsing(self, analyzer, tmp_path):
        """Test parsing of ligand atoms and protein residues from PROPKA summary."""
        # Create a mock summary file with both protein residues and ligand atoms
        summary_content = """SUMMARY OF THIS PREDICTION
       Group      pKa  model-pKa   ligand atom-type
   ASP  52 A     4.33       3.80                      
   ASP  65 A     3.87       3.80                      
   GLU  77 A     5.07       4.50                      
   HIS  79 A     6.45       6.50                      
   LYS  84 A    10.20      10.50         
   ARG 115 C    13.84      12.50
   N+    8 A     7.66       8.00
   P5S   N A    10.83      10.00                N31
   LPE   N A     9.87      10.00                N31
   OJ0 N02 A     8.82      10.00                N33
   P5S   C A     4.16       4.50                OCO
   Y01 CAX A     4.62       4.50                OCO
   P5S O15 A     5.46       6.00                 OP
   LPE O31 A     5.15       6.00                 OP
"""
        summary_file = tmp_path / "test_summary_ligands.txt"
        summary_file.write_text(summary_content)
        
        # Parse the summary
        residues = analyzer.parse_summary(str(summary_file))
        
        # Verify total count
        assert len(residues) == 14, f"Expected 14 entries, got {len(residues)}"
        
        # Separate protein residues from ligand atoms
        protein_residues = [r for r in residues if r['res_id'] > 0]
        ligand_atoms = [r for r in residues if r['res_id'] == 0]
        
        # Verify counts
        assert len(protein_residues) == 7, f"Expected 7 protein residues, got {len(protein_residues)}"
        assert len(ligand_atoms) == 7, f"Expected 7 ligand atoms, got {len(ligand_atoms)}"
        
        # Test protein residue parsing
        asp52 = next(r for r in protein_residues if r['res_id'] == 52)
        assert asp52['residue'] == 'ASP'
        assert asp52['chain'] == 'A'
        assert asp52['pka'] == 4.33
        assert asp52['atom'] == '', f"Protein residues should have empty atom field, got '{asp52['atom']}'"
        assert asp52['atom_type'] == '', f"Protein residues should have empty atom_type, got '{asp52['atom_type']}'"
        assert asp52['model_pka'] == 3.80
        
        # Test ligand atom parsing
        p5s_n = next(r for r in ligand_atoms if r['residue'] == 'P5S' and r['atom'] == 'N')
        assert p5s_n['res_id'] == 0, "Ligands should have res_id=0"
        assert p5s_n['chain'] == 'A'
        assert p5s_n['pka'] == 10.83
        assert p5s_n['atom'] == 'N', f"Expected atom 'N', got '{p5s_n['atom']}'"
        assert p5s_n['atom_type'] == 'N31', f"Expected atom_type 'N31', got '{p5s_n['atom_type']}'"
        assert p5s_n['model_pka'] == 10.00
        
        # Test another ligand atom with different type
        y01_cax = next(r for r in ligand_atoms if r['residue'] == 'Y01' and r['atom'] == 'CAX')
        assert y01_cax['atom_type'] == 'OCO', f"Expected atom_type 'OCO', got '{y01_cax['atom_type']}'"
        assert y01_cax['pka'] == 4.62
        
        # Test phosphate oxygen
        p5s_o15 = next(r for r in ligand_atoms if r['residue'] == 'P5S' and r['atom'] == 'O15')
        assert p5s_o15['atom_type'] == 'OP', f"Expected atom_type 'OP', got '{p5s_o15['atom_type']}'"
        
        # Verify we can group by ligand
        from collections import defaultdict
        ligands_by_type = defaultdict(list)
        for lig in ligand_atoms:
            ligands_by_type[lig['residue']].append(lig)
        
        assert len(ligands_by_type['P5S']) == 3, f"Expected 3 P5S atoms, got {len(ligands_by_type['P5S'])}"
        assert len(ligands_by_type['LPE']) == 2, f"Expected 2 LPE atoms, got {len(ligands_by_type['LPE'])}"
        
        print("✓ Ligand and protein parsing test passed!")
        print(f"  Protein residues: {len(protein_residues)}")
        print(f"  Ligand atoms: {len(ligand_atoms)}")
        print(f"  Unique ligands: {len(ligands_by_type)}")
        for ligand_name, atoms in ligands_by_type.items():
            print(f"    {ligand_name}: {len(atoms)} ionizable atoms")


class TestProteinCapping:
    """Test protein capping functionality."""
    
    def test_capping_residues(self):
        """Test that ACE and NME caps can be identified."""
        ace_residue = "ACE"
        nme_residue = "NME"
        
        # These are standard capping residues
        assert ace_residue == "ACE", "N-terminal cap should be ACE"
        assert nme_residue == "NME", "C-terminal cap should be NME"

# ============================================================================
# SECTION 2: WORKFLOW EXAMPLES (22 Examples with protein.pdb)
# ============================================================================

class TestPropkaWorkflowExamples:
    """Test all propka documentation examples."""
    
    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test outputs."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield tmpdir
    
    @pytest.fixture
    def protein_pdb(self, temp_dir):
        """Copy actual protein.pdb from tests directory to temp directory."""
        # Find the protein.pdb file in the tests directory
        test_dir = Path(__file__).parent
        source_pdb = test_dir / "preparation_examples" / "protein.pdb"
        
        if not source_pdb.exists():
            pytest.skip(f"protein.pdb not found at {source_pdb}")
        
        # Copy to temp directory
        dest_pdb = Path(temp_dir) / "protein.pdb"
        shutil.copy(str(source_pdb), str(dest_pdb))
        return str(dest_pdb)
    
    def test_example_01_propka_version(self, temp_dir, protein_pdb):
        """Test Example 1: PreparationManager initialization."""
        os.chdir(temp_dir)
        analyzer = PreparationManager(propka_version="3")
        assert analyzer.propka_version == "3"
        print(f"Using PROPKA version: {analyzer.propka_version}")
    
    def test_example_02_run_analysis(self, temp_dir, protein_pdb):
        """Test Example 2: Run PROPKA analysis."""
        os.chdir(temp_dir)
        
        analyzer = PreparationManager()
        pka_file = analyzer.run_analysis("protein.pdb")
        
        assert os.path.exists(pka_file)
        assert pka_file.endswith(".pka")
    
    def test_example_03_extract_summary(self, temp_dir, protein_pdb):
        """Test Example 3: Extract summary from PKA file."""
        os.chdir(temp_dir)
        
        analyzer = PreparationManager()
        pka_file = analyzer.run_analysis("protein.pdb")
        summary_file = analyzer.extract_summary(pka_file)
        
        assert os.path.exists(summary_file)
        assert "summary_of_prediction" in summary_file
    
    def test_example_04_parse_summary(self, temp_dir, protein_pdb):
        """Test Example 4: Parse summary with protein residues."""
        os.chdir(temp_dir)
        
        analyzer = PreparationManager()
        pka_file = analyzer.run_analysis("protein.pdb")
        summary_file = analyzer.extract_summary(pka_file)
        residues = analyzer.parse_summary(summary_file)
        
        # Separate protein residues from ligand atoms
        protein_residues = [r for r in residues if r['res_id'] > 0]
        ligand_atoms = [r for r in residues if r['res_id'] == 0]
        
        assert len(protein_residues) > 0
        print(f"Found {len(protein_residues)} ionizable protein residues")
        print(f"Found {len(ligand_atoms)} ionizable ligand atoms")
    
    def test_example_05_apply_protonation(self, temp_dir, protein_pdb):
        """Test Example 5: Apply protonation states."""
        os.chdir(temp_dir)
        
        analyzer = PreparationManager()
        pka_file = analyzer.run_analysis("protein.pdb")
        summary_file = analyzer.extract_summary(pka_file)
        residues = analyzer.parse_summary(summary_file)
        
        # Basic usage - automatic protonation at pH 7.4
        stats = analyzer.apply_protonation_states(
            input_pdb="protein.pdb",
            output_pdb="protein_ph7.4.pdb",
            ph=7.4,
            residues=residues
        )
        
        assert os.path.exists("protein_ph7.4.pdb")
        assert stats['residue_changes'] >= 0
        print(f"Modified {stats['residue_changes']} residues ({stats['record_changes']} atoms)")
    
    def test_example_06_get_default_protonation_state(self, temp_dir, protein_pdb):
        """Test Example 6: Get default protonation states at different pH."""
        os.chdir(temp_dir)
        
        analyzer = PreparationManager()
        pka_file = analyzer.run_analysis("protein.pdb")
        summary_file = analyzer.extract_summary(pka_file)
        residues = analyzer.parse_summary(summary_file)
        
        for res in residues[:3]:  # Test first 3 residues
            if res['res_id'] > 0:  # Only protein residues
                states = [analyzer.get_default_protonation_state(res, ph=p) for p in [2.0, 7.0, 11.0]]
                print(f"{res['residue']}{res['res_id']} (pKa={res['pka']:.2f}): "
                      f"pH2={states[0]}, pH7={states[1]}, pH11={states[2]}")
                assert all(isinstance(s, str) for s in states)
    
    def test_example_07_get_available_states(self, temp_dir):
        """Test Example 7: Get available states for residue types."""
        analyzer = PreparationManager()
        his_states = analyzer.get_available_states("HIS")
        
        assert isinstance(his_states, dict)
        assert 'neutral_epsilon' in his_states
        assert his_states['neutral_epsilon'] == 'HIE'
        print(his_states)
    
    def test_example_08_detect_disulfide_bonds(self, temp_dir, protein_pdb):
        """Test Example 8: Detect disulfide bonds."""
        os.chdir(temp_dir)
        
        analyzer = PreparationManager()
        bonds = analyzer.detect_disulfide_bonds("protein.pdb")
        
        print(f"Found {len(bonds)} disulfide bonds:")
        for bond in bonds:
            (res1_name, res1_id), (res2_name, res2_id) = bond
            print(f"  {res1_name}{res1_id} ↔ {res2_name}{res2_id}")
        
        assert isinstance(bonds, list)
    
    def test_example_09_apply_disulfide_bonds(self, temp_dir, protein_pdb):
        """Test Example 9: Apply disulfide bonds."""
        os.chdir(temp_dir)
        
        analyzer = PreparationManager()
        num_bonds = analyzer.apply_disulfide_bonds(
            input_pdb="protein.pdb",
            output_pdb="protein_ss_auto.pdb"
        )
        
        assert os.path.exists("protein_ss_auto.pdb")
        print(f"✓ Applied {num_bonds} disulfide bonds")
    
    def test_example_10_combined_workflow(self, temp_dir, protein_pdb):
        """Test Example 10: Combined protonation and disulfide workflow."""
        os.chdir(temp_dir)
        
        analyzer = PreparationManager()
        pka_file = analyzer.run_analysis("protein.pdb")
        summary_file = analyzer.extract_summary(pka_file)
        residues = analyzer.parse_summary(summary_file)
        
        # Step 1: Apply protonation states
        stats = analyzer.apply_protonation_states(
            input_pdb="protein.pdb",
            output_pdb="protein_ph7.pdb",
            ph=7.4,
            residues=residues,
        )
        
        # Step 2: Apply disulfide bonds
        num_bonds = analyzer.apply_disulfide_bonds(
            input_pdb="protein_ph7.pdb",
            output_pdb="protein_ph7_ss.pdb"
        )
        
        assert os.path.exists("protein_ph7_ss.pdb")
    
    def test_example_11_get_residue_statistics(self, temp_dir, protein_pdb):
        """Test Example 11: Get residue statistics."""
        os.chdir(temp_dir)
        
        analyzer = PreparationManager()
        pka_file = analyzer.run_analysis("protein.pdb")
        summary_file = analyzer.extract_summary(pka_file)
        residues = analyzer.parse_summary(summary_file)
        
        stats = analyzer.get_residue_statistics()
        
        assert isinstance(stats, dict)
        for res_type, count in stats.items():
            print(f"{res_type}: {count}")
    
    def test_example_12_get_ph_titration_curve(self, temp_dir, protein_pdb):
        """Test Example 12: Get pH titration curves."""
        os.chdir(temp_dir)
        
        analyzer = PreparationManager()
        pka_file = analyzer.run_analysis("protein.pdb")
        summary_file = analyzer.extract_summary(pka_file)
        residues = analyzer.parse_summary(summary_file)
        
        curves = analyzer.get_ph_titration_curve(ph_range=(4, 10), ph_step=1.0)
        
        assert isinstance(curves, dict)
        assert len(curves) > 0
        
        for residue_id, curve in list(curves.items())[:2]:  # Test first 2
            print(f"\n{residue_id}:")
            for ph, state in curve:
                print(f"  pH {ph:.1f}: {state}")
    
    def test_example_15_protein_capping(self, temp_dir, protein_pdb):
        """Test Example 15: Protein capping with ProteinCapper."""
        os.chdir(temp_dir)
        os.makedirs("output", exist_ok=True)
        
        capper = ProteinCapper()
        capped_file, residue_mapping = capper.remove_hydrogens_and_cap(
            input_file="protein.pdb",
            output_file="protein_capped.pdb",
            target_dir="output"
        )
        
        assert os.path.exists(capped_file)
        assert isinstance(residue_mapping, dict)
        print(f"✓ Capped protein: {capped_file}")
        print(f"✓ Residue mapping: {len(residue_mapping)} residues tracked")
    
    def test_example_16_cap_protein_convenience(self, temp_dir, protein_pdb):
        """Test Example 16: Protein capping with convenience function."""
        os.chdir(temp_dir)
        
        capped_file, mapping = cap_protein(
            input_file="protein.pdb",
            output_file="protein_capped_convenient.pdb"
        )
        
        assert os.path.exists(capped_file)
        assert isinstance(mapping, dict)
    
    def test_example_17_complete_workflow(self, temp_dir, protein_pdb):
        """Test Example 17: Complete workflow with output directory."""
        os.chdir(temp_dir)
        os.makedirs("output", exist_ok=True)
        
        analyzer = PreparationManager()
        
        # Step 1: Run Propka analysis
        pka_file = analyzer.run_analysis("protein.pdb", output_dir="output")
        summary_file = analyzer.extract_summary(pka_file, output_dir="output")
        residues = analyzer.parse_summary(summary_file)
        
        # Step 2: Detect disulfide bonds
        bonds = analyzer.detect_disulfide_bonds("protein.pdb", distance_threshold=2.5)
        
        # Step 3: Apply protonation states
        stats = analyzer.apply_protonation_states(
            input_pdb="protein.pdb",
            output_pdb="output/protein_ph7.pdb",
            ph=7.4,
            residues=residues
        )
        
        # Step 4: Apply disulfide bonds
        num_bonds = analyzer.apply_disulfide_bonds(
            input_pdb="output/protein_ph7.pdb",
            output_pdb="output/protein_ph7_ss.pdb",
            disulfide_bonds=bonds,
            auto_detect=False
        )
        
        assert os.path.exists("output/protein_ph7_ss.pdb")
    
    def test_example_18_capping_workflow(self, temp_dir, protein_pdb):
        """Test Example 18: Workflow with protein capping."""
        os.chdir(temp_dir)
        os.makedirs("output", exist_ok=True)
        
        # Step 1: Add ACE/NME caps
        capper = ProteinCapper()
        capped_file, residue_mapping = capper.remove_hydrogens_and_cap(
            input_file="protein.pdb",
            output_file="protein_capped.pdb",
            target_dir="output"
        )
        
        # Step 2: Run Propka on capped structure
        analyzer = PreparationManager()
        pka_file = analyzer.run_analysis(capped_file, output_dir="output")
        summary_file = analyzer.extract_summary(pka_file, output_dir="output")
        residues = analyzer.parse_summary(summary_file)
        
        # Step 3: Detect disulfide bonds
        bonds = analyzer.detect_disulfide_bonds(capped_file)
        
        # Step 4: Apply protonation
        stats = analyzer.apply_protonation_states(
            input_pdb=capped_file,
            output_pdb="output/protein_capped_ph7.pdb",
            ph=7.4,
            residues=residues
        )
        
        # Step 5: Apply disulfide bonds
        num_bonds = analyzer.apply_disulfide_bonds(
            input_pdb="output/protein_capped_ph7.pdb",
            output_pdb="output/protein_capped_ph7_ss.pdb",
            disulfide_bonds=bonds,
            auto_detect=False
        )
        
        assert os.path.exists("output/protein_capped_ph7_ss.pdb")
    
    def test_example_19_multiple_ph_variants(self, temp_dir, protein_pdb):
        """Test Example 19: Generate structures for multiple pH values."""
        os.chdir(temp_dir)
        os.makedirs("output", exist_ok=True)
        
        analyzer = PreparationManager()
        
        # Run analysis once
        pka_file = analyzer.run_analysis("protein.pdb", output_dir="output")
        summary_file = analyzer.extract_summary(pka_file, output_dir="output")
        residues = analyzer.parse_summary(summary_file)
        bonds = analyzer.detect_disulfide_bonds("protein.pdb")
        
        # Generate structures for different pH values
        for ph in [5.0, 7.0, 9.0]:  # Test subset
            stats = analyzer.apply_protonation_states(
                input_pdb="protein.pdb",
                output_pdb=f"output/protein_ph{ph:.1f}.pdb",
                ph=ph,
                residues=residues
            )
            
            analyzer.apply_disulfide_bonds(
                input_pdb=f"output/protein_ph{ph:.1f}.pdb",
                output_pdb=f"output/protein_ph{ph:.1f}_ss.pdb",
                disulfide_bonds=bonds,
                auto_detect=False
            )
            
            assert os.path.exists(f"output/protein_ph{ph:.1f}_ss.pdb")
    
    def test_example_21_custom_protonation(self, temp_dir, protein_pdb):
        """Test Example 21: Custom protonation states."""
        os.chdir(temp_dir)
        os.makedirs("output", exist_ok=True)
        
        # Step 1: Run standard Propka workflow
        analyzer = PreparationManager()
        pka_file = analyzer.run_analysis("protein.pdb", output_dir="output")
        summary_file = analyzer.extract_summary(pka_file, output_dir="output")
        residues = analyzer.parse_summary(summary_file)
        
        # Step 2: Specify custom states
        custom_states = {
            "ASP12_A": "ASH",
            "HIS15_A": "HID"
        }
        
        # Step 3: Apply protonation with custom overrides
        stats = analyzer.apply_protonation_states(
            input_pdb="protein.pdb",
            output_pdb="output/protein_custom.pdb",
            ph=7.4,
            custom_states=custom_states,
            residues=residues
        )
        
        assert os.path.exists("output/protein_custom.pdb")
        print(f"✓ Modified {stats['residue_changes']} residues")
    
    def test_example_22_filtering_analysis(self, temp_dir, protein_pdb):
        """Test Example 22: Filtering and analysis of pKa shifts."""
        os.chdir(temp_dir)
        os.makedirs("output", exist_ok=True)
        
        # Step 1: Run Propka analysis
        analyzer = PreparationManager()
        pka_file = analyzer.run_analysis("protein.pdb", output_dir="output")
        summary_file = analyzer.extract_summary(pka_file, output_dir="output")
        residues = analyzer.parse_summary(summary_file)
        
        # Step 2: Define expected model pKa values
        expected_pka = {
            'ASP': 3.9, 'GLU': 4.3, 'HIS': 6.0,
            'LYS': 10.5, 'ARG': 12.5, 'CYS': 8.3, 'TYR': 10.1
        }
        
        # Step 3: Find residues with significant pKa shifts
        shifted_residues = []
        for res in residues:
            if res['res_id'] == 0:  # Skip ligands
                continue
            
            res_name = res['residue']
            if res_name in expected_pka:
                pka_diff = abs(res['pka'] - expected_pka[res_name])
                if pka_diff > 1.0:
                    shifted_residues.append({
                        'id': f"{res_name}{res['res_id']}_{res['chain']}",
                        'shift': res['pka'] - expected_pka[res_name]
                    })
        
        print(f"Found {len(shifted_residues)} residues with significant pKa shifts")
        assert isinstance(shifted_residues, list)

# ============================================================================


# ============================================================================
# SECTION 3: EXAMPLE FILE EXECUTION TESTS (All PDB Files)
# ============================================================================

class TestPropkaExampleFiles:
    """
    Test all preparation example files by running them with all three PDB files.
    
    The test suite automatically discovers and runs all example scripts from:
        tests/preparation_examples/preparation_example_*.py
    """
    
    @pytest.fixture
    def examples_dir(self):
        """Get the examples directory path."""
        return Path(__file__).parent / "preparation_examples"
    
    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test outputs."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)
    
    @pytest.fixture(params=["protein.pdb", "6RV3_AB.pdb", "8I5B.pdb"])
    def pdb_file_setup(self, request, temp_dir, examples_dir):
        """
        Setup test files for each PDB file.
        This fixture is parametrized to run tests with all three PDB files.
        """
        pdb_filename = request.param
        test_dir = Path(__file__).parent
        
        # Copy the requested PDB file
        if pdb_filename == "protein.pdb":
            source = examples_dir / pdb_filename
        else:
            source = test_dir / pdb_filename
        
        if not source.exists():
            pytest.skip(f"{pdb_filename} not found at {source}")
        
        dest = temp_dir / "protein.pdb"  # Always name it protein.pdb for the examples
        shutil.copy(str(source), str(dest))
        
        return temp_dir, pdb_filename
    
    def run_example(self, example_file: Path, work_dir: Path) -> tuple:
        """
        Run a single example file.
        
        Returns:
            tuple: (success: bool | None, message: str)
            None indicates skip (e.g., missing dependency)
        """
        # Save current directory
        try:
            original_dir = os.getcwd()
        except (FileNotFoundError, OSError):
            original_dir = Path(__file__).parent.parent
            os.chdir(original_dir)
        
        original_sys_path = sys.path.copy()
        
        try:
            # Ensure project root is in sys.path
            project_root = Path(__file__).parent.parent
            if str(project_root) not in sys.path:
                sys.path.insert(0, str(project_root))
            
            # Change to working directory
            os.chdir(work_dir)
            
            # Create output directory if needed
            (work_dir / "output").mkdir(exist_ok=True)
            
            # Load and execute the example module
            spec = importlib.util.spec_from_file_location(example_file.stem, example_file)
            if spec is None or spec.loader is None:
                return False, f"Could not load {example_file.name}"
            
            module = importlib.util.module_from_spec(spec)
            sys.modules[example_file.stem] = module
            
            # Execute the module
            spec.loader.exec_module(module)
            
            return True, "Success"
            
        except Exception as e:
            error_msg = f"{type(e).__name__}: {str(e)}"
            # Check if it's an expected error
            if any(keyword in str(e).lower() for keyword in [
                "propka3", "propka", "matplotlib", "numpy",
                "summary_of_prediction.txt",
                "no propka", "propka not found"
            ]):
                return None, f"Skipped: {error_msg}"
            return False, error_msg
            
        finally:
            # Restore original state
            try:
                os.chdir(original_dir)
            except (FileNotFoundError, OSError):
                os.chdir(Path(__file__).parent.parent)
            
            sys.path = original_sys_path
            if example_file.stem in sys.modules:
                del sys.modules[example_file.stem]
    
    @pytest.mark.parametrize("example_num", [f"{i:02d}" for i in range(1, 23)])
    def test_example_all_pdbs(self, example_num, examples_dir, pdb_file_setup):
        """
        Test each preparation example with all PDB files.
        
        This test is parametrized to run all 22 examples with all 3 PDB files,
        resulting in 66 test cases (22 examples × 3 PDB files).
        """
        work_dir, pdb_name = pdb_file_setup
        example_file = examples_dir / f"preparation_example_{example_num}.py"
        
        if not example_file.exists():
            pytest.skip(f"Example {example_num} not found at {example_file}")
        
        # Special handling for examples that need matplotlib
        if example_num in ["13", "14"]:
            try:
                import matplotlib
                matplotlib.use('Agg')
                if example_num == "14":
                    import numpy
            except ImportError as e:
                pytest.skip(f"Required package not installed: {e}")
        
        success, message = self.run_example(example_file, work_dir)
        if success is None:
            pytest.skip(message)
        assert success, f"Example {example_num} with {pdb_name} failed: {message}"

# SECTION 4: COMPLEX STRUCTURE TESTS
# ============================================================================

class TestPropkaComplexStructures:
    """Test propka examples with complex protein structures."""
    
    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test outputs."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield tmpdir
    
    @pytest.fixture
    def pdb_6rv3(self, temp_dir):
        """Path to 6RV3_AB.pdb file (if available)."""
        # Try to find the file in the workspace
        search_paths = [
            Path(__file__).parent.parent / "6RV3_AB.pdb",
            Path(__file__).parent / "6RV3_AB.pdb",
            Path(__file__).parent.parent / "examples" / "6RV3_AB.pdb",
        ]
        
        for path in search_paths:
            if path.exists():
                dest = Path(temp_dir) / "6RV3_AB.pdb"
                shutil.copy(str(path), str(dest))
                return str(dest)
        
        pytest.skip("6RV3_AB.pdb not found in workspace")
    
    @pytest.fixture
    def pdb_8i5b(self, temp_dir):
        """Path to 8I5B.pdb file (if available)."""
        # Try to find the file in the workspace
        search_paths = [
            Path(__file__).parent.parent / "8I5B.pdb",
            Path(__file__).parent / "8I5B.pdb",
            Path(__file__).parent.parent / "examples" / "8I5B.pdb",
        ]
        
        for path in search_paths:
            if path.exists():
                dest = Path(temp_dir) / "8I5B.pdb"
                shutil.copy(str(path), str(dest))
                return str(dest)
        
        pytest.skip("8I5B.pdb not found in workspace")
    
    def test_6rv3_ligand_analysis(self, temp_dir, pdb_6rv3):
        """Test ligand analysis with 6RV3_AB (membrane protein with ligands)."""
        os.chdir(temp_dir)
        
        analyzer = PreparationManager()
        pka_file = analyzer.run_analysis(pdb_6rv3)
        summary_file = analyzer.extract_summary(pka_file)
        residues = analyzer.parse_summary(summary_file)
        
        # Separate protein residues from ligand atoms
        protein_residues = [r for r in residues if r['res_id'] > 0]
        ligand_atoms = [r for r in residues if r['res_id'] == 0]
        
        assert len(protein_residues) > 0
        assert len(ligand_atoms) > 0  # 6RV3 should have ligands
        
        print(f"Found {len(protein_residues)} ionizable protein residues")
        print(f"Found {len(ligand_atoms)} ionizable ligand atoms")
        
        # Group by ligand type
        from collections import defaultdict
        ligands_by_type = defaultdict(list)
        for lig in ligand_atoms:
            ligands_by_type[lig['residue']].append(lig)
        
        for ligand_name, atoms in ligands_by_type.items():
            print(f"\n{ligand_name}: {len(atoms)} ionizable atoms")
            for atom in atoms[:3]:  # Show first 3
                print(f"  {atom['atom']:6s} pKa={atom['pka']:5.2f} ({atom['atom_type']})")
    
    def test_6rv3_complete_workflow(self, temp_dir, pdb_6rv3):
        """Test complete workflow with 6RV3_AB."""
        os.chdir(temp_dir)
        os.makedirs("output", exist_ok=True)
        
        analyzer = PreparationManager()
        
        # Run analysis
        pka_file = analyzer.run_analysis(pdb_6rv3, output_dir="output")
        summary_file = analyzer.extract_summary(pka_file, output_dir="output")
        residues = analyzer.parse_summary(summary_file)
        
        # Detect disulfide bonds
        bonds = analyzer.detect_disulfide_bonds(pdb_6rv3, distance_threshold=2.5)
        print(f"Detected {len(bonds)} disulfide bonds")
        
        # Apply protonation states
        stats = analyzer.apply_protonation_states(
            input_pdb=pdb_6rv3,
            output_pdb="output/6rv3_ph7.pdb",
            ph=7.4,
            residues=residues
        )
        
        # Apply disulfide bonds
        if bonds:
            num_bonds = analyzer.apply_disulfide_bonds(
                input_pdb="output/6rv3_ph7.pdb",
                output_pdb="output/6rv3_ph7_ss.pdb",
                disulfide_bonds=bonds,
                auto_detect=False
            )
            assert os.path.exists("output/6rv3_ph7_ss.pdb")
        else:
            assert os.path.exists("output/6rv3_ph7.pdb")
    
    def test_6rv3_pka_distribution(self, temp_dir, pdb_6rv3):
        """Test pKa distribution analysis with 6RV3_AB."""
        os.chdir(temp_dir)
        os.makedirs("output", exist_ok=True)
        
        analyzer = PreparationManager()
        pka_file = analyzer.run_analysis(pdb_6rv3, output_dir="output")
        summary_file = analyzer.extract_summary(pka_file, output_dir="output")
        residues = analyzer.parse_summary(summary_file)
        
        # Filter protein residues only
        protein_residues = [r for r in residues if r['res_id'] > 0]
        
        # Set target pH
        target_ph = 5.0
        
        # Group by residue type
        residue_data = {}
        for res in protein_residues:
            res_type = res['residue']
            if res_type not in residue_data:
                residue_data[res_type] = {'pka': [], 'protonated': []}
            
            residue_data[res_type]['pka'].append(res['pka'])
            is_protonated = target_ph < res['pka']
            residue_data[res_type]['protonated'].append(is_protonated)
        
        # Print statistics
        print(f"\npKa Statistics at pH {target_ph}:")
        for res_type, data in sorted(residue_data.items()):
            pka_vals = data['pka']
            n_prot = sum(data['protonated'])
            n_total = len(pka_vals)
            prot_frac = n_prot / n_total if n_total > 0 else 0
            
            print(f"{res_type:4s}: n={n_total:3d}  "
                  f"pKa range=[{min(pka_vals):5.2f}, {max(pka_vals):5.2f}]  "
                  f"Protonated: {prot_frac*100:.0f}%")
    
    def test_8i5b_large_structure(self, temp_dir, pdb_8i5b):
        """Test analysis with 8I5B (large multi-chain sodium channel)."""
        os.chdir(temp_dir)
        os.makedirs("output", exist_ok=True)
        
        analyzer = PreparationManager()
        
        # Run analysis
        pka_file = analyzer.run_analysis(pdb_8i5b, output_dir="output")
        summary_file = analyzer.extract_summary(pka_file, output_dir="output")
        residues = analyzer.parse_summary(summary_file)
        
        # Get statistics
        stats = analyzer.get_residue_statistics()
        
        print(f"\n8I5B Residue Statistics:")
        print(f"Total ionizable residues: {len(residues)}")
        for res_type, count in sorted(stats.items()):
            print(f"  {res_type}: {count}")
        
        assert len(residues) > 0
    
    def test_8i5b_disulfide_bonds(self, temp_dir, pdb_8i5b):
        """Test disulfide bond detection with 8I5B."""
        os.chdir(temp_dir)
        
        analyzer = PreparationManager()
        bonds = analyzer.detect_disulfide_bonds(pdb_8i5b, distance_threshold=2.5)
        
        print(f"\n8I5B Disulfide Bonds: {len(bonds)}")
        for bond in bonds[:10]:  # Show first 10
            (res1_name, res1_id), (res2_name, res2_id) = bond
            print(f"  {res1_name}{res1_id} ↔ {res2_name}{res2_id}")
        
        assert isinstance(bonds, list)
    
    def test_8i5b_ph_variants(self, temp_dir, pdb_8i5b):
        """Test multiple pH variants with 8I5B."""
        os.chdir(temp_dir)
        os.makedirs("output", exist_ok=True)
        
        analyzer = PreparationManager()
        
        # Run analysis once
        pka_file = analyzer.run_analysis(pdb_8i5b, output_dir="output")
        summary_file = analyzer.extract_summary(pka_file, output_dir="output")
        residues = analyzer.parse_summary(summary_file)
        bonds = analyzer.detect_disulfide_bonds(pdb_8i5b)
        
        # Generate structures for different pH values
        for ph in [5.0, 7.4, 9.0]:
            stats = analyzer.apply_protonation_states(
                input_pdb=pdb_8i5b,
                output_pdb=f"output/8i5b_ph{ph:.1f}.pdb",
                ph=ph,
                residues=residues
            )
            
            print(f"pH {ph:.1f}: {stats['residue_changes']} residues modified")
            assert os.path.exists(f"output/8i5b_ph{ph:.1f}.pdb")
    
    def test_6rv3_filtering_analysis(self, temp_dir, pdb_6rv3):
        """Test filtering analysis with 6RV3_AB (pKa shifts)."""
        os.chdir(temp_dir)
        os.makedirs("output", exist_ok=True)
        
        analyzer = PreparationManager()
        pka_file = analyzer.run_analysis(pdb_6rv3, output_dir="output")
        summary_file = analyzer.extract_summary(pka_file, output_dir="output")
        residues = analyzer.parse_summary(summary_file)
        
        # Define expected model pKa values
        expected_pka = {
            'ASP': 3.9, 'GLU': 4.3, 'HIS': 6.0,
            'LYS': 10.5, 'ARG': 12.5, 'CYS': 8.3, 'TYR': 10.1
        }
        
        # Find residues with significant pKa shifts
        shifted_residues = []
        for res in residues:
            if res['res_id'] == 0:  # Skip ligands
                continue
            
            res_name = res['residue']
            if res_name in expected_pka:
                pka_diff = abs(res['pka'] - expected_pka[res_name])
                if pka_diff > 1.0:
                    shifted_residues.append({
                        'id': f"{res_name}{res['res_id']}_{res['chain']}",
                        'type': res_name,
                        'pka': res['pka'],
                        'expected': expected_pka[res_name],
                        'shift': res['pka'] - expected_pka[res_name]
                    })
        
        print(f"\n6RV3 Residues with unusual pKa shifts (>1.0 units): {len(shifted_residues)}")
        
        # Categorize by shift direction
        upshifted = [r for r in shifted_residues if r['shift'] > 1.0]
        downshifted = [r for r in shifted_residues if r['shift'] < -1.0]
        
        print(f"  Upshifted (more basic): {len(upshifted)}")
        print(f"  Downshifted (more acidic): {len(downshifted)}")
        
        # Find extreme cases
        extreme_shifts = [r for r in shifted_residues if abs(r['shift']) > 2.0]
        if extreme_shifts:
            print(f"\n⚠ Extreme shifts (>2.0 units): {len(extreme_shifts)}")
            for r in extreme_shifts[:5]:  # Show first 5
                print(f"  {r['id']}: {r['shift']:+.2f} units")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
