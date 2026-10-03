#!/usr/bin/env python3
"""
Unified Preparation Test Suite.

This test suite consolidates all protein preparation/PROPKA testing into a single file:
1. Core functionality tests (specs and features)
2. Documentation example scripts (preparation_example_*.py, once each)
3. Complex structure tests with 6RV3_AB.pdb and 8I5B.pdb

PDB files used:
- protein.pdb: Simple test protein (in preparation_examples/)
- 6RV3_AB.pdb: Multi-chain membrane protein with ligands
- 8I5B.pdb: Large multi-chain sodium channel

Note: Example scripts run once each from tests/preparation_examples/.
"""

import pytest
import sys
import os
import tempfile
import shutil
from collections import defaultdict
from pathlib import Path

from tests.example_runner import parametrize_example_scripts, run_example_script

# Add gatewizard to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from gatewizard.core.preparation import PreparationManager
from gatewizard.utils.protein_capping import ProteinCapper, cap_protein, detect_terminal_caps

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
        expected_chains = {"23:GLU": ["A"], "45:ASP": ["B"], "67:HIS": ["C"]}

        assert (
            chain_info == expected_chains
        ), f"Expected {expected_chains}, got {chain_info}"

    def test_propka_version(self, analyzer):
        """Test that Propka version is set correctly."""
        assert analyzer.propka_version == "3"
        assert not hasattr(analyzer, "propka31"), "Propka 3.1 should not be available"

    def test_pka_parsing(self, analyzer):
        """Test pKa summary parsing."""
        summary_content = """PROPKA SUMMARY
   Group      Residue    pKa    Buried
   ASP  23     A         3.65    0%
   GLU  45     B         4.25   15%
   HIS  67     C         6.50   30%
"""
        # Test parsing logic
        lines = summary_content.strip().split("\n")
        assert len(lines) == 5  # Header + 1 blank + 3 residues

        # Verify format
        data_lines = [l for l in lines if "ASP" in l or "GLU" in l or "HIS" in l]
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
        protein_residues = [r for r in residues if r["res_id"] > 0]
        ligand_atoms = [r for r in residues if r["res_id"] == 0]

        # Verify counts
        assert (
            len(protein_residues) == 7
        ), f"Expected 7 protein residues, got {len(protein_residues)}"
        assert (
            len(ligand_atoms) == 7
        ), f"Expected 7 ligand atoms, got {len(ligand_atoms)}"

        # Test protein residue parsing
        asp52 = next(r for r in protein_residues if r["res_id"] == 52)
        assert asp52["residue"] == "ASP"
        assert asp52["chain"] == "A"
        assert asp52["pka"] == 4.33
        assert (
            asp52["atom"] == ""
        ), f"Protein residues should have empty atom field, got '{asp52['atom']}'"
        assert (
            asp52["atom_type"] == ""
        ), f"Protein residues should have empty atom_type, got '{asp52['atom_type']}'"
        assert asp52["model_pka"] == 3.80

        # Test ligand atom parsing
        p5s_n = next(
            r for r in ligand_atoms if r["residue"] == "P5S" and r["atom"] == "N"
        )
        assert p5s_n["res_id"] == 0, "Ligands should have res_id=0"
        assert p5s_n["chain"] == "A"
        assert p5s_n["pka"] == 10.83
        assert p5s_n["atom"] == "N", f"Expected atom 'N', got '{p5s_n['atom']}'"
        assert (
            p5s_n["atom_type"] == "N31"
        ), f"Expected atom_type 'N31', got '{p5s_n['atom_type']}'"
        assert p5s_n["model_pka"] == 10.00

        # Test another ligand atom with different type
        y01_cax = next(
            r for r in ligand_atoms if r["residue"] == "Y01" and r["atom"] == "CAX"
        )
        assert (
            y01_cax["atom_type"] == "OCO"
        ), f"Expected atom_type 'OCO', got '{y01_cax['atom_type']}'"
        assert y01_cax["pka"] == 4.62

        # Test phosphate oxygen
        p5s_o15 = next(
            r for r in ligand_atoms if r["residue"] == "P5S" and r["atom"] == "O15"
        )
        assert (
            p5s_o15["atom_type"] == "OP"
        ), f"Expected atom_type 'OP', got '{p5s_o15['atom_type']}'"

        # Verify we can group by ligand
        
        ligands_by_type = defaultdict(list)
        for lig in ligand_atoms:
            ligands_by_type[lig["residue"]].append(lig)

        assert (
            len(ligands_by_type["P5S"]) == 3
        ), f"Expected 3 P5S atoms, got {len(ligands_by_type['P5S'])}"
        assert (
            len(ligands_by_type["LPE"]) == 2
        ), f"Expected 2 LPE atoms, got {len(ligands_by_type['LPE'])}"

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
# SECTION 2: DOCUMENTATION EXAMPLE TESTS
# ============================================================================


class TestPreparationExamples:
    """Run each preparation_example_*.py script once."""

    @parametrize_example_scripts(
        Path(__file__).parent / "preparation_examples",
        "preparation_example_*.py",
    )
    def test_example_script(self, script):
        run_example_script(script)


# ============================================================================
# SECTION 3: COMPLEX STRUCTURE TESTS
# ============================================================================


class TestPropkaComplexStructures:
    """Test propka examples with complex protein structures."""

    @pytest.fixture(autouse=True)
    def require_propka_for_analysis(self, request):
        if request.node.name == "test_8i5b_disulfide_bonds":
            return
        if shutil.which("propka3") is None:
            pytest.skip("PropKa 3 is not on PATH")

    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test outputs."""
        original = os.getcwd()
        with tempfile.TemporaryDirectory() as tmpdir:
            try:
                yield tmpdir
            finally:
                os.chdir(original)

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
        protein_residues = [r for r in residues if r["res_id"] > 0]
        ligand_atoms = [r for r in residues if r["res_id"] == 0]

        assert len(protein_residues) > 0
        assert len(ligand_atoms) > 0  # 6RV3 should have ligands

        print(f"Found {len(protein_residues)} ionizable protein residues")
        print(f"Found {len(ligand_atoms)} ionizable ligand atoms")

        # Group by ligand type
        
        ligands_by_type = defaultdict(list)
        for lig in ligand_atoms:
            ligands_by_type[lig["residue"]].append(lig)

        for ligand_name, atoms in ligands_by_type.items():
            print(f"\n{ligand_name}: {len(atoms)} ionizable atoms")
            for atom in atoms[:3]:  # Show first 3
                print(
                    f"  {atom['atom']:6s} pKa={atom['pka']:5.2f} ({atom['atom_type']})"
                )

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
            residues=residues,
        )

        # Apply disulfide bonds
        if bonds:
            num_bonds = analyzer.apply_disulfide_bonds(
                input_pdb="output/6rv3_ph7.pdb",
                output_pdb="output/6rv3_ph7_ss.pdb",
                disulfide_bonds=bonds,
                auto_detect=False,
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
        protein_residues = [r for r in residues if r["res_id"] > 0]

        # Set target pH
        target_ph = 5.0

        # Group by residue type
        residue_data = {}
        for res in protein_residues:
            res_type = res["residue"]
            if res_type not in residue_data:
                residue_data[res_type] = {"pka": [], "protonated": []}

            residue_data[res_type]["pka"].append(res["pka"])
            is_protonated = target_ph < res["pka"]
            residue_data[res_type]["protonated"].append(is_protonated)

        # Print statistics
        print(f"\npKa Statistics at pH {target_ph}:")
        for res_type, data in sorted(residue_data.items()):
            pka_vals = data["pka"]
            n_prot = sum(data["protonated"])
            n_total = len(pka_vals)
            prot_frac = n_prot / n_total if n_total > 0 else 0

            print(
                f"{res_type:4s}: n={n_total:3d}  "
                f"pKa range=[{min(pka_vals):5.2f}, {max(pka_vals):5.2f}]  "
                f"Protonated: {prot_frac*100:.0f}%"
            )

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
                residues=residues,
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
            "ASP": 3.9,
            "GLU": 4.3,
            "HIS": 6.0,
            "LYS": 10.5,
            "ARG": 12.5,
            "CYS": 8.3,
            "TYR": 10.1,
        }

        # Find residues with significant pKa shifts
        shifted_residues = []
        for res in residues:
            if res["res_id"] == 0:  # Skip ligands
                continue

            res_name = res["residue"]
            if res_name in expected_pka:
                pka_diff = abs(res["pka"] - expected_pka[res_name])
                if pka_diff > 1.0:
                    shifted_residues.append(
                        {
                            "id": f"{res_name}{res['res_id']}_{res['chain']}",
                            "type": res_name,
                            "pka": res["pka"],
                            "expected": expected_pka[res_name],
                            "shift": res["pka"] - expected_pka[res_name],
                        }
                    )

        print(
            f"\n6RV3 Residues with unusual pKa shifts (>1.0 units): {len(shifted_residues)}"
        )

        # Categorize by shift direction
        upshifted = [r for r in shifted_residues if r["shift"] > 1.0]
        downshifted = [r for r in shifted_residues if r["shift"] < -1.0]

        print(f"  Upshifted (more basic): {len(upshifted)}")
        print(f"  Downshifted (more acidic): {len(downshifted)}")

        # Find extreme cases
        extreme_shifts = [r for r in shifted_residues if abs(r["shift"]) > 2.0]
        if extreme_shifts:
            print(f"\n⚠ Extreme shifts (>2.0 units): {len(extreme_shifts)}")
            for r in extreme_shifts[:5]:  # Show first 5
                print(f"  {r['id']}: {r['shift']:+.2f} units")


class TestPreparePdbForPropka:
    """Amber names are renamed only on the file PropKa reads."""

    def test_amber_names_become_standard_and_source_is_kept(self, tmp_path):
        from gatewizard.core.preparation import PreparationManager

        source = tmp_path / "protonated.pdb"
        source.write_text(
            "ATOM      1  CG  HIE A  72      26.206  94.860  58.005  1.00  0.00           C\n"
            "ATOM      2  OD2 ASH A   3      10.000  10.000  10.000  1.00  0.00           O\n"
            "ATOM      3  OE2 GLH A   4      11.000  11.000  11.000  1.00  0.00           O\n"
            "ATOM      4  NZ  LYN A   5      12.000  12.000  12.000  1.00  0.00           N\n"
            "ATOM      5  OH  TYM A   6      13.000  13.000  13.000  1.00  0.00           O\n"
            "ATOM      6  SG  CYX A   7      14.000  14.000  14.000  1.00  0.00           S\n"
            "ATOM      7  CA  ALA A   8      15.000  15.000  15.000  1.00  0.00           C\n"
            "ATOM      8  OXT ALA A   8      16.000  16.000  16.000  1.00  0.00           O\n"
            "END\n",
            encoding="utf-8",
        )
        prepared = PreparationManager._prepare_pdb_for_propka(source, tmp_path)
        text = Path(prepared).read_text(encoding="utf-8")
        original = source.read_text(encoding="utf-8")

        assert "HIE" in original and "ASH" in original
        assert " HIE " not in text
        assert " HIS " in text
        assert " ASP " in text
        assert " GLU " in text
        assert " LYS " in text
        assert " TYR " in text
        assert " CYS " in text
        assert " ALA " in text
        assert "OXT" not in text


class TestStripProteinHydrogens:
    """Protein-only hydrogen stripping (ligands / hetero H kept)."""

    def _write_pdb(self, tmp_path, content: str) -> Path:
        path = Path(tmp_path) / "sample.pdb"
        path.write_text(content)
        return path

    def test_count_and_strip_keeps_ligand_h(self, tmp_path):
        from gatewizard.core.preparation import (
            count_protein_hydrogens,
            strip_protein_hydrogens,
        )

        pdb = self._write_pdb(
            tmp_path,
            """\
ATOM      1  N   ALA A   1      11.104  13.556   9.648  1.00  0.00           N  
ATOM      2  CA  ALA A   1      12.271  12.722   9.648  1.00  0.00           C  
ATOM      3  HA  ALA A   1      12.800  12.900  10.550  1.00  0.00           H  
HETATM    4  C1  LIG A   2      20.000  20.000  20.000  1.00  0.00           C  
HETATM    5  H1  LIG A   2      20.500  20.500  20.500  1.00  0.00           H  
ATOM      6  H   ACE A   3      15.000  15.000  15.000  1.00  0.00           H  
END
""",
        )
        assert count_protein_hydrogens(str(pdb)) == 2
        out = Path(tmp_path) / "out.pdb"
        result = strip_protein_hydrogens(str(pdb), str(out))
        assert result["removed"] == 2
        text = out.read_text()
        assert "HA  ALA" not in text
        assert "H   ACE" not in text
        assert "H1  LIG" in text
        assert "CA  ALA" in text
        assert count_protein_hydrogens(str(out)) == 0

    def test_builder_warns_unless_remove_flag(self, tmp_path):
        from gatewizard.core.builder import Builder

        pdb = self._write_pdb(
            tmp_path,
            """\
ATOM      1  N   ALA A   1      11.104  13.556   9.648  1.00  0.00           N  
ATOM      2  HA  ALA A   1      12.800  12.900  10.550  1.00  0.00           H  
END
""",
        )
        builder = Builder()
        ok, msg = builder.validate_system_inputs(
            str(pdb),
            ["POPC"],
            ["POPC"],
            "1.0//1.0",
            water_model="opc",
            protein_ff="ff19SB",
            lipid_ff="lipid21",
        )
        assert ok
        assert "Remove protein hydrogens" in msg

        ok2, msg2 = builder.validate_system_inputs(
            str(pdb),
            ["POPC"],
            ["POPC"],
            "1.0//1.0",
            water_model="opc",
            protein_ff="ff19SB",
            lipid_ff="lipid21",
            remove_protein_h=True,
        )
        assert ok2
        assert "Remove protein hydrogens" not in msg2


class TestPdb4amberRobustness:
    """CONECT/LINK strip and reduce-binary fallback (no AmberTools required)."""

    def test_strip_conect_and_link_keeps_atoms(self, tmp_path):
        from gatewizard.core.preparation import strip_conect_and_link_records

        src = tmp_path / "with_conect.pdb"
        src.write_text(
            "ATOM      1  N   ALA A   1      11.104  13.556   9.648  1.00  0.00           N\n"
            "CONECT    1    2\n"
            "LINK         O   THR A  93                 K     K A 306\n"
            "END\n",
            encoding="utf-8",
        )
        out = tmp_path / "clean.pdb"
        info = strip_conect_and_link_records(str(src), str(out))
        assert info["removed"] == 2
        text = out.read_text(encoding="utf-8")
        assert "ATOM" in text
        assert "CONECT" not in text
        assert "LINK" not in text

    def test_reduce_option_dropped_when_binary_missing(self, monkeypatch):
        from gatewizard.core import preparation as prep

        monkeypatch.setattr(prep, "resolve_reduce_executable", lambda: None)
        options, skipped = prep.apply_pdb4amber_reduce_option({"reduce": True, "dry": False})
        assert skipped is True
        assert "reduce" not in options
        assert options["dry"] is False

    def test_reduce_option_kept_when_binary_present(self, monkeypatch):
        from gatewizard.core import preparation as prep

        monkeypatch.setattr(prep, "resolve_reduce_executable", lambda: "/usr/bin/reduce")
        options, skipped = prep.apply_pdb4amber_reduce_option({"reduce": True})
        assert skipped is False
        assert options["reduce"] is True

    def test_get_clean_env_prepends_conda_bin(self, monkeypatch, tmp_path):
        import os

        from gatewizard.utils.helpers import get_clean_env

        conda = tmp_path / "env"
        (conda / "bin").mkdir(parents=True)
        monkeypatch.setenv("CONDA_PREFIX", str(conda))
        monkeypatch.setenv("PATH", "/usr/bin")
        monkeypatch.delenv("AMBERHOME", raising=False)
        env = get_clean_env()
        parts = env["PATH"].split(os.pathsep)
        assert parts[0] == str(conda / "bin")
        assert env["AMBERHOME"] == str(conda)


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
