from pathlib import Path

from gatewizard.core.preparation import PreparationManager

here = Path(__file__).resolve().parent
out_dir = here / "output"
out_dir.mkdir(exist_ok=True)

sequential = out_dir / "preserve_numbers_sequential.pdb"
sequential.write_text(
    "ATOM      1  CA  ALA A   1       1.000   1.000   1.000  1.00  0.00           C\n"
    "ATOM      2  CA  GLY A   2       2.000   2.000   2.000  1.00  0.00           C\n"
    "HETATM    3  C1  LIG A   3       3.000   3.000   3.000  1.00  0.00           C\n"
    "ATOM      4  C   ACE B   4       4.000   4.000   4.000  1.00  0.00           C\n"
    "ATOM      5  CA  VAL B   5       5.000   5.000   5.000  1.00  0.00           C\n"
    "ATOM      6  N   NME B   6       6.000   6.000   6.000  1.00  0.00           N\n"
    "END\n",
    encoding="utf-8",
)
remum = out_dir / "preserve_numbers_sequential_renum.txt"
remum.write_text(
    "ALA A      1    ALA    10\n"
    "GLY A      2    GLY    20\n"
    "LIG A      3    LIG   900\n"
    "ACE B      4    ACE     1\n"
    "VAL B      5    VAL     2\n"
    "NME B      6    NME     3\n",
    encoding="utf-8",
)
cap_map = out_dir / "preserve_numbers_gatewizard_residue_mapping.txt"
cap_map.write_text(
    "# Format: ORIGINAL_RESNAME CHAIN ORIGINAL_ID FINAL_RESNAME FINAL_ID\n"
    "  - B      -    ACE     1\n"
    "VAL B     21    VAL     2\n"
    "  - B      -    NME     3\n",
    encoding="utf-8",
)

restored = out_dir / "preserve_numbers_restored.pdb"
manager = PreparationManager()
result = manager.restore_original_residue_numbers(
    str(sequential),
    remum_path=str(remum),
    output_pdb=str(restored),
    cap_mapping_path=str(cap_map),
)

print(f"✓ Restored residue numbers: {result['output_file']}")
print(f"  Cap assignments: {result['cap_assignments']}")
assert result["residue_numbers_preserved"]
assert restored.is_file()
text = restored.read_text(encoding="utf-8")
assert "ALA A  10" in text
assert "GLY A  20" in text
assert "LIG A 900" in text
assert "ACE B  20" in text
assert "VAL B  21" in text
assert "NME B  22" in text
