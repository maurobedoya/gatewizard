"""Force-field and charge-model pairing for ligand parametrization."""

import math

from rdkit import Chem

from gatewizard.tools.ligand_parametrization import (
    _bond_orders_for_depiction,
    _load_gaff_mol2,
    build_tleap_ligand_lines,
    forcefield_charge_warning,
    tleap_quote_path,
)


def test_recommended_pairs_have_no_warning():
    assert forcefield_charge_warning("gaff", "bcc") is None
    assert forcefield_charge_warning("gaff2", "abcg2") is None


def test_non_recommended_pairs_warn_and_name_the_amber_pairs():
    for atom_type, charge_method in (("gaff2", "bcc"), ("gaff", "abcg2")):
        warning = forcefield_charge_warning(atom_type, charge_method)
        assert warning is not None
        assert "not a recommended" in warning
        assert "GAFF with AM1-BCC" in warning
        assert "GAFF2 with ABCG2" in warning


def test_mixed_ligand_force_fields_source_both_leaprc():
    lines = build_tleap_ligand_lines(
        {
            "AAA": {
                "frcmod": "AAA.frcmod",
                "lib": "AAA.lib",
                "atom_type": "gaff2",
            },
            "BBB": {
                "frcmod": "BBB.frcmod",
                "lib": "BBB.lib",
                "atom_type": "gaff",
            },
        }
    )
    assert "source leaprc.gaff2" in lines
    assert "source leaprc.gaff" in lines
    assert lines.index("source leaprc.gaff2") < lines.index('loadamberparams "AAA.frcmod"')
    assert 'loadamberparams "BBB.frcmod"' in lines
    assert 'loadoff "BBB.lib"' in lines


def test_tleap_paths_with_spaces_are_quoted():
    spaced = (
        "/mnt/c/Users/kcoru/OneDrive - Universidad Catolica del Maule/"
        "project/ligand_params/NFX/NFX.frcmod"
    )
    lines = build_tleap_ligand_lines(
        {
            "NFX": {
                "frcmod": spaced,
                "lib": spaced.replace(".frcmod", ".lib"),
                "atom_type": "gaff2",
            }
        }
    )
    assert f'loadamberparams "{spaced}"' in lines
    assert f'loadoff "{spaced.replace(".frcmod", ".lib")}"' in lines
    assert "loadamberparams /mnt/c" not in lines
    assert tleap_quote_path("AAA/AAA.frcmod") == '"AAA/AAA.frcmod"'
    assert tleap_quote_path('weird"name.frcmod') == '"weird\\"name.frcmod"'


def _hexagon(radius, origin):
    ox, oy, oz = origin
    return [
        (
            ox + radius * math.cos(i * math.pi / 3),
            oy + radius * math.sin(i * math.pi / 3),
            oz,
        )
        for i in range(6)
    ]


def _mol2_with_bond_type_3(path):
    """Two rings and a real triple bond, all stored as mol2 bond type 3."""
    aromatic = _hexagon(1.40, (0.0, 0.0, 0.0))
    saturated = _hexagon(1.53, (10.0, 0.0, 0.0))
    triple = [(20.0, 0.0, 0.0), (21.20, 0.0, 0.0)]
    carbonyl = [(30.0, 0.0, 0.0), (31.22, 0.0, 0.0)]
    coords = aromatic + saturated + triple + carbonyl
    lines = [
        "@<TRIPOS>MOLECULE",
        "TEST",
        f"{len(coords)} 14",
        "SMALL",
        "NO_CHARGES",
        "",
        "@<TRIPOS>ATOM",
    ]
    for index, (x, y, z) in enumerate(coords, start=1):
        lines.append(f"{index:6d} C{index:<4} {x:10.4f} {y:10.4f} {z:10.4f} ca")
    lines.append("@<TRIPOS>BOND")
    bond_id = 1
    for start in (1, 7):
        for offset in range(6):
            a = start + offset
            b = start + (offset + 1) % 6
            # Alternate 3 and 1 so each carbon has valence 4 and the file loads.
            # The length, not this number, decides the drawn order.
            order = 3 if offset % 2 == 0 else 1
            lines.append(f"{bond_id:6d} {a:5d} {b:5d} {order}")
            bond_id += 1
    lines.append(f"{bond_id:6d} {13:5d} {14:5d} 3")
    lines.append(f"{bond_id + 1:6d} {15:5d} {16:5d} 2")
    path.write_text("\n".join(lines) + "\n")


def test_mol2_bond_type_3_is_drawn_from_the_bond_length(tmp_path):
    mol2 = tmp_path / "rings.mol2"
    _mol2_with_bond_type_3(mol2)
    mol = _bond_orders_for_depiction(_load_gaff_mol2(mol2))

    def kind(a, b):
        return mol.GetBondBetweenAtoms(a, b).GetBondType()

    aromatic = [kind(i, (i + 1) % 6) for i in range(6)]
    assert aromatic.count(Chem.BondType.DOUBLE) == 3
    assert aromatic.count(Chem.BondType.SINGLE) == 3
    assert all(
        kind(6 + i, 6 + (i + 1) % 6) == Chem.BondType.SINGLE for i in range(6)
    )
    assert kind(12, 13) == Chem.BondType.TRIPLE
    assert kind(14, 15) == Chem.BondType.DOUBLE
