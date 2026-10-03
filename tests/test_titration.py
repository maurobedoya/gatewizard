"""Titration figure filtering matches the GUI panels, without running PropKa."""

from gatewizard.core.titration import (
    assemble_titration,
    disulfide_bridges,
    disulfide_report,
    titration_figure,
)


class _Manager:
    def get_ph_titration_curve(self, ph_range, ph_step):
        self.called = (ph_range, ph_step)
        return {}

    def get_default_protonation_state(self, residue, ph):
        kind = residue["residue"]
        pka = residue["pka"]
        if kind == "CYS":
            return "CYS" if ph < pka else "CYM"
        if kind == "ASP":
            return "ASH" if ph < pka else "ASP"
        return kind


def _payload():
    residues = [
        {"residue": "CYS", "res_id": 77, "chain": "A", "pka": 99.99},
        {"residue": "CYS", "res_id": 449, "chain": "A", "pka": 99.99},
        {"residue": "CYS", "res_id": 12, "chain": "A", "pka": 8.4},
        {"residue": "ASP", "res_id": 10, "chain": "A", "pka": 4.0},
    ]
    return assemble_titration(
        _Manager(),
        residues,
        target_ph=7.0,
        ph_min=6.0,
        ph_max=8.0,
        ph_step=1.0,
    )


def test_exclude_drops_residues_from_panels_a_and_b_only():
    payload = _payload()
    figure = titration_figure(payload, exclude_ids=["CYS77", "CYS449:A"])
    names = {row["name"]: row for row in figure["types"]}
    assert names["CYS"]["pka"] == [8.4]
    assert "CYS77:A" in figure["excluded_ids"]
    assert "CYS449:A" in figure["excluded_ids"]
    curve_ids = [row["id"] for row in figure["curves"]]
    assert "CYS77:A" in curve_ids
    assert "CYS12:A" in curve_ids


def test_curve_ids_and_type_filter_select_panel_c():
    payload = _payload()
    figure = titration_figure(
        payload,
        curve_ids=["ASP10"],
        type_filter=["ASP"],
    )
    assert [row["name"] for row in figure["types"]] == ["ASP"]
    assert [row["id"] for row in figure["curves"]] == ["ASP10:A"]


def test_empty_curve_ids_draws_no_curves():
    payload = _payload()
    figure = titration_figure(payload, curve_ids=[])
    assert figure["curves"] == []
    assert {row["name"] for row in figure["types"]} == {"ASP", "CYS"}


def _sg(serial, resname, chain, resseq, x, y, z):
    return (
        f"ATOM  {serial:5d}  SG  {resname:>3s} {chain}{resseq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           S\n"
    )


def test_disulfide_bridges_use_distance_and_ssbond(tmp_path):
    pdb = tmp_path / "bridges.pdb"
    pdb.write_text(
        "SSBOND   1 CYS A   50    CYS A   51                          1555   1555  2.03\n"
        + _sg(1, "CYS", "A", 77, 0.0, 0.0, 0.0)
        + _sg(2, "CYS", "A", 449, 2.0, 0.0, 0.0)
        + _sg(3, "CYS", "A", 12, 20.0, 0.0, 0.0)
        + _sg(4, "CYX", "B", 3, 0.0, 10.0, 0.0)
        + _sg(5, "CYX", "B", 4, 0.0, 12.0, 0.0)
        + "END\n",
        encoding="utf-8",
    )
    bonds = disulfide_bridges(str(pdb))
    pairs = {(b["res_id_a"], b["res_id_b"], b["chain_a"]) for b in bonds}
    assert (77, 449, "A") in pairs
    assert (3, 4, "B") in pairs
    assert (50, 51, "A") in pairs
    assert all(12 not in (b["res_id_a"], b["res_id_b"]) for b in bonds)


def test_disulfide_report_hides_bridged_cysteines_and_explains_propka():
    payload = _payload()
    bridges = [
        {
            "residue_a": "CYS",
            "chain_a": "A",
            "res_id_a": 77,
            "residue_b": "CYS",
            "chain_b": "A",
            "res_id_b": 449,
            "distance": 2.02,
            "source": "distance",
        }
    ]
    report = disulfide_report(payload["residues"], bridges)
    hidden = {row["id"] for row in report["hidden"]}
    assert hidden == {"CYS77:A", "CYS449:A"}
    assert report["warning"].count("99.99") == 1
    assert "disulfide" in report["warning"]
    assert "CYS12" not in report["warning"]
    figure = titration_figure(
        payload, exclude_ids=[row["id"] for row in report["hidden"]]
    )
    names = {row["name"]: row for row in figure["types"]}
    assert names["CYS"]["pka"] == [8.4]
