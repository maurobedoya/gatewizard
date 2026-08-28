"""Tests for AmberEquilibrationManager and amber_analysis helpers."""

import re
from pathlib import Path

import pytest

from gatewizard.tools.equilibration import AmberEquilibrationManager, EquilibrationStage
from gatewizard.utils import amber_analysis
from gatewizard.utils.optional_deps import list_md_engine_candidates, parse_engine_variant

POPC_DIR = Path(__file__).parent / "equilibration_examples" / "popc_membrane"
PRMTOP = POPC_DIR / "system.prmtop"
INPCRD = POPC_DIR / "system.inpcrd"
BILAYER_PDB = POPC_DIR / "bilayer_protein_protonated_prepared_lipid.pdb"
SYSTEM_PDB = POPC_DIR / "system.pdb"

TEMPLATES_DIR = Path(__file__).parent.parent / "equilibration" / "amber"

MDA_AVAILABLE = True
try:
    import MDAnalysis  # noqa: F401
except ImportError:
    MDA_AVAILABLE = False


def _make_manager(tmp_path: Path) -> AmberEquilibrationManager:
    return AmberEquilibrationManager(tmp_path)


# ---------------------------------------------------------------------------
# Templates / constants
# ---------------------------------------------------------------------------


class TestClassConstants:
    def test_scheme_mapping(self):
        assert set(AmberEquilibrationManager.SCHEME_MAPPING.keys()) == {
            "NVT",
            "NPT",
            "NPAT",
            "NPgT",
        }
        assert AmberEquilibrationManager.SCHEME_MAPPING["NPT"] == "NPT"

    def test_stage_index_to_key(self):
        m = AmberEquilibrationManager.STAGE_INDEX_TO_KEY
        assert m[0] == "step0_minimization"
        assert m[7] == "step7_production"


class TestTemplateExistence:
    EQ_FILES = [
        "step0_minimization.mdin",
        "step1_equilibration.mdin",
        "step2_equilibration.mdin",
        "step3_equilibration.mdin",
        "step4_equilibration.mdin",
        "step5_equilibration.mdin",
        "step6_equilibration.mdin",
    ]
    PROD_ENSEMBLES = ["NVT", "NPT", "NPAT", "NPgT"]

    @pytest.mark.parametrize("filename", EQ_FILES)
    def test_eq_template_exists(self, filename):
        p = TEMPLATES_DIR / "eq" / filename
        assert p.exists(), f"Missing template: {p}"

    @pytest.mark.parametrize("ensemble", PROD_ENSEMBLES)
    def test_production_template_exists(self, ensemble):
        p = TEMPLATES_DIR / "production" / ensemble / "step7_production.mdin"
        assert p.exists(), f"Missing template: {p}"

    def test_placeholders_present(self):
        text = (TEMPLATES_DIR / "eq" / "step1_equilibration.mdin").read_text()
        for token in ("{TEMPERATURE}", "{NSTLIM}", "{DT}", "{NTR}", "{RESTRAINT_BLOCK}"):
            assert token in text

    def test_nvt_early_stages_have_no_barostat(self):
        """True NVT: heating stages have no pressure coupling."""
        for fname in ("step1_equilibration.mdin", "step2_equilibration.mdin"):
            text = (TEMPLATES_DIR / "eq" / fname).read_text()
            assert "ntp=" not in text
            assert "barostat=" not in text
            assert "csurften=" not in text

    def test_nvt_production_is_true_nvt(self):
        text = (TEMPLATES_DIR / "production" / "NVT" / "step7_production.mdin").read_text()
        assert "ntp=" not in text
        assert "barostat=" not in text
        assert "csurften=" not in text

    def test_packing_eq6_is_npgt(self):
        text = (TEMPLATES_DIR / "eq" / "step6_equilibration.mdin").read_text()
        assert "ntp=3" in text
        assert "csurften=3" in text

    def test_npat_production_flags(self):
        text = (TEMPLATES_DIR / "production" / "NPAT" / "step7_production.mdin").read_text()
        assert "ntp=2" in text
        assert "baroscalingdir=3" in text

    def test_npgt_packing_flags(self):
        text = (TEMPLATES_DIR / "eq" / "step3_equilibration.mdin").read_text()
        assert "ntp=3" in text
        assert "{GAMMA_TEN}" in text
        assert "csurften=3" in text

    def test_no_dihedral_nmropt(self):
        """Templates must not enable nmropt/DISANG (comments may mention them)."""
        for fname in self.EQ_FILES:
            text = (TEMPLATES_DIR / "eq" / fname).read_text()
            assert "nmropt=" not in text.lower()
            assert "DISANG=" not in text
            assert "dihe.restraint" not in text
        for ens in self.PROD_ENSEMBLES:
            text = (TEMPLATES_DIR / "production" / ens / "step7_production.mdin").read_text()
            assert "nmropt=" not in text.lower()
            assert "DISANG=" not in text
            assert "dihe.restraint" not in text


class TestGetDefaultStageParams:
    def test_default_length_and_timestep_ladder(self):
        stages = AmberEquilibrationManager.get_default_stage_params()
        assert len(stages) == 7
        assert stages[0].name == "Minimization"
        assert stages[0].minimize_steps == 10000
        # 1 fs through eq4, then 2 fs
        for s in stages[1:5]:
            assert s.timestep == 1.0
        assert stages[5].timestep == 2.0
        assert stages[6].timestep == 2.0

    def test_with_production(self):
        stages = AmberEquilibrationManager.get_default_stage_params(include_production=True)
        assert len(stages) == 8
        assert stages[-1].name == "Production"
        assert stages[-1].timestep == 2.0


# ---------------------------------------------------------------------------
# Setup / mdin generation
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not PRMTOP.is_file(), reason="POPC example files missing")
class TestSetupAmberEquilibration:
    def test_setup_writes_mdins_and_run_script(self, tmp_path):
        work = tmp_path / "work"
        work.mkdir()
        for src in (PRMTOP, INPCRD):
            (work / src.name).write_bytes(src.read_bytes())
        pdb_src = SYSTEM_PDB if SYSTEM_PDB.is_file() else BILAYER_PDB
        if pdb_src.is_file():
            (work / "system.pdb").write_bytes(pdb_src.read_bytes())

        mgr = _make_manager(work)
        stages = AmberEquilibrationManager.get_default_stage_params("NPT")
        result = mgr.setup_amber_equilibration(
            stage_params_list=stages,
            output_name="eq_amber",
            amber_executable="pmemd.cuda",
        )
        out = result["amber_dir"]
        assert (out / "run_equilibration.sh").is_file()
        assert (out / "step0_minimization.mdin").is_file()
        assert (out / "step6_equilibration.mdin").is_file()
        script = (out / "run_equilibration.sh").read_text()
        assert "pmemd.cuda" in script
        assert "RESUME" in script
        assert "_gw_amber_stage_done" in script
        assert "-ref $REF" in script

    def test_inpcrd_gets_box_from_bilayer_cryst1(self, tmp_path):
        """Packmol inpcrd often lacks a box line; Amber needs it from CRYST1."""
        work = tmp_path / "work"
        work.mkdir()
        for src in (PRMTOP, INPCRD):
            (work / src.name).write_bytes(src.read_bytes())
        assert BILAYER_PDB.is_file()
        (work / BILAYER_PDB.name).write_bytes(BILAYER_PDB.read_bytes())
        if SYSTEM_PDB.is_file():
            (work / "system.pdb").write_bytes(SYSTEM_PDB.read_bytes())

        # Confirm fixture inpcrd has no box (coords only)
        assert not AmberEquilibrationManager._inpcrd_has_box(work / "system.inpcrd")

        mgr = _make_manager(work)
        stages = AmberEquilibrationManager.get_default_stage_params("NPT")[:2]
        result = mgr.setup_amber_equilibration(
            stage_params_list=stages,
            output_name="eq_box",
            scheme_type="NPT",
            system_files={
                "prmtop": str(work / "system.prmtop"),
                "inpcrd": str(work / "system.inpcrd"),
                "pdb": str(work / "system.pdb") if (work / "system.pdb").is_file() else None,
                "bilayer_pdb": str(work / BILAYER_PDB.name),
            },
        )
        out_inpcrd = result["amber_dir"] / "system.inpcrd"
        assert AmberEquilibrationManager._inpcrd_has_box(out_inpcrd)
        cryst = AmberEquilibrationManager._read_cryst1_cell(work / BILAYER_PDB.name)
        assert cryst is not None
        last = out_inpcrd.read_text().splitlines()[-1]
        vals = [float(x) for x in last.split()]
        assert abs(vals[0] - cryst[0]) < 1e-3
        assert abs(vals[1] - cryst[1]) < 1e-3
        assert abs(vals[2] - cryst[2]) < 1e-3
        out_prmtop = result["amber_dir"] / "system.prmtop"
        assert AmberEquilibrationManager._prmtop_has_box(out_prmtop)

    def test_prmtop_ifbox_text_patch(self, tmp_path):
        """Text fallback sets IFBOX when ParmEd is unavailable."""
        src = tmp_path / "system.prmtop"
        src.write_bytes(PRMTOP.read_bytes())
        cell = AmberEquilibrationManager._read_cryst1_cell(BILAYER_PDB)
        assert cell is not None
        a, b, c, _alpha, beta, _gamma = cell
        patched = AmberEquilibrationManager._patch_prmtop_ifbox_text(
            src.read_text(), a, b, c, beta=beta
        )
        dst = tmp_path / "patched.prmtop"
        dst.write_text(patched)
        assert AmberEquilibrationManager._prmtop_has_box(dst)
        assert "%FLAG BOX_DIMENSIONS" in patched

    def test_nvt_generated_is_true_nvt(self, tmp_path):
        """Generated Amber NVT inputs must not enable a barostat."""
        work = tmp_path / "work"
        work.mkdir()
        for src in (PRMTOP, INPCRD):
            (work / src.name).write_bytes(src.read_bytes())
        mgr = _make_manager(work)
        stages = AmberEquilibrationManager.get_default_stage_params(
            "NVT", include_production=True
        )
        result = mgr.setup_amber_equilibration(
            stage_params_list=stages, output_name="eq_nvt"
        )
        for stem in ("step1_equilibration", "step7_production"):
            text = (result["amber_dir"] / f"{stem}.mdin").read_text()
            assert "ntp=" not in text
            assert "barostat=" not in text
            assert "Generated by GateWizard v" in text
            assert "Templates version: testing" in text
            assert "{GW_VERSION}" not in text
            assert "{GW_GENERATED_ON}" not in text
        packing = (result["amber_dir"] / "step3_equilibration.mdin").read_text()
        assert "ntp=3" in packing
        assert "NPgT SCHEME" in packing

    def test_setup_npgt_scheme_type(self, tmp_path):
        """NPgT must stay NPgT — not NPGT from .upper()."""
        work = tmp_path / "work"
        work.mkdir()
        for src in (PRMTOP, INPCRD):
            (work / src.name).write_bytes(src.read_bytes())
        mgr = _make_manager(work)
        stages = AmberEquilibrationManager.get_default_stage_params("NPgT")[:4]
        result = mgr.setup_amber_equilibration(
            stage_params_list=stages,
            output_name="eq_npgt",
            scheme_type="NPgT",
        )
        eq1 = (result["amber_dir"] / "step1_equilibration.mdin").read_text()
        assert "NVT SCHEME" in eq1
        assert "NPGT SCHEME" not in eq1
        eq3 = (result["amber_dir"] / "step3_equilibration.mdin").read_text()
        assert "NPgT SCHEME" in eq3
        assert "NPGT SCHEME" not in eq3

    def test_setup_npgt_from_npgt_alias(self, tmp_path):
        """All-caps NPGT from legacy callers maps to NPgT."""
        work = tmp_path / "work"
        work.mkdir()
        for src in (PRMTOP, INPCRD):
            (work / src.name).write_bytes(src.read_bytes())
        mgr = _make_manager(work)
        stages = AmberEquilibrationManager.get_default_stage_params("NPgT")[:2]
        result = mgr.setup_amber_equilibration(
            stage_params_list=stages,
            output_name="eq_npgt2",
            scheme_type="NPGT",
        )
        assert (result["amber_dir"] / "run_equilibration.sh").is_file()

    def test_generate_mdin_substitutes_nstlim_and_dt(self, tmp_path):
        mgr = _make_manager(tmp_path)
        content = mgr.generate_mdin_file(
            stage_name="Equilibration 1",
            stage_params={
                "temperature": 303.15,
                "timestep": 1.0,
                "time_ns": 0.125,
                "constraints": {},
            },
            stage_index=1,
            scheme_type="NPT",
        )
        assert "temp0=303.1500" in content
        assert "dt=0.001" in content
        assert "nstlim=125000" in content
        assert "ntr=0" in content
        assert "ntwx=5000" in content

    def test_generate_mdin_eq6_ntwx_from_dcd_freq(self, tmp_path):
        mgr = _make_manager(tmp_path)
        content = mgr.generate_mdin_file(
            stage_name="Equilibration 6",
            stage_params={
                "temperature": 303.15,
                "timestep": 2.0,
                "time_ns": 47.625,
                "dcd_freq": 50000,
                "ensemble": "NPgT",
                "constraints": {"protein_backbone": 0.1},
            },
            stage_index=6,
            scheme_type="NPT",
        )
        assert "ntwx=50000" in content
        assert "nstlim=23812500" in content

    @pytest.mark.skipif(not MDA_AVAILABLE, reason="MDAnalysis not installed")
    def test_group_restraints_with_and_without_custom(self, tmp_path):
        work = tmp_path / "work"
        work.mkdir()
        for src in (PRMTOP, INPCRD):
            (work / src.name).write_bytes(src.read_bytes())
        pdb_src = SYSTEM_PDB if SYSTEM_PDB.is_file() else BILAYER_PDB
        assert pdb_src.is_file()
        (work / "system.pdb").write_bytes(pdb_src.read_bytes())

        mgr = _make_manager(work)
        stages = [
            EquilibrationStage(
                name="Minimization",
                ensemble="NPT",
                time_ns=0.0,
                timestep=1.0,
                temperature=310.15,
                minimize_steps=100,
                constraints={
                    "protein_backbone": 10.0,
                    "protein_sidechain": 0.0,
                    "lipid_head": 2.5,
                    "lipid_tail": 0.0,
                    "water": 0.0,
                    "ions": 0.0,
                    "other": 0.0,
                },
            ),
            EquilibrationStage(
                name="Equilibration 1",
                ensemble="NPT",
                time_ns=0.01,
                timestep=1.0,
                temperature=310.15,
                constraints={
                    "protein_backbone": 5.0,
                    "protein_sidechain": 0.0,
                    "lipid_head": 0.0,
                    "lipid_tail": 0.0,
                    "water": 0.0,
                    "ions": 0.0,
                    "other": 0.0,
                    "ligand_CUSTOM": 1.0,
                },
            ),
        ]
        selections = AmberEquilibrationManager.get_default_selections(str(work / "system.pdb"))
        selections["ligand_CUSTOM"] = "name CA"  # reuse CA as stand-in custom key
        result = mgr.setup_amber_equilibration(
            stage_params_list=stages,
            selections=selections,
            output_name="eq_rest",
        )
        mini = (result["amber_dir"] / "step0_minimization.mdin").read_text()
        assert "ntr=1" in mini
        assert "ATOM " in mini
        assert not re.search(r"ATOM\s+\d+-\d+", mini)
        atom_lines = [ln for ln in mini.splitlines() if ln.startswith("ATOM ")]
        assert atom_lines
        for ln in atom_lines:
            nums = ln.split()[1:]
            assert nums and len(nums) % 2 == 0, ln
            assert all(n.isdigit() and int(n) > 0 for n in nums), ln
        assert "protein backbone" in mini.lower() or "Protein" in mini or "protein" in mini

        eq1 = (result["amber_dir"] / "step1_equilibration.mdin").read_text()
        assert "ntr=1" in eq1
        assert "ligand CUSTOM" in eq1.lower() or "ligand_CUSTOM" in eq1 or "ligand" in eq1.lower()

    def test_zero_fc_omits_ntr(self, tmp_path):
        mgr = _make_manager(tmp_path)
        content = mgr.generate_mdin_file(
            stage_name="Production",
            stage_params={
                "temperature": 310.15,
                "timestep": 2.0,
                "time_ns": 1.0,
                "constraints": {"protein_backbone": 0.0},
            },
            stage_index=7,
            scheme_type="NPT",
            restraint_block="",
        )
        assert "ntr=0" in content


class TestAmberAtomCardFormat:
    def test_pairs_not_hyphens(self):
        lines = AmberEquilibrationManager._format_atom_card(
            [(1, 1), (3, 6), (10, 10)]
        )
        assert lines == ["ATOM 1 1 3 6 10 10"]
        assert all("-" not in line for line in lines)

    def test_wraps_at_seven_pairs(self):
        ranges = [(i, i) for i in range(1, 9)]
        lines = AmberEquilibrationManager._format_atom_card(ranges)
        assert len(lines) == 2
        assert len(lines[0].split()) == 15  # ATOM + 7 pairs
        assert lines[1] == "ATOM 8 8"

    @pytest.mark.skipif(not MDA_AVAILABLE, reason="MDAnalysis not installed")
    def test_contiguous_backbone_is_start_end_pair(self, tmp_path):
        pdb = tmp_path / "ala.pdb"
        pdb.write_text(
            "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
            "ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00  0.00           C\n"
            "ATOM      3  C   ALA A   1       2.000   1.400   0.000  1.00  0.00           C\n"
            "ATOM      4  O   ALA A   1       1.200   2.400   0.000  1.00  0.00           O\n"
            "ATOM      5  CB  ALA A   1       2.000  -0.800   1.200  1.00  0.00           C\n"
            "END\n"
        )
        mgr = _make_manager(tmp_path)
        block = mgr.build_group_restraint_block(
            system_pdb=pdb,
            constraints={"protein_backbone": 0.1},
        )
        assert "0.1000" in block
        assert "ATOM 1 4" in block
        assert "ATOM 1-4" not in block

    @pytest.mark.skipif(not PRMTOP.is_file(), reason="POPC example files missing")
    @pytest.mark.skipif(not MDA_AVAILABLE, reason="MDAnalysis not installed")
    def test_group_restraints_without_pdb_use_prmtop(self, tmp_path):
        work = tmp_path / "work"
        work.mkdir()
        for src in (PRMTOP, INPCRD):
            (work / src.name).write_bytes(src.read_bytes())

        mgr = _make_manager(work)
        stages = [
            EquilibrationStage(
                name="Equilibration 6",
                ensemble="NPgT",
                time_ns=0.01,
                timestep=2.0,
                temperature=310.15,
                constraints={"protein_backbone": 0.1},
            )
        ]
        result = mgr.setup_amber_equilibration(
            stage_params_list=stages,
            output_name="eq_nopdb",
        )
        mdin = (result["amber_dir"] / "step1_equilibration.mdin").read_text()
        assert "ntr=1" in mdin
        assert "ATOM " in mdin
        assert "0.1000" in mdin
        assert not re.search(r"ATOM\s+\d+-\d+", mdin)


class TestRunScriptResources:
    def test_cuda_visible_devices(self, tmp_path):
        mgr = _make_manager(tmp_path)
        path = mgr.generate_run_script(
            amber_dir=tmp_path,
            prmtop_name="system.prmtop",
            inpcrd_name="system.inpcrd",
            stage_stems=["step0_minimization", "step1_equilibration"],
            amber_executable="pmemd.cuda",
            use_gpu=True,
            gpu_id=1,
            num_gpus=2,
        )
        text = path.read_text()
        assert "CUDA_VISIBLE_DEVICES" in text
        assert "1,2" in text
        assert 'MINI_AMBER="pmemd"' in text
        assert 'AMBER="pmemd.cuda"' in text
        assert "$AMBER" in text
        assert "exit code ${ec}" in text

    def test_cpu_eq_gpu_prod_uses_per_stage_executable(self, tmp_path):
        mgr = _make_manager(tmp_path)
        stage_resources = [
            {
                "stage_kind": "minimization",
                "cpu_cores": 6,
                "use_gpu": False,
                "num_gpus": 0,
            },
            {
                "stage_kind": "equilibration",
                "cpu_cores": 6,
                "use_gpu": False,
                "num_gpus": 0,
            },
            {
                "stage_kind": "production",
                "cpu_cores": 1,
                "use_gpu": True,
                "num_gpus": 1,
                "gpu_id": 0,
            },
        ]
        path = mgr.generate_run_script(
            amber_dir=tmp_path,
            prmtop_name="system.prmtop",
            inpcrd_name="system.inpcrd",
            stage_stems=[
                "step0_minimization",
                "step1_equilibration",
                "step7_production",
            ],
            amber_executable="pmemd.cuda",
            use_gpu=True,
            stage_resources=stage_resources,
        )
        text = path.read_text()
        assert "# --- step0_minimization ---" in text
        assert "# --- step1_equilibration ---" in text
        assert "# --- step7_production ---" in text
        mini_block = text.split("# --- step1_equilibration ---")[0]
        eq_block = text.split("# --- step1_equilibration ---")[1].split(
            "# --- step7_production ---"
        )[0]
        prod_block = text.split("# --- step7_production ---")[1]
        assert "$MINI_AMBER" in mini_block
        assert "$MINI_AMBER" in eq_block
        assert "$AMBER" in prod_block
        assert "CUDA_VISIBLE_DEVICES" in prod_block
        assert "CPU stages use $MINI_AMBER" in text

    def test_cpu_fallback_keeps_cpu_minimization(self, tmp_path):
        mgr = _make_manager(tmp_path)
        path = mgr.generate_run_script(
            amber_dir=tmp_path,
            prmtop_name="system.prmtop",
            inpcrd_name="system.inpcrd",
            stage_stems=["step0_minimization"],
            amber_executable="pmemd.cuda",
            use_gpu=False,
        )
        text = path.read_text()
        assert 'MINI_AMBER="pmemd"' in text
        assert 'AMBER="pmemd.cuda"' in text


# ---------------------------------------------------------------------------
# Discovery / analysis
# ---------------------------------------------------------------------------


class TestAmberDiscovery:
    def test_variant_labels(self):
        assert parse_engine_variant("", "amber", "/x/pmemd.cuda") == "CUDA"
        assert parse_engine_variant("", "amber", "pmemd") == "CPU"
        assert parse_engine_variant("", "amber", "pmemd.MPI") == "MPI"
        assert parse_engine_variant("", "amber", "pmemd.cuda.MPI") == "CUDA+MPI"
        assert parse_engine_variant("", "amber", "sander") == "CPU"

    def test_list_candidates_returns_list(self):
        results = list_md_engine_candidates("amber")
        assert isinstance(results, list)
        for item in results:
            assert "executable" in item
            assert "variant" in item


SAMPLE_MDOUT = """
  NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =   300.00  PRESS =     1.0
 Etot   =     -1000.0000  EKtot   =       500.0000  EPtot      =     -1500.0000
 BOND   =        10.0000  ANGLE   =        20.0000  DIHED      =         5.0000
 ------------------------------------------------------------------------------
  NSTEP =     1000   TIME(PS) =       1.000  TEMP(K) =   301.00  PRESS =     1.1
 Etot   =      -990.0000  EKtot   =       505.0000  EPtot      =     -1495.0000
 BOND   =        11.0000  ANGLE   =        21.0000  DIHED      =         6.0000
 ------------------------------------------------------------------------------

|  Final Results
|  TIMINGS
|     ns/day =  12.5
"""

SAMPLE_MIN_MDOUT = """
   NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER
      1       5.7760E+09     4.6877E+08     1.1352E+11     H8R     11288

 BOND    =     7448.2046  ANGLE   =    25219.0067  DIHED      =    16870.4361
 VDWAALS = *************  EEL     =   -28209.9767  HBOND      =        0.0000
 1-4 VDW =     6650.0882  1-4 EEL =   -17611.4169  RESTRAINT  =        0.0000


   NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER
    100      -1.0420E+05     4.0344E+00     7.1559E+02     C12     33338

 BOND    =    11472.8704  ANGLE   =    10601.2054  DIHED      =    16832.7527
 VDWAALS =    43102.3105  EEL     =  -174989.3218  HBOND      =        0.0000
 1-4 VDW =     6682.6720  1-4 EEL =   -18037.0376  RESTRAINT  =      137.0205
 EAMBER  =  -104334.5483
"""


class TestAmberAnalysis:
    def test_parse_completed(self, tmp_path):
        mdout = tmp_path / "step1_equilibration.mdout"
        mdin = tmp_path / "step1_equilibration.mdin"
        mdin.write_text(" &cntrl\n nstlim=1000,\n dt=0.001,\n /\n")
        mdout.write_text(SAMPLE_MDOUT)
        info = amber_analysis.parse_amber_mdout(mdout, mdin_file=mdin)
        assert info.steps_completed == 1000
        assert info.total_steps == 1000
        assert info.completed is True
        assert info.ns_per_day == pytest.approx(12.5)

    def test_parse_interrupted_not_completed(self, tmp_path):
        mdout = tmp_path / "step1_equilibration.mdout"
        mdin = tmp_path / "step1_equilibration.mdin"
        mdin.write_text(" &cntrl\n nstlim=5000,\n dt=0.001,\n /\n")
        # Has Final Results but only reached 1000 / 5000
        mdout.write_text(SAMPLE_MDOUT)
        info = amber_analysis.parse_amber_mdout(mdout, mdin_file=mdin)
        assert info.steps_completed == 1000
        assert info.interrupted is True
        assert info.completed is False

    def test_energetic_properties(self, tmp_path):
        mdout = tmp_path / "run.mdout"
        mdout.write_text(SAMPLE_MDOUT)
        props = amber_analysis.list_amber_energy_properties([mdout])
        assert "EPtot" in props
        assert "TEMP" in props
        result = amber_analysis.run_amber_energetic_analysis(
            [mdout], properties=["EPtot", "TEMP"]
        )
        assert result["series"]
        assert len(result["x"]) == 2

    def test_file_times_override_time_axis(self, tmp_path):
        """UI per-file durations (ns) must rescale X like NAMD/GROMACS."""
        eq1 = tmp_path / "step1_equilibration.mdout"
        eq2 = tmp_path / "step2_equilibration.mdout"
        # mdout TIME spans only 1 ps, but UI assigns 1 ns each
        eq1.write_text(SAMPLE_MDOUT)
        eq2.write_text(SAMPLE_MDOUT)
        result = amber_analysis.run_amber_energetic_analysis(
            [eq1, eq2],
            properties=["TEMP"],
            file_times={
                "step1_equilibration.mdout": 1.0,
                "step2_equilibration.mdout": 1.0,
            },
        )
        x = result["x"]
        assert len(x) == 4
        assert x[0] == pytest.approx(0.0)
        assert x[1] == pytest.approx(1.0)  # end of file 1
        assert x[2] == pytest.approx(1.0)  # start of file 2
        assert x[3] == pytest.approx(2.0)  # end of file 2

    def test_minimization_mdout_respects_file_times(self, tmp_path):
        """Min ENERGY tables have no TIME(PS); still consume assigned ns."""
        mini = tmp_path / "step0_minimization.mdout"
        eq1 = tmp_path / "step1_equilibration.mdout"
        mini.write_text(SAMPLE_MIN_MDOUT)
        eq1.write_text(SAMPLE_MDOUT)
        props = amber_analysis.list_amber_energy_properties([mini])
        assert "ENERGY" in props or "Etot" in props
        assert "BOND" in props
        result = amber_analysis.run_amber_energetic_analysis(
            [mini, eq1],
            properties=["Etot", "BOND"],
            file_times={
                "step0_minimization.mdout": 1.0,
                "step1_equilibration.mdout": 1.0,
            },
        )
        x = result["x"]
        assert len(x) == 4  # 2 min + 2 eq frames
        assert x[0] == pytest.approx(0.0)
        assert x[1] == pytest.approx(1.0)
        assert x[2] == pytest.approx(1.0)
        assert x[3] == pytest.approx(2.0)

    def test_minimization_reads_mdinfo_for_live_steps(self, tmp_path):
        """mdinfo is rewritten each ntpr; mdout may lag — use max of both."""
        mdout = tmp_path / "step0_minimization.mdout"
        mdinfo = tmp_path / "step0_minimization.mdinfo"
        mdin = tmp_path / "step0_minimization.mdin"
        mdin.write_text(" &cntrl\n imin=1,\n maxcyc=5000,\n ntpr=50,\n /\n")
        mdout.write_text(
            "   NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER\n"
            "    100      -1.8000E+05     5.0000E+00     1.0000E+03     C19     39162\n"
        )
        mdinfo.write_text(
            "   NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER\n"
            "    200      -1.8756E+05     4.0209E+00     1.1031E+03     C19     39162\n"
            "\n"
            " BOND    =    12438.4496  ANGLE   =     6302.9898  DIHED      =    16270.1949\n"
        )
        info = amber_analysis.parse_amber_mdout(
            mdout, is_minimization=True, mdin_file=mdin, mdinfo_file=mdinfo
        )
        assert info.steps_completed == 200
        assert info.total_steps == 5000
        assert info.completed is False

        prog = amber_analysis.get_equilibration_progress(tmp_path)
        assert prog["minimization"].status == "running"
        assert prog["minimization"].log_file == mdinfo
        assert prog["minimization"].timing.steps_completed == 200

    def test_running_production_uses_mdinfo_cumulative_timings(self, tmp_path):
        """Running production must not estimate ns/day from mdinfo mtime (~1 min)."""
        mdout = tmp_path / "step7_production.mdout"
        mdinfo = tmp_path / "step7_production.mdinfo"
        mdin = tmp_path / "step7_production.mdin"
        mdin.write_text(" &cntrl\n nstlim=100000000,\n dt=0.002,\n /\n")
        # Partial mdout without TIMINGS — typical while production is still running.
        mdout.write_text(
            "  NSTEP =  45000000   TIME(PS) =   90000.000  TEMP(K) =   300.00\n"
            " Etot   = -1000.0000  EPtot      =     -1500.0000\n"
        )
        mdinfo.write_text(
            "| Current Timing Info\n"
            "| Total steps : 100000000 | Completed : 46500000 | Remaining : 53500000\n"
            "|\n"
            "| Average timings for last    5000 steps:\n"
            "|     Elapsed(s) =      64.4 Per Step(ms) =      12.9\n"
            "|         ns/day =     130.0   seconds/ns =    663.5\n"
            "|\n"
            "| Average timings for all steps:\n"
            "|     Elapsed(s) =   62300.0 Per Step(ms) =      13.0\n"
            "|         ns/day =     129.0   seconds/ns =    669.0\n"
        )
        info = amber_analysis.parse_amber_mdout(
            mdout, mdin_file=mdin, mdinfo_file=mdinfo
        )
        assert info.steps_completed == 46500000
        assert info.ns_per_day == pytest.approx(129.0)
        assert info.wall_elapsed_seconds == pytest.approx(62300.0)
        assert info.completed is False
        sim_ns = info.steps_completed * info.timestep_fs * 1e-6
        assert sim_ns == pytest.approx(93.0)

    def test_completed_stage_prefers_wall_derived_ns_per_day(self, tmp_path):
        """Long CPU stages: last-window ns/day must not override total wall rate."""
        mdout = tmp_path / "step3_equilibration.mdout"
        mdin = tmp_path / "step3_equilibration.mdin"
        # 0.25 ns = 125000 steps * 2 fs
        mdin.write_text(" &cntrl\n nstlim=125000,\n dt=0.002,\n /\n")
        wall_s = 38 * 3600 + 59 * 60 + 14  # ~39 h
        mdout.write_text(
            "  NSTEP =    125000   TIME(PS) =     250.000  TEMP(K) =   300.00\n"
            " Etot   = -1000.0000  EPtot      =     -1500.0000\n"
            "\n"
            "| Average timings for last    5000 steps:\n"
            "|     Elapsed(s) =     200.0 Per Step(ms) =      40.0\n"
            "|         ns/day =     210.2   seconds/ns =    411.0\n"
            "\n"
            "|  Final Results\n"
            "|  TIMINGS\n"
            f"|     Elapsed(wallclock) = {wall_s:.1f} seconds\n"
            "|     ns/day =     210.2\n"  # bogus leftover / window rate in TIMINGS
        )
        info = amber_analysis.parse_amber_mdout(mdout, mdin_file=mdin)
        assert info.completed is True
        assert info.wall_elapsed_seconds == pytest.approx(wall_s)
        # 0.25 ns / (wall_s/86400) ≈ 0.154 ns/day — not 210.2
        expected = 0.25 / (wall_s / 86400.0)
        assert info.ns_per_day == pytest.approx(expected, rel=1e-3)
        assert info.ns_per_day < 1.0

    def test_mdinfo_last_window_does_not_set_wall(self, tmp_path):
        """Last-window Elapsed(s) is not stage wall time."""
        mdout = tmp_path / "step7_production.mdout"
        mdinfo = tmp_path / "step7_production.mdinfo"
        mdin = tmp_path / "step7_production.mdin"
        mdin.write_text(" &cntrl\n nstlim=1000000,\n dt=0.002,\n /\n")
        mdout.write_text(
            "  NSTEP =   500000   TIME(PS) =    1000.000  TEMP(K) =   300.00\n"
        )
        mdinfo.write_text(
            "| Current Timing Info\n"
            "| Total steps : 1000000 | Completed : 500000 | Remaining : 500000\n"
            "|\n"
            "| Average timings for last    5000 steps:\n"
            "|     Elapsed(s) =      64.4 Per Step(ms) =      12.9\n"
            "|         ns/day =     130.0   seconds/ns =    663.5\n"
        )
        info = amber_analysis.parse_amber_mdout(
            mdout, mdin_file=mdin, mdinfo_file=mdinfo
        )
        assert info.ns_per_day == pytest.approx(130.0)
        # No cumulative Elapsed — wall is estimated from live rate × simulated ns.
        assert info.wall_elapsed_seconds == pytest.approx((1.0 / 130.0) * 86400.0, rel=1e-3)
        assert info.completed is False

    def test_parse_large_mdout_uses_head_and_tail(self, tmp_path):
        """Multi-MB ENERGY dumps must not hide nstlim in the head or TIMINGS in the tail."""
        mdout = tmp_path / "step7_production.mdout"
        mdin = tmp_path / "step7_production.mdin"
        mdin.write_text(" &cntrl\n nstlim=1000,\n dt=0.001,\n /\n")
        header = "  nstlim = 1000\n  dt     = 0.001\n"
        footer = (
            "  NSTEP =     1000   TIME(PS) =       1.000  TEMP(K) =   301.00\n"
            "|  Final Results\n"
            "|  TIMINGS\n"
            "|     Elapsed(wallclock) = 100.0 seconds\n"
            "|     ns/day =  12.5\n"
        )
        junk = b"  NSTEP =      123   TIME(PS) =       0.123  TEMP(K) =   300.00\n"
        with mdout.open("wb") as handle:
            handle.write(header.encode("utf-8"))
            target = 6 * 1024 * 1024
            while handle.tell() < target:
                handle.write(junk)
            handle.write(b"\n")
            handle.write(footer.encode("utf-8"))
        info = amber_analysis.parse_amber_mdout(mdout, mdin_file=mdin)
        assert info.steps_completed == 1000
        assert info.total_steps == 1000
        assert info.completed is True
        assert info.ns_per_day > 0
