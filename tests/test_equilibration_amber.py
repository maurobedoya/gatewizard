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
        assert AmberEquilibrationManager.SCHEME_MAPPING["NPT"] == "02_NPT"

    def test_stage_index_to_key(self):
        m = AmberEquilibrationManager.STAGE_INDEX_TO_KEY
        assert m[0] == "step0_minimization"
        assert m[7] == "step7_production"


class TestTemplateExistence:
    ENSEMBLES = ["01_NVT", "02_NPT", "03_NPAT", "04_NPgT"]
    FILENAMES = [
        "step0_minimization.mdin",
        "step1_equilibration.mdin",
        "step2_equilibration.mdin",
        "step3_equilibration.mdin",
        "step4_equilibration.mdin",
        "step5_equilibration.mdin",
        "step6_equilibration.mdin",
        "step7_production.mdin",
    ]

    @pytest.mark.parametrize("ensemble", ENSEMBLES)
    @pytest.mark.parametrize("filename", FILENAMES)
    def test_template_exists(self, ensemble, filename):
        p = TEMPLATES_DIR / ensemble / filename
        assert p.exists(), f"Missing template: {p}"

    def test_placeholders_present(self):
        text = (TEMPLATES_DIR / "02_NPT" / "step1_equilibration.mdin").read_text()
        for token in ("{TEMPERATURE}", "{NSTLIM}", "{DT}", "{NTR}", "{RESTRAINT_BLOCK}"):
            assert token in text

    def test_nvt_early_stages_have_no_barostat(self):
        """True NVT: heating stages have no pressure coupling."""
        for fname in ("step1_equilibration.mdin", "step2_equilibration.mdin"):
            text = (TEMPLATES_DIR / "01_NVT" / fname).read_text()
            assert "ntp=" not in text
            assert "barostat" not in text
            assert "csurften" not in text

    def test_nvt_all_stages_are_true_nvt(self):
        """GateWizard NVT pack keeps constant volume through production."""
        for fname in self.FILENAMES:
            if fname == "step0_minimization.mdin":
                continue
            text = (TEMPLATES_DIR / "01_NVT" / fname).read_text()
            assert "ntp=" not in text, fname
            assert "barostat=" not in text, fname
            assert "csurften=" not in text, fname

    def test_npat_flags(self):
        text = (TEMPLATES_DIR / "03_NPAT" / "step3_equilibration.mdin").read_text()
        assert "ntp=2" in text
        assert "baroscalingdir=3" in text

    def test_npgt_flags(self):
        text = (TEMPLATES_DIR / "04_NPgT" / "step3_equilibration.mdin").read_text()
        assert "ntp=3" in text
        assert "{GAMMA_TEN}" in text
        assert "csurften=3" in text

    def test_no_dihedral_nmropt(self):
        """Templates must not enable nmropt/DISANG (comments may mention them)."""
        for ens in self.ENSEMBLES:
            for fname in self.FILENAMES:
                text = (TEMPLATES_DIR / ens / fname).read_text()
                assert "nmropt=" not in text.lower()
                assert "DISANG=" not in text
                assert "dihe.restraint" not in text


class TestGetDefaultStageParams:
    def test_default_length_and_timestep_ladder(self):
        stages = AmberEquilibrationManager.get_default_stage_params()
        assert len(stages) == 7
        assert stages[0].name == "Minimization"
        assert stages[0].minimize_steps == 5000
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
        for stem in ("step1_equilibration", "step3_equilibration", "step7_production"):
            text = (result["amber_dir"] / f"{stem}.mdin").read_text()
            assert "ntp=" not in text
            assert "barostat=" not in text
            assert "Generated by GateWizard v" in text
            assert "Templates version: v1" in text
            assert "{GW_VERSION}" not in text
            assert "{GW_GENERATED_ON}" not in text

    def test_setup_npgt_scheme_type(self, tmp_path):
        """NPgT must stay NPgT — not NPGT from .upper()."""
        work = tmp_path / "work"
        work.mkdir()
        for src in (PRMTOP, INPCRD):
            (work / src.name).write_bytes(src.read_bytes())
        mgr = _make_manager(work)
        stages = AmberEquilibrationManager.get_default_stage_params("NPgT")[:2]
        result = mgr.setup_amber_equilibration(
            stage_params_list=stages,
            output_name="eq_npgt",
            scheme_type="NPgT",
        )
        mdin = (result["amber_dir"] / "step1_equilibration.mdin").read_text()
        assert "NPgT SCHEME" in mdin
        assert "NPGT SCHEME" not in mdin

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
        assert "MINI_AMBER=" in text


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
