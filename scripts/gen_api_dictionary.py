#!/usr/bin/env python3
"""Generate docs/api/dictionary.md from public classes and functions.

Uses the AST of known source files so docs CI does not need the scientific
stack. Re-run after adding a public method:

    python scripts/gen_api_dictionary.py
"""

from __future__ import annotations

import ast
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/api/dictionary.md"

# (section title, markdown page, source path, classes, module-level functions)
INVENTORY: list[tuple[str, str, str, list[str], list[str]]] = [
    (
        "Structure Manager",
        "structure_manager.md",
        "gatewizard/core/structure_manager.py",
        ["StructureManager"],
        ["parse_pdb", "assign_secondary_structure_map"],
    ),
    (
        "Preparation",
        "preparation.md",
        "gatewizard/core/preparation.py",
        ["PreparationManager"],
        [
            "run_propka",
            "extract_summary_section",
            "parse_summary_section",
            "modify_pdb_based_on_summary",
            "count_protein_hydrogens",
            "strip_protein_hydrogens",
            "complete_missing_heavy_atoms",
        ],
    ),
    (
        "Residue mapping",
        "preparation.md",
        "gatewizard/utils/residue_mapping.py",
        [],
        [
            "parse_residue_mapping_by_chain",
            "parse_residue_mapping_file",
        ],
    ),
    (
        "Protein capping",
        "preparation.md",
        "gatewizard/utils/protein_capping.py",
        ["ProteinCapper"],
        ["cap_protein", "detect_terminal_caps", "is_already_capped"],
    ),
    (
        "Builder",
        "builder.md",
        "gatewizard/core/builder.py",
        ["Builder"],
        [],
    ),
    (
        "Job monitor",
        "builder.md",
        "gatewizard/core/job_monitor.py",
        ["JobMonitor"],
        [],
    ),
    (
        "Force fields",
        "builder.md",
        "gatewizard/tools/force_fields.py",
        ["ForceFieldManager"],
        [],
    ),
    (
        "MemPrO",
        "mempro.md",
        "gatewizard/core/mempro.py",
        ["MemPrO"],
        ["compute_orientation_transform", "apply_orientation_transform"],
    ),
    (
        "Hydration",
        "hydration.md",
        "gatewizard/tools/packmol_hydration.py",
        [],
        [
            "check_packmol_available",
            "detect_hydrogen_status",
            "estimate_cavity_volume",
            "build_hydrate_inp_text",
            "prepare_hydration_job",
            "hydrate_cavity",
            "preview_hydrate_inp",
            "run_custom_packmol",
        ],
    ),
    (
        "Equilibration",
        "equilibration.md",
        "gatewizard/tools/equilibration.py",
        [
            "NAMDEquilibrationManager",
            "OpenMMEquilibrationManager",
            "GROMACSEquilibrationManager",
            "AmberEquilibrationManager",
            "EquilibrationStage",
        ],
        ["pdb_has_protein", "adjust_protocol_for_missing_protein"],
    ),
    (
        "Equilibration extras",
        "equilibration.md",
        "gatewizard/utils/equilibration_cluster_script.py",
        [],
        [
            "cluster_engine_executable",
            "write_cluster_run_script",
            "resolve_cluster_launch_script",
        ],
    ),
    (
        "NAMD / energy analysis",
        "analysis.md",
        "gatewizard/utils/namd_analysis.py",
        ["EnergyAnalyzer"],
        ["parse_namd_log", "run_energetic_analysis"],
    ),
    (
        "Trajectory analysis",
        "analysis.md",
        "gatewizard/utils/trajectory_analysis.py",
        ["TrajectoryAnalyzer"],
        ["prepare_structural_inputs", "run_structural_analysis"],
    ),
    (
        "Bilayer analysis",
        "analysis.md",
        "gatewizard/utils/lipid_bilayer_analysis.py",
        ["BilayerTrajectoryAnalyzer"],
        ["run_bilayer_analysis"],
    ),
    (
        "OpenMM analysis",
        "analysis.md",
        "gatewizard/utils/openmm_analysis.py",
        ["OpenMMLogAnalyzer"],
        ["run_openmm_energetic_analysis"],
    ),
    (
        "GROMACS analysis",
        "analysis.md",
        "gatewizard/utils/gromacs_analysis.py",
        ["GROMACSLogEnergyAnalyzer"],
        [],
    ),
    (
        "Amber analysis",
        "analysis.md",
        "gatewizard/utils/amber_analysis.py",
        ["AmberLogEnergyAnalyzer"],
        ["run_amber_energetic_analysis"],
    ),
    (
        "FATSLiM / GridMAT",
        "analysis.md",
        "gatewizard/utils/fatslim_apl.py",
        [],
        ["run_fatslim_apl"],
    ),
    (
        "Cluster helpers",
        "equilibration.md",
        "gatewizard/utils/cluster/ssh.py",
        ["SSHSession"],
        ["connect_ssh", "run_remote", "rsync_to_remote", "rsync_from_remote"],
    ),
]


def _public_methods(cls: ast.ClassDef) -> list[str]:
    names: list[str] = []
    for node in cls.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            if node.name.startswith("_") and node.name != "__init__":
                continue
            names.append(node.name)
    return names


def _public_functions(tree: ast.Module, wanted: list[str]) -> list[str]:
    found: set[str] = set()
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            if node.name in wanted:
                found.add(node.name)
    return [name for name in wanted if name in found]


def _classes_by_name(tree: ast.Module) -> dict[str, ast.ClassDef]:
    return {
        node.name: node for node in tree.body if isinstance(node, ast.ClassDef)
    }


def _anchor(name: str) -> str:
    return name.lower()


def _has_heading(page_text: str, name: str) -> bool:
    return re.search(rf"^### {re.escape(name)}\s*$", page_text, re.M) is not None


def render() -> str:
    lines = [
        "# Method dictionary",
        "",
        "Scan the public Python API by name. Each entry jumps to the heading on",
        "the hand-written module page. Examples live in `tests/<module>_examples/`",
        "and are included on those pages.",
        "",
        "Regenerate with `python scripts/gen_api_dictionary.py`.",
        "",
        '<div class="api-dictionary" markdown="1">',
        "",
    ]
    for title, page, source, class_names, func_names in INVENTORY:
        path = ROOT / source
        if not path.is_file():
            raise SystemExit(f"missing source: {source}")
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=source)
        classes = _classes_by_name(tree)
        page_text = (ROOT / "docs/api" / page).read_text(encoding="utf-8")
        lines.append(f"## {title}")
        lines.append("")
        rows: list[tuple[str, str, str]] = []
        for cls_name in class_names:
            cls = classes.get(cls_name)
            if cls is None:
                raise SystemExit(f"{source}: class {cls_name} not found")
            for method in _public_methods(cls):
                label = f"{cls_name}.{method}"
                if method == "__init__":
                    target = "constructor"
                elif _has_heading(page_text, method):
                    target = _anchor(method)
                else:
                    target = ""
                href = f"{page}#{target}" if target else page
                rows.append((label.lower(), label, href))
        for func in _public_functions(tree, func_names):
            target = _anchor(func) if _has_heading(page_text, func) else ""
            href = f"{page}#{target}" if target else page
            rows.append((func.lower(), func, href))
        extra_source = None
        if title == "FATSLiM / GridMAT":
            extra_source = ROOT / "gatewizard/utils/gridmat_md_apl.py"
        if extra_source and extra_source.is_file():
            extra_tree = ast.parse(extra_source.read_text(encoding="utf-8"))
            extra_page = (ROOT / "docs/api" / page).read_text(encoding="utf-8")
            for func in _public_functions(extra_tree, ["run_gridmat_md_apl"]):
                target = _anchor(func) if _has_heading(extra_page, func) else ""
                href = f"{page}#{target}" if target else page
                rows.append((func.lower(), func, href))
        rows.sort(key=lambda item: item[0])
        if not rows:
            lines.append("_No public names._")
            lines.append("")
            continue
        for _, label, href in rows:
            lines.append(f"- [`{label}`]({href})")
        lines.append("")
    lines.append("</div>")
    lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def main(argv: list[str] | None = None) -> int:
    argv = argv if argv is not None else sys.argv[1:]
    text = render()
    if argv == ["--check"]:
        current = OUT.read_text(encoding="utf-8") if OUT.is_file() else ""
        if current != text:
            print("docs/api/dictionary.md is stale. Run:")
            print("  python scripts/gen_api_dictionary.py")
            return 1
        print("docs/api/dictionary.md is up to date")
        return 0
    OUT.write_text(text, encoding="utf-8")
    print(f"wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
