"""
Builder Example 25: Blocking (synchronous) preparation with wait=True

Demonstrates how to run prepare_system() in blocking mode so that the
Python script does not continue until the job finishes (or errors).
This is especially useful for scripting multiple sequential preparations.

Two approaches are shown:
  A) Inline ``wait=True`` inside prepare_system()
  B) Launch in background + call wait_for_completion() later

NOTE: This example requires AmberTools and packmol-memgen.
      It uses the real test PDB file tests/2MVJ_2ligs.pdb.
"""

from gatewizard.core.builder import Builder
from gatewizard.tools.ligand_parametrization import (
    detect_ligands,
    parametrize_all_ligands,
)

pdb_file = "tests/2MVJ_2ligs.pdb"
working_dir = "./systems"

# ── Step 1: Detect and parametrize ligands (same as example 23) ──────
print("Step 1: Detecting ligands...")
ligands = detect_ligands(pdb_file)
for lig in ligands:
    print(f"  Found: {lig.name} ({lig.num_atoms} atoms, {lig.formula})")

print("\nStep 2: Parametrizing ligands...")
ligand_results = parametrize_all_ligands(
    pdb_file=pdb_file,
    output_dir=f"{working_dir}/ligand_params",
    charges={"AAA": 0, "BBB": 0},
    charge_method="bcc",
    atom_type="gaff2",    # GAFF2 atom types
)
print(f"  Parametrized: {list(ligand_results.keys())}")

# ── Step 3: Configure builder ────────────────────────────────────────
builder = Builder()
builder.set_configuration(
    water_model="tip3p",
    protein_ff="ff14SB",
    lipid_ff="lipid21",
    preoriented=True,
    parametrize=True,
    salt_concentration=0.15,
    dist_wat=17.5,
    notprotonate=True,
    ligand_params=ligand_results,
)

# ── Approach A: inline wait ──────────────────────────────────────────
# With wait=True the call blocks until the job finishes.
# prepare_system returns (success, message, job_dir) where
#   success reflects the *final* outcome, not just "launched OK".
#
# NOTE: The actual preparation is long-running.  Set a short timeout
#       so the example returns quickly during tests; remove or increase
#       the timeout for real production runs.
print("\n── Approach A: prepare_system with wait=True ──")
success_a, message_a, job_dir_a = builder.prepare_system(
    pdb_file=pdb_file,
    working_dir=working_dir,
    upper_lipids=["POPC"],
    lower_lipids=["POPC"],
    lipid_ratios="1//1",
    output_folder_name="membrane_2MVJ_wait",
    wait=True,                # ← block until done or error
    wait_timeout=10,          # short timeout for demo (use 3600+ for real runs)
    wait_poll_interval=2,     # check every 2 s
    wait_verbose=True,        # print elapsed time
)
print(f"Result : {success_a}")
print(f"Message: {message_a}")

# ── Approach B: launch then wait later ───────────────────────────────
print("\n── Approach B: launch + wait_for_completion() ──")
success_b, message_b, job_dir_b = builder.prepare_system(
    pdb_file=pdb_file,
    working_dir=working_dir,
    upper_lipids=["POPC"],
    lower_lipids=["POPC"],
    lipid_ratios="1//1",
    output_folder_name="membrane_2MVJ_bg",
)

if success_b and job_dir_b is not None:
    # Do other work here …
    print("Doing other work while system builds …")

    # Then block until that specific job is done
    completed, wait_msg = builder.wait_for_completion(
        job_dir_b,
        poll_interval=2,
        timeout=10,          # short timeout for demo
        verbose=True,
    )
    print(f"Completed: {completed}")
    print(f"Message  : {wait_msg}")
else:
    print(f"Not launched: {message_b}")

print("\nExample 25 finished.")
