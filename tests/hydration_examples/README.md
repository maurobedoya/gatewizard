# Packmol Hydration Examples

Example scripts corresponding to **Example 1–6** in `docs/api/hydration.md`.

Each file `hydration_example_XX.py` maps to **Example XX** in the API documentation.

## Examples

| File | Description |
|------|-------------|
| `hydration_example_01.py` | Check PACKMOL availability |
| `hydration_example_02.py` | Detect hydrogen status on 6RV3_AB |
| `hydration_example_03.py` | Estimate cavity volume in a fixed box |
| `hydration_example_04.py` | Build PACKMOL input text (preview) |
| `hydration_example_05.py` | Prepare hydration job files |
| `hydration_example_06.py` | Full hydrate_cavity workflow (requires PACKMOL) |

## Test PDB

All structure examples use [`tests/6RV3_AB.pdb`](../6RV3_AB.pdb) (multi-chain membrane protein).

Fixed box (Å) for examples 3–6 — 20 Å cube around the structure centroid:

- `BOX_MIN = (-13.03, -49.98, 18.67)`
- `BOX_MAX = (6.97, -29.98, 38.67)`

## Running

```bash
conda activate gatewizard
cd /path/to/gatewizard

python tests/hydration_examples/hydration_example_01.py
pytest tests/test_packmol_hydration.py -v
pytest tests/test_packmol_hydration.py::TestHydrationExamples -v
```

Example 6 skips the PACKMOL run when the executable is not installed (CI-safe).
