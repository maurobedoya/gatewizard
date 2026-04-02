# MemPrO Test Examples

This directory contains example scripts corresponding to the documentation
examples in `docs/api/mempro.md`.

Each example file (`mempro_example_XX.py`) maps to **Example XX** in the
API documentation.

## Examples

| File | Description |
|------|-------------|
| `mempro_example_01.py` | Create MemPrO instance and check availability |
| `mempro_example_02.py` | Check if MemPrO is installed |
| `mempro_example_03.py` | Run MemPrO with default settings |
| `mempro_example_04.py` | Run with custom parameters |
| `mempro_example_05.py` | Run with dual membrane mode |
| `mempro_example_06.py` | Run with peripheral mode |
| `mempro_example_07.py` | Parse results from existing Orient directory |
| `mempro_example_08.py` | Get best oriented PDB path |
| `mempro_example_09.py` | Build and inspect command |
| `mempro_example_10.py` | Inspect OrientationResult attributes |
| `mempro_example_11.py` | Full orientation workflow |
| `mempro_example_12.py` | Parse results and load into viewer |

## Running Examples

```bash
# Examples 1-2, 7-10 work without mempro installed (unit-testable)
python tests/mempro_examples/mempro_example_01.py

# Examples 3-6, 11 require mempro + a real PDB file
# Examples 7-10, 12 require a pre-existing Orient/ folder
```
