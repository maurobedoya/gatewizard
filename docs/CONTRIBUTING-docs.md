# Contributing to the docs

Hand-written Markdown in `docs/api/` plus numbered example scripts. Do not paste script bodies into Markdown.

## New public method

1. Add a stable heading on the module page: `### method_name` (anchor `#method_name`).
2. Short signature, params, returns. Link an existing example or add `tests/<module>_examples/<module>_example_NN.py`.
3. Include the script with a snippet, not a copy:

   ````markdown
   ```python
   --8<-- "tests/analysis_examples/analysis_example_01.py"
   ```
   ````

4. Regenerate the dictionary:

   ```bash
   python scripts/gen_api_dictionary.py
   ```

5. Pytest picks up new `*_example_*.py` files automatically (`tests/example_runner.py`). One run per file.

## Rules for examples

- Paths relative to the example folder or `Path(__file__).parent`.
- Input data in the repo. No machine-local paths.
- Optional tools (GROMACS, FATSLiM, PACKMOL, MemPrO, display) **skip** with a reason.
- Assert public results (keys, shapes, files next to the example).
