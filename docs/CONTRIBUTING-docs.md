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

## Manual PDF

From the repository root, with `requirements-docs.txt` installed:

```bash
python scripts/build_docs_pdf.py
```

`scripts/build_docs_pdf.sh` runs the same command. The file is `docs/manual/GateWizard-manual.pdf`, which the README links to. The cover uses the black GUI emblem and wordmark (`app-window-light.png` and `wordmark-light.png`). Versions come from `pyproject.toml` and the GUI `package.json`.

Contributors are read from `contributors.yml` in each repository listed in `docs/manual/repositories.yml`. Add a name to that repository's file. The same page also lists GitHub usernames, with profile links, from each repository's public contributor list. Add a repository by listing its path there. The page order is `nav` in `mkdocs.yml`, and `--8<--` snippets are expanded. Pass `-o` for another path.
