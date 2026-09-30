# Preparation examples

Scripts included in `docs/api/preparation.md` and run once by
`pytest tests/test_preparation.py::TestPreparationExamples`.

```bash
cd tests/preparation_examples
python preparation_example_01.py
```

Uses `protein.pdb` in this folder. `preparation_example_23.py` restores
original residue numbers (gaps and ACE/NME) without PropKa. Complex-structure
tests in `tests/test_preparation.py` use `tests/6RV3_AB.pdb` and `tests/8I5B.pdb`.
