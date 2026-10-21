# Contributing to GateWizard

Thanks for helping with the API. The desktop app lives in a separate repository: [gatewizard-gui](https://github.com/franciscoadasme/gatewizard-gui).

Please read the [Code of Conduct](CODE_OF_CONDUCT.md) before participating.

## Where to start

| You want to | Go to |
|-------------|--------|
| Ask how something works | [Discussions](https://github.com/maurobedoya/gatewizard/discussions) (label `api`) |
| Report a bug or request a feature in this library | [Issues](https://github.com/maurobedoya/gatewizard/issues/new/choose) |
| Report a desktop-app bug | [gatewizard-gui issues](https://github.com/franciscoadasme/gatewizard-gui/issues/new/choose) |
| Report a security vulnerability | [SECURITY.md](SECURITY.md) — email [gatewizard.project@gmail.com](mailto:gatewizard.project@gmail.com), not a public issue |

Search existing discussions and issues before opening a new one.

## Set up a development copy

Linux or WSL:

```bash
git clone https://github.com/maurobedoya/gatewizard.git
cd gatewizard
conda env create -f environment.yml
conda activate gatewizard
```

`environment.yml` installs the package in editable mode (`pip install -e ".[full]"`). On macOS, create the environment without `cudatoolkit`.

Check the install:

```bash
python -c "import gatewizard; print(gatewizard.__version__)"
```

## Tests

Run the suite from the repository root:

```bash
pytest tests/
```

Add or update tests for behavior you change. Do not open a pull request with failing tests. Optional tools (GROMACS, FATSLiM, PACKMOL, MemPrO, a display) should skip with a reason when they are missing, not fail.

Details are in [tests/README.md](tests/README.md).

## Documentation

API docs are Markdown in `docs/`. Follow [docs/CONTRIBUTING-docs.md](docs/CONTRIBUTING-docs.md): link example scripts instead of pasting them, and regenerate the dictionary when you add a public method.

## Pull requests

1. Fork the repository and create a branch from `main`.
2. Keep the change focused on one fix or feature.
3. Describe what changed and why in the pull request.
4. Run the relevant tests.

For a large change, open an issue or a discussion first so we can agree on the approach.

This project is licensed under the [MIT License](LICENSE). Contributions are accepted under the same license.
