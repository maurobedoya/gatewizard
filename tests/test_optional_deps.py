"""Tests for optional dependency utilities."""

import sys
from importlib import metadata
from pathlib import Path
from unittest.mock import patch

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from gatewizard import __version__ as GATEWIZARD_VERSION
from gatewizard.utils.optional_deps import (
    EXTERNAL_TOOL_REGISTRY,
    get_dependency_versions,
    get_external_tool_versions,
    get_optional_dependencies_status,
    get_package_version,
    parse_engine_variant,
    parse_tool_version,
    resolve_executable,
)


class TestGetPackageVersion:
    @patch("gatewizard.utils.optional_deps.metadata.version", return_value="9.9.9")
    def test_uses_distribution_metadata(self, mock_version):
        assert get_package_version("numpy", distribution_name="numpy") == "9.9.9"
        mock_version.assert_called_once_with("numpy")

    @patch(
        "gatewizard.utils.optional_deps.metadata.version",
        side_effect=metadata.PackageNotFoundError,
    )
    @patch("gatewizard.utils.optional_deps.safe_import")
    def test_falls_back_to_module_version(self, mock_import, _mock_metadata):
        mock_import.return_value = type("Module", (), {"__version__": "1.2.3"})()
        assert get_package_version("fakepkg") == "1.2.3"

    @patch(
        "gatewizard.utils.optional_deps.metadata.version",
        side_effect=metadata.PackageNotFoundError,
    )
    @patch("gatewizard.utils.optional_deps.safe_import", return_value=None)
    def test_returns_none_when_unavailable(self, _mock_import, _mock_metadata):
        assert get_package_version("missingpkg") is None


class TestParseToolVersion:
    def test_parse_gromacs_version(self):
        text = ":-) GROMACS - gmx, 2026.2 (-:"
        assert parse_tool_version(text, "gromacs") == "2026.2"

    def test_parse_namd_skips_charmpp_banner(self):
        text = (
            "Charm++> No provisioning arguments specified. Running with a single PE.\n"
            "Charm++> Charm++ version 7.5.0\n"
            "Info: NAMD 3.0.1 for Linux-x86_64-multicore"
        )
        assert parse_tool_version(text, "namd") == "3.0.1"

    def test_parse_namd_uses_canonical_tool_name_for_namd3(self):
        text = (
            "Charm++> Charm++ version 7.5.0\n"
            "Info: NAMD 3.0.2 for Linux-x86_64-multicore"
        )
        assert parse_tool_version(text, "namd3") == "3.0.2"

    def test_namd_version_from_install_path(self):
        from gatewizard.utils.optional_deps import _version_from_install_path

        path = "/mnt/c/software/namd/NAMD_3.0.2_Linux-x86_64-multicore/namd3"
        assert _version_from_install_path(path, "namd3") == "3.0.2"

    def test_rejects_teLeap_error_as_version(self):
        text = "/home/user/mamba-env/bin/teLeap: invalid option -- 'L'"
        assert parse_tool_version(text, "tleap") is None

    def test_rejects_single_digit_version(self):
        from gatewizard.utils import optional_deps

        assert optional_deps._is_plausible_version("1") is False
        assert optional_deps._is_plausible_version("24.8") is True

    def test_parse_packmol_version(self):
        text = "Version 20.14.3"
        assert parse_tool_version(text, "packmol") == "20.14.3"

    def test_ignores_jax_warning_for_mempro(self):
        text = (
            "WARNING:2026-06-05 18:35:07,599:jax._src.xla_bridge:909: "
            "An NVIDIA GPU may be present on this machine"
        )
        assert parse_tool_version(text, "mempro") is None


class TestParseEngineVariant:
    def test_gromacs_cuda_from_gpu_support_line(self):
        text = "GROMACS version:    2024.4\nGPU support:             CUDA\n"
        assert parse_engine_variant(text, "gromacs") == "CUDA"

    def test_gromacs_cpu_when_gpu_disabled(self):
        text = "GROMACS version:    2024.4\nGPU support:             disabled\n"
        assert parse_engine_variant(text, "gromacs") == "CPU"

    def test_gromacs_cuda_from_build_string(self):
        text = "gromacs-2024.4-nompi_cuda_h123_0"
        assert parse_engine_variant(text, "gromacs", "/env/bin/gmx") == "CUDA"

    def test_namd_cuda_from_install_path(self):
        path = "/opt/NAMD_3.0.1_Linux-x86_64-multicore-CUDA/namd3"
        assert parse_engine_variant("", "namd", path) == "CUDA"

    def test_namd_cpu_from_multicore_path(self):
        path = "/opt/NAMD_3.0.1_Linux-x86_64-multicore/namd3"
        assert parse_engine_variant("", "namd", path) == "CPU"

    def test_parse_ambertools_from_version_file(self):
        text = "24.8"
        assert parse_tool_version(text, "ambertools") == "24.8"


class TestGetDependencyVersions:
    def test_includes_gatewizard_version(self):
        report = get_dependency_versions(include_platform=False)
        assert "dependencies" in report
        assert report["dependencies"]["gatewizard"]["version"] == GATEWIZARD_VERSION
        assert report["dependencies"]["gatewizard"]["available"] is True

    def test_required_packages_present(self):
        report = get_dependency_versions(include_optional=False, include_platform=False)
        deps = report["dependencies"]
        assert "numpy" in deps
        assert "MDAnalysis" in deps
        assert deps["numpy"]["install_group"] == "core"
        assert "parmed" not in deps

    def test_optional_packages_have_install_groups(self):
        report = get_dependency_versions(include_platform=False)
        deps = report["dependencies"]
        assert deps["parmed"]["install_group"] == "md"
        assert deps["openmm"]["install_group"] == "md"
        assert deps["mempro"]["install_group"] == "orientation"

    def test_platform_metadata_when_requested(self):
        report = get_dependency_versions(include_platform=True)
        assert "platform" in report
        assert "python_version" in report["platform"]

    def test_external_tools_when_requested(self):
        report = get_dependency_versions(include_external_tools=True)
        assert "executables" in report
        names = {item["name"] for item in report["executables"]}
        assert names == set(EXTERNAL_TOOL_REGISTRY.keys())

    @pytest.mark.parametrize("package", ["numpy", "MDAnalysis"])
    def test_installed_required_packages_have_versions(self, package):
        report = get_dependency_versions(include_platform=False)
        info = report["dependencies"][package]
        assert info["available"] is True
        assert info["version"]


class TestGetOptionalDependenciesStatus:
    def test_only_reports_optional_packages(self):
        status = get_optional_dependencies_status()
        assert "parmed" in status
        assert "openmm" in status
        assert "mempro" in status
        assert "numpy" not in status

    def test_entries_include_version_field(self):
        status = get_optional_dependencies_status()
        for info in status.values():
            assert "available" in info
            assert "description" in info
            assert "version" in info
            assert "install_group" in info


class TestExternalTools:
    @patch("gatewizard.utils.optional_deps._resolve_ambertools_version", return_value="24.8")
    @patch("gatewizard.utils.optional_deps.resolve_executable", return_value="/usr/bin/packmol-memgen")
    @patch("gatewizard.utils.optional_deps.probe_executable_version", return_value=None)
    def test_packmol_memgen_falls_back_to_ambertools(
        self, _mock_probe, _mock_resolve, _mock_ambertools
    ):
        tools = get_external_tool_versions()
        memgen = next(item for item in tools if item["name"] == "packmol-memgen")
        assert memgen["version"] == "24.8"

    @patch("gatewizard.utils.optional_deps._resolve_ambertools_version", return_value="24.8")
    def test_ambertools_uses_conda_or_env_resolvers(self, mock_resolve):
        tools = get_external_tool_versions()
        amber = next(item for item in tools if item["name"] == "ambertools")
        assert amber["version"] == "24.8"
        mock_resolve.assert_called_once()

    @patch("gatewizard.utils.optional_deps.probe_executable_version", return_value=None)
    @patch(
        "gatewizard.utils.optional_deps.resolve_executable",
        return_value="/mnt/c/software/namd/NAMD_3.0.2_Linux-x86_64-multicore/namd3",
    )
    def test_namd_falls_back_to_install_path(self, _mock_resolve, _mock_probe):
        tools = get_external_tool_versions()
        namd = next(item for item in tools if item["name"] == "namd")
        assert namd["version"] == "3.0.2"

    @patch("gatewizard.utils.optional_deps.resolve_executable", return_value="/usr/bin/packmol")
    @patch(
        "gatewizard.utils.optional_deps.probe_executable_version",
        return_value="20.14.3",
    )
    def test_get_external_tool_versions_uses_probe(self, mock_probe, mock_resolve):
        tools = get_external_tool_versions()
        packmol = next(item for item in tools if item["name"] == "packmol")
        assert packmol["available"] is True
        assert packmol["version"] == "20.14.3"
        mock_probe.assert_called()

    @patch("gatewizard.utils.optional_deps.get_package_version", return_value="0.1.0")
    @patch("gatewizard.utils.optional_deps.resolve_executable", return_value="/usr/bin/mempro")
    @patch("gatewizard.utils.optional_deps.probe_executable_version")
    def test_mempro_prefers_python_package_version(
        self, mock_probe, _mock_resolve, _mock_pkg_version
    ):
        tools = get_external_tool_versions()
        mempro = next(item for item in tools if item["name"] == "mempro")
        assert mempro["version"] == "0.1.0"
        mock_probe.assert_not_called()

    @patch("gatewizard.utils.optional_deps.shutil.which", return_value="/usr/bin/gmx")
    def test_resolve_executable(self, _mock_which):
        assert resolve_executable(("gmx",)) == "/usr/bin/gmx"
