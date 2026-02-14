from pathlib import Path


def test_runtime_dependencies_exclude_dev_and_test_packages():
    pyproject = Path(__file__).resolve().parents[1] / "pyproject.toml"
    content = pyproject.read_text(encoding="utf-8")

    project_block = content.split("[project.optional-dependencies]")[0]
    assert '"pytest>=' not in project_block
    assert '"requests>=' not in project_block
    assert '"ipykernel>=' not in project_block


def test_python_requirement_matches_supported_dependency_range():
    pyproject = Path(__file__).resolve().parents[1] / "pyproject.toml"
    content = pyproject.read_text(encoding="utf-8")
    assert 'requires-python = ">=3.8,<3.9"' in content
