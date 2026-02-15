from pathlib import Path
import re


def test_runtime_dependencies_exclude_dev_and_test_packages():
    pyproject = Path(__file__).resolve().parents[1] / "pyproject.toml"
    content = pyproject.read_text(encoding="utf-8")

    assert "[project.optional-dependencies]" in content
    project_block = content.split("[project.optional-dependencies]", 1)[0]
    assert '"pytest>=' not in project_block
    assert '"requests>=' not in project_block
    assert '"ipykernel>=' not in project_block


def test_python_requirement_matches_supported_dependency_range():
    pyproject = Path(__file__).resolve().parents[1] / "pyproject.toml"
    content = pyproject.read_text(encoding="utf-8")
    match = re.search(r'^requires-python = "([^"]+)"$', content, re.MULTILINE)
    assert match is not None
    assert match.group(1) == ">=3.8,<3.9"


def test_runtime_dependency_uses_triangle_not_tri():
    pyproject = Path(__file__).resolve().parents[1] / "pyproject.toml"
    setup_py = Path(__file__).resolve().parents[1] / "setup.py"
    content = pyproject.read_text(encoding="utf-8")
    setup_content = setup_py.read_text(encoding="utf-8")

    assert '"triangle"' in content
    assert '"tri @ git+https://github.com/bmmeijers/tri/"' not in content
    assert 'REQUIREMENTS = ["geompreds", "triangle"]' in setup_content
    assert 'REQUIREMENTS = ["geompreds", "tri"]' not in setup_content
