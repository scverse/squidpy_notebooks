from pathlib import Path
import runpy

import pytest

ROOT = Path(__file__).parent.parent
TUTORIALS = tuple((ROOT / "tutorials").resolve().glob("*.py"))
EXAMPLES = tuple((ROOT / "examples").resolve().glob("*/*.py"))


@pytest.mark.parametrize("tutorial", TUTORIALS)
def test_tutorials(tutorial):
    runpy.run_path(tutorial)


@pytest.mark.parametrize("example", EXAMPLES)
def test_examples(example):
    runpy.run_path(example)
