from pathlib import Path
import runpy

import pytest

ROOT = Path(__file__).parent.parent
TUTORIALS = tuple(map(str, ((ROOT / "tutorials").resolve().glob("*.py"))))
EXAMPLES = tuple(map(str, ((ROOT / "examples").resolve().glob("*/*.py"))))


@pytest.mark.parametrize("tutorial", TUTORIALS)
def test_tutorials(tutorial):
    runpy.run_path(tutorial)


@pytest.mark.parametrize("example", EXAMPLES)
def test_examples(example):
    runpy.run_path(example)
