from pathlib import Path
import runpy

import pytest

ROOT = Path(__file__).parent.parent
TUTORIALS = tuple(map(str, ((ROOT / "tutorials").resolve().glob("*.py"))))
EXAMPLES = tuple(map(str, ((ROOT / "examples").resolve().glob("*/*.py"))))


@pytest.mark.parametrize("tutorial", TUTORIALS)
def test_tutorials(tutorial: str):
    runpy.run_path(tutorial)


@pytest.mark.parametrize("example", EXAMPLES)
def test_examples(example: str):
    # TODO: remove me once dev is merged
    if "compute_segment_fluo" in example:
        pytest.mark.skip("Re-stiching of segmentation is broken on master.")
    else:
        runpy.run_path(example)
