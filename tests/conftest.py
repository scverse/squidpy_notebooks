import gc

import pytest

import matplotlib

matplotlib.use("Agg")


@pytest.fixture(scope="session", autouse=True)
def _garbage_collect():
    gc.collect()
