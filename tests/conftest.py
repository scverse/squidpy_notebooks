import gc

import pytest

import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("Agg")


@pytest.fixture(autouse=True)
def clear_figures():
    yield

    plt.close("all")

    gc.collect()
    gc.collect()
    gc.collect()
