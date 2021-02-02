from abc import ABCMeta
from glob import glob
from logging import warning
from pathlib import Path
from textwrap import indent
import os

import matplotlib

matplotlib.use("Agg")
_test_template = """
def _test(self):
{code}
"""


class MetaTester(ABCMeta):
    def __new__(cls, clsname: str, superclasses, attributedict):
        print(clsname, attributedict)
        if not clsname.startswith("Test"):
            raise ValueError(f"Class `{clsname}` does not start with `Test`.")

        dirname = Path(clsname[4:].lower())
        if not dirname.is_dir():
            raise ValueError(f"Path `{dirname}` is not a directory.")

        files = [p for p in glob(str(dirname / "**/*.py"), recursive=True) if Path(p).is_file()]
        if not len(files):
            warning(f"Class `{clsname}` does not contain any automatically generated tests.")

        for file in files:
            tmp = {}
            with open(file) as fin:
                exec(
                    compile(
                        _test_template.format(code=indent(fin.read(), "    ")),
                        "",
                        "exec",
                    ),
                    tmp,
                )
            test_key = "test_" + file.rstrip(".py").replace(os.sep, "_").lstrip("_")
            if test_key in attributedict:
                raise KeyError(f"Test `{test_key}` is already present: `{attributedict.keys()}`.")
            attributedict[test_key] = tmp["_test"]

        return super().__new__(cls, clsname, superclasses, attributedict)
