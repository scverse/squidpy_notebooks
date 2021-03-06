from typing import Any, Dict
import os
import re
import stat
import codecs

from sphinx_gallery.utils import _replace_md5
from sphinx_gallery.binder import gen_binder_rst, check_binder_conf
from sphinx_gallery.gen_rst import CODE_DOWNLOAD, TIMING_CONTENT, replace_py_ipynb

EXAMPLE_HEADER = """
.. DO NOT EDIT.
.. THIS FILE WAS AUTOMATICALLY GENERATED BY SPHINX-GALLERY.
.. TO MAKE CHANGES, EDIT THE SOURCE PYTHON FILE:
.. "{0}"
.. LINE NUMBERS ARE GIVEN BELOW.

.. only:: html
{2}
.. rst-class:: sphx-glr-example-title

.. _sphx_glr_{1}:
"""


def save_rst_example(
    example_rst: str, example_file: str, time_elapsed: float, memory_used: float, gallery_conf: Dict[str, Any]
) -> None:
    """
    Saves the rst notebook to example_file including header & footer.

    Parameters
    ----------
    example_rst
        RST containing the executed file content.
    example_file
        Filename with full path of python example file in documentation folder.
    time_elapsed
        Time elapsed in seconds while executing file.
    memory_used
        Additional memory used during the run.
    gallery_conf
        Sphinx-Gallery configuration dictionary.

    Returns
    -------
    None
    """
    example_fname = os.path.relpath(example_file, gallery_conf["src_dir"])
    ref_fname = example_fname.replace(os.path.sep, "_")

    binder_conf = check_binder_conf(gallery_conf.get("binder"))
    binder_badge_rst = ""
    if len(binder_conf) > 0:
        binder_badge_rst = gen_binder_rst(example_file, binder_conf, gallery_conf)

    example_rst = EXAMPLE_HEADER.format(example_fname, ref_fname, binder_badge_rst) + example_rst

    if time_elapsed >= gallery_conf["min_reported_time"]:
        time_m, time_s = divmod(time_elapsed, 60)
        example_rst += TIMING_CONTENT.format(time_m, time_s)
    if gallery_conf["show_memory"]:
        example_rst += f"**Estimated memory usage:** {memory_used: .0f} MB\n\n"

    # Generate a binder URL if specified
    fname = os.path.basename(example_file)
    example_rst += CODE_DOWNLOAD.format(fname, replace_py_ipynb(fname), "", ref_fname)

    write_file_new = re.sub(r"\.py$", ".rst.new", example_file)
    with codecs.open(write_file_new, "w", encoding="utf-8") as f:
        f.write(example_rst)
    # make it read-only so that people don't try to edit it
    mode = os.stat(write_file_new).st_mode
    ro_mask = 0x777 ^ (stat.S_IWRITE | stat.S_IWGRP | stat.S_IWOTH)
    os.chmod(write_file_new, mode & ro_mask)
    # in case it wasn't in our pattern, only replace the file if it's
    # still stale.
    _replace_md5(write_file_new, mode="t")