#!/usr/bin/env python3

"""Pandoc filter that removes `interpreted-text` from `Code` in Markdown."""
from pandocfilters import Code, toJSONFilter


def _filter(key, value, format, meta):
    if format != "markdown":
        return
    if key != "Code":
        return
    if "interpreted-text" in value[0][1]:
        return Code([value[0][0], [], []], value[1])


if __name__ == "__main__":
    toJSONFilter(_filter)
