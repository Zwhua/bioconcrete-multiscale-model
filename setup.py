"""Compatibility shim for legacy pip editable installs on Python 3.8.

Project metadata remains canonical in ``pyproject.toml``; modern installers use
PEP 517/660 while older pip can fall back to ``setup.py develop``.
"""

from setuptools import setup

setup()
