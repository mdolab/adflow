from sphinx_mdolab_theme.config import *
import subprocess

# -- Path setup --------------------------------------------------------------
# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import sys

sys.path.insert(0, os.path.abspath("../"))
sys.path.insert(0, os.path.abspath("./"))  # to import custom Sphinx extension

from ext.costfunctionslist import TEMP_FILE

# build doxygen
read_the_docs_build = os.environ.get("READTHEDOCS", None) == "True"
if read_the_docs_build:
    subprocess.call("doxygen", shell=True)

# -- Project information -----------------------------------------------------
project = "ADflow"

# -- General configuration ---------------------------------------------------
# Built-in Sphinx extensions are already contained in the imported variable
# here we add external extensions, which must also be added to requirements.txt
# so that RTD can import and use them
extensions.extend(["numpydoc", "ext.costfunctionslist"])

html_static_path = ['_static']

# mock import for autodoc
autodoc_mock_imports = ["numpy", "mpi4py", "baseclasses"]

# bibtex sources
bibtex_bibfiles.extend(["citations.bib"])

# Ignore the temp file from the cost functions list extension
exclude_patterns.extend([TEMP_FILE])
