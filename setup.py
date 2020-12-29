from setuptools import setup
import re
import os

__version__ = re.findall(
    r"""__version__ = ["']+([0-9\.]*)["']+""",
    open('adflow/__init__.py').read(),
)[0]

this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(name='adflow',
      version=__version__,
      description="ADflow is a multi-block structured flow solver developed by the MDO Lab at the University of Michigan",
      long_description=long_description,
      long_description_content_type="text/markdown",
      keywords='RANS adjoint fast optimization',
      author='',
      author_email='',
      url='https://github.com/mdolab/adflow',
      license='LGPL version 2.1',
      packages=[
          'adflow',
      ],
      package_data={
          'adflow': ['*.so']
      },
      install_requires=[
            'numpy>=1.16',
            'mdolab-baseclasses>=1.2',
            'mpi4py>=3.0',
            'petsc4py>=3.11',

      ],
      extras_require={
          'testing': ["parameterized", "testflo"]
      },
      classifiers=[
        "Operating System :: Linux",
        "Programming Language :: Python, Fortran"]
      )
