from setuptools import setup
import re

__version__ = re.findall(
    r"""__version__ = ["']+([0-9\.]*)["']+""",
    open('adflow/__init__.py').read(),
)[0]

setup(name='adflow',
      version=__version__,


      description="ADflow is a multi-block structured flow solver developed by the MDO Lab at the University of Michigan",
      long_description="""ADflow is a multi-block structured flow solver developed by the MDO Lab at the University of Michigan.
      It solves the compressible Euler, laminar Navier-Stokes and Reynolds-Averaged Navier-Stokes equations.
      ADflow's features include the following:

      - Discrete adjoint implementation
      - "Complexified" code for complex-step derivative verification
      - Massively parallel (both CPU and memory scalable) implementation using MPI.

      ## Documentation
      Please see the [documentation](https://mdolab-adflow.readthedocs-hosted.com/en/latest/) for installation details and API documentation.

      To locally build the documentation, enter the `doc` folder and enter `make html` in terminal.
      You can then view the built documentation in the `_build` folder.

      """,
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
            'numpy>=1.16.4',
            'baseclasses>=1.2.0',
            'mpi4py>=3.0.2',
            'petsc4py>=3.11.0',

      ],
      classifiers=[
        "Operating System :: Linux",
        "Programming Language :: Python, Fortran"]
      )

