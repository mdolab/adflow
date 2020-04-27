# ADflow
[![Build Status](https://travis-ci.com/mdolab/adflow.svg?branch=master)](https://travis-ci.com/mdolab/adflow) 
[![Documentation Status](https://readthedocs.com/projects/mdolab-adflow/badge/?version=latest)](https://mdolab-adflow.readthedocs-hosted.com/?badge=latest)

ADflow is a multi-block structured flow solver developed by the MDO Lab at the University of Michigan.
It solves the compressible Euler, laminar Navier-Stokes and Reynolds-Averaged Navier-Stokes equations.
ADflow's features include the following:

- Discrete adjoint implementation
- "Complexified" code for complex-step derivative verification
- Massively parallel (both CPU and memory scalable) implementation using MPI.

## Documentation
Please see the [documentation](https://mdolab-adflow.readthedocs-hosted.com/en/latest/) for installation details and API documentation.

To locally build the documentation, enter the `doc` folder and enter `make html` in terminal.
You can then view the built documentation in the `_build` folder.

## Citation

## License
Copyright 2019 MDO Lab

Distributed using the GNU Lesser General Public License (LGPL), verstion 2.1; see
the LICENSE file for details.
