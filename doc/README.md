# Building the documentation

There are two types of documentation:
- Fortran code documentation (Doxygen)
- General installation, guide and python API documentation (Sphinx)

To build:
- Fortran documentation: `doxygen Doxyfile`
- General documentation: `make html`


To view documentation:
- open `_build/html/index.html`

## Note
Its possible to build either documentation independently. General documentation links to the Fortran documentation meaning that the Doxygen should be built before Sphinx.

## Dependencies
- Sphinx >= 1.6.7
- Doxygen >= 1.8.6