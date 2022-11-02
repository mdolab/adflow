from .libadflow import version


def fortranVersion():
    """
    Return the version of the compiled Fortran code

    Returns
    -------
    ver : str
            the value of `git describe --dirty --always --tags` at compile time of the libadflow object.
    """

    ver = version()
    return ver.decode("utf-8").strip()
