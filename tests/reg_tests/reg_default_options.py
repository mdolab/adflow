"""these are a set of default options for the regression tests"""

from baseclasses.utils import CaseInsensitiveDict


defaultFuncList = [
    "lift",
    "drag",
    "cl",
    "cd",
    "fx",
    "fy",
    "fz",
    "cfx",
    "cfy",
    "cfz",
    "mx",
    "my",
    "mz",
    "cmx",
    "cmy",
    "cmz",
    "sepsensor",
    "sepsensoravgx",
    "sepsensoravgy",
    "sepsensoravgz",
    "cofxx",
    "cofxy",
    "cofxz",
    "cofyx",
    "cofyy",
    "cofyz",
    "cofzx",
    "cofzy",
    "cofzz",
]

defaultAeroDVs = ["alpha", "beta", "mach", "P", "T", "xRef", "yRef", "zRef"]

adflowDefOpts = CaseInsensitiveDict(
    {
        # Common Parameters
        "outputDirectory": "../output_files",
        # Physics Parameters
        "smoother": "Runge-Kutta",
        "equationType": "Euler",
    }
)

IDWarpDefOpts = CaseInsensitiveDict({})
