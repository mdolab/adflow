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
]

defaultAeroDVs = ["alpha", "beta", "mach", "P", "T", "xRef", "yRef", "zRef"]

adflowDefOpts = CaseInsensitiveDict(
    {
        # Common Parameters
        "outputDirectory": "../output_files",
        # Physics Parameters
        "smoother": "Runge-Kutta",
        "equationType": "Euler",
        "sepSensorModel": "heaviside",
    }
)

IDWarpDefOpts = CaseInsensitiveDict({})
