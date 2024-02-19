# Analytical benchmark for annulus and spherical shell

- The analytical solutions for freeslip-freeslip are taken from https://doi.org/10.5194/gmd-14-1899-2021, through the `assess` python package which is a requirement and must be installed in python through pip(3).
    - The solutions for `mixed` case were derived and available through a fork of `assess` package here https://github.com/suganth1997/assess
    - So the python command to install this version of `assess` package must be
        - `pip3 install git+https://github.com/suganth1997/assess`
    - The correct version of Python to which the pip3 installs the package must be linked, which should happen automatically if only one Python version is on the system.
        - https://cmake.org/cmake/help/latest/module/FindPython3.html
        - Maybe use `pyenv` to manage the versions of python in the system
            - https://github.com/pyenv/pyenv
<br/>

- The type of forcing can be controlled with the `delta` parameter in the `.prm` files

- The type of coundary condition can also be controlled with the `freeslip` or `mixed` parameter in the `.prm` files
    - `freeslip` - Freeslip BC on CMB and surface
    - `mixed` - Freeslip BC on CMB and No-slip on surface
    - If both are set to false, `no-slip` is used for CMB and surface
