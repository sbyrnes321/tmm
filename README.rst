Author homepage: http://sjbyrnes.com

Package home at PyPI: http://pypi.python.org/pypi/tmm

Package home at Github: https://github.com/sbyrnes321/tmm

This is a group of programs written in Python / NumPy for simulating light propagation in planar multilayer thin films, including the effects of multiple internal reflections and interference, using the "Transfer Matrix Method". It can also simulate combinations of thin and thick films (e.g. a thick piece of glass with a multi-layer antireflection coating on one side and a mirror on the other side), or purely thick films.

In addition to calculating how much light is transmitted and reflected, the program can calculate, at any given point in the structure, how much light is being absorbed there. This is a very important feature for solar-cell modeling, for example.

It can also calculate the parameters measured in ellipsometry. It can also calculate the RGB or xyY color of a multilayer thin film (this requires colorpy, https://pypi.python.org/pypi/colorpy).

For more information, see manual.pdf . Get it directly at http://sjbyrnes.com/fresnel_manual.pdf . For the list of all functions and how to call them, browse the source code or go to https://pythonhosted.org/tmm/

Tested in Python 2.7 and 3.3. (If you want to do color calculations in Python 3, you need to use the Python-3-compatible version of colorpy `here <https://github.com/fish2000/ColorPy/>`_.)
