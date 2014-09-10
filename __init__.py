# -*- coding: utf-8 -*-
"""
tmm - A transfer-matrix method optics package, for calculating
reflection, transmission, absorption, and other relevant aspects of thin
and thick multilayer (or single-layer) films.

Written by Steven Byrnes, http://sjbyrnes.com . Package lives at
https://pypi.python.org/pypi/tmm

Released under MIT license (Expat).

For detailed derivations of the formulas, see manual.pdf. (If you can't
find it, go to http://sjbyrnes.com/fresnel_manual.pdf ) For the
various functions and their explanations and syntaxes, browse tmm_core,
particularly the docstrings of all the functions. Most important of
these are:

tmm.coh_tmm(...) -- the transfer-matrix-method calculation in the coherent
case (i.e. thin films)

tmm.inc_tmm(...) -- the transfer-matrix-method calculation in the incoherent
case (i.e. films tens or hundreds of wavelengths thick, or whose
thickness is not very uniform.

In tmm.examples you will find some sample calculations and plots.

In tmm.tests you will find various tests that the package is coded
correctly.

In tmm.colors you will find extra functions for calculating the color of a
multilayer thin film under reflected light.
"""
from __future__ import division, print_function, absolute_import

from .tmm_core import *
