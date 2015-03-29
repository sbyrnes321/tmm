# -*- coding: utf-8 -*-
"""
Functions to calculate the color of a multilayer thin film under reflected
light. A perfect mirror will look white, because we imagine seeing the white
light source ("illuminant") reflected in it. A half-reflective mirror will be
gray, a non-reflective surface will be black, etc. See tmm.examples.sample5()
for a few example calculations for how this is used.

For functions that require an illuminant, the most common choice would be to
use colorpy.illuminants.get_illuminant_D65(), which approximates a phase of
natural daylight. See http://en.wikipedia.org/wiki/Illuminant_D65 .
"""

from __future__ import division, print_function, absolute_import

from numpy import arange, array

import numpy as np

from .tmm_core import coh_tmm

try:
    import colorpy
    import colorpy.illuminants
    import colorpy.ciexyz
except ImportError:
    print('Warning: Colorpy not detected (or perhaps an error occurred when',
          'loading it). Film color calculations (in tmm.color)',
          'will not work. http://pypi.python.org/pypi/colorpy')

inf = float('inf')

def calc_reflectances(n_fn_list, d_list, th_0, pol='s', spectral_range='narrow'):
    """
    Calculate the reflection spectrum of a thin-film stack.
    
    n_fn_list[m] should be a function that inputs wavelength in nm and
    outputs refractive index of the m'th layer. In other words,
    n_fn_list[2](456) == 1.53 + 0.4j mans that layer #2 has a refractive index
    of 1.53 + 0.4j at 456nm. These functions could be defined with
    scipy.interpolate.interp1d() for example.
    
    pol, d_list and th_0 are defined as in tmm.coh_tmm ... but d_list
    MUST be in units of nanometers
    
    spectral_range can be 'full' if all the functions in n_fn_list can take
    wavelength arguments between 360-830nm; or 'narrow' if some or all require
    arguments only in the range 400-700nm. The wavelengths outside the
    'narrow' range make only a tiny difference to the color, because they are
    almost invisible to the eye. If spectral_range is 'narrow', then the n(400)
    values are used for 360-400 and n(700) for 700-830nm
    
    Returns a 2-column array where the first column is wavelength in nm
    (360,361,362,...,830) and the second column is reflectivity (from 0
    to 1, where 1 is a perfect mirror). This range is chosen to be
    consistent with colorpy.illuminants. See  colorpy.ciexyz.start_wl_nm etc.
    """
    
    lam_vac_list = arange(360, 831)
    
    num_layers = len(n_fn_list)
    
    def extend_spectral_range(n_fn):
        """
        Starting with a narrow-spectrum refractive index function
        n_fn(wavelength), create then return the corresponding full-spectrum
        refractive index function
        """
        def extended_n_fn(lam):
            if lam < 400:
                return n_fn(400)
            elif lam > 700:
                return n_fn(700)
            else:
                return n_fn(lam)
        return extended_n_fn
    
    if spectral_range == 'narrow':
        n_fn_list = [extend_spectral_range(n_fn) for n_fn in n_fn_list]
    
    final_answer = []
    
    for lam_vac in lam_vac_list:
        n_list = [n_fn_list[i](lam_vac) for i in range(num_layers)]
        R = coh_tmm(pol, n_list, d_list, th_0, lam_vac)['R']
        final_answer.append([lam_vac,R])
    final_answer = array(final_answer)

    return final_answer

def calc_spectrum(reflectances, illuminant):
    """
    * reflectances is the output of calc_reflec_spec()
    * illuminant is a 2D numpy arrays, with one row for each wavelength,
      with the first column holding the wavelength in nm, and the
      second column the intensity. This is the form returned by the
      functions in colorpy.illuminants.  It is normally assumed that
      illuminant is normalized so that Y=1.
    """
    #Both colorpy.illuminants and calc_reflec_spec should go from
    #colorpy.ciexyz.start_wl_nm etc, so they should have matching
    #wavelength specifications
    if not np.all(reflectances[:,0] == illuminant[:,0]):
        raise ValueError('Wavelength range is inconsistent...Both should be 360,361,...,830.\n'
        + 'reflectances[0]=' + str(reflectances[0]) + ', reflectances[-1]=' + str(reflectances[-1])
        + '\nilluminant[0]=' + str(illuminant[0]) + ', illuminant[-1]=' + str(reflectances[-1]))
    
    final_answer = []
    for i,lam in enumerate(reflectances[:,0]):
        final_answer.append([lam, reflectances[i,1] * illuminant[i,1]])
    return array(final_answer)

def calc_color(spectrum, scale=None, show_warnings=True):
    """
    Calculate the color in various representations.
    
    spectrum is the output of calc_spectrum.
    
    scale is the scaling method. Possibilities are:

    * scale=None means don't scale. This is usually what you want, bucause
      the illuminant should be pre-scaled in an appropriate way.
      (Specifically, it's scaled to get Y=1 for a perfect reflector.)
    * scale='Y1' means that the intensity is increased or decreased in
      order to set Y (the luminance) to 1. So you can get white but not gray,
      you can get orange but not brown, etc.
    * scale=0.789 multiplies X,Y,Z by 0.789. Any number > 0 is OK.
    
    Returns a dictionary with rgb, irgb, xy, xyY, and XYZ. Definitions:

    * xy, xyY and XYZ are defined as in
        http://en.wikipedia.org/wiki/CIE_1931_color_space
    * rgb is the linear (i.e., proportional to intensity, not
      gamma-corrected) version of sRGB.
    * irgb is ready-to-display sRGB, i.e. it is clipped to the range 0-1,
      and gamma-corrected, and rounded to three integers in the range 0-255.

    (sRGB is the standard RGB used in modern displays and printers.)
    """
    assert (scale is None or scale == 'Y1'
            or (type(scale) is float and scale > 0))
    XYZ = colorpy.ciexyz.xyz_from_spectrum(spectrum)
    assert min(XYZ) >= 0
    if scale == 'Y1' or type(scale) is float:
        factor = (1.0 / XYZ[1] if scale == 'Y1' else scale)
        XYZ[0] *= factor
        XYZ[1] *= factor
        XYZ[2] *= factor
    X,Y,Z = XYZ
    if show_warnings:
        if Y > 1:
            print('Warning: Oversaturated color! XYZ = ', XYZ)
    xy = [X / (X + Y + Z), Y / (X + Y + Z)]
    xyY = [xy[0], xy[1], Y]
    rgb = colorpy.colormodels.rgb_from_xyz(XYZ)
    irgb = colorpy.colormodels.irgb_from_rgb(rgb)
    return {'xy':xy, 'xyY':xyY, 'XYZ':XYZ, 'rgb':rgb, 'irgb':irgb}
        
def plot_reflectances(reflectances, filename='temp_plot.png', title='Reflectance', ylabel='Fraction reflected'):
    """
    Makes nice colored plot of reflectances. reflectances is the output of
    calc_reflectances(...)
    """
    colorpy.plots.spectrum_plot(reflectances, title, filename, ylabel=ylabel)

def plot_spectrum(spectrum, filename='temp_plot.png', title='Reflected light under illumination', ylabel='Intensity (a.u.)'):
    """
    Makes nice colored plot of the reflected color spectrum you see under a
    certain illuminant. spectrum is the output of
    calc_spectrum(...)
    """
    colorpy.plots.spectrum_plot(spectrum, title, filename, ylabel=ylabel)

