# -*- coding: utf-8 -*-
"""
Tests to ensure tmm package was coded correctly. Use run_all() to
run them all in order.
"""

from __future__ import division, print_function, absolute_import
import unittest
from tmm import (
    coh_tmm,
    inc_tmm,
    ellips,
    position_resolved,
    absorp_in_each_layer,
    snell,
    absorp_analytic_fn,
    interface_r,
    inc_absorp_in_each_layer,
    interface_R,
    interface_T,
    power_entering_from_r,
)

import numpy as np

# "5 * degree" is 5 degrees expressed in radians
# "1.2 / degree" is 1.2 radians expressed in degrees
degree = np.pi / 180


class TestBasic(unittest.TestCase):
    n_list = [1, 2 + 4j, 3 + 0.3j, 1 + 0.1j]
    d_list = [np.inf, 2, 3, np.inf]
    th_0 = 0.1
    lam_vac = 100

    def test_coh_stack_rt(self):
        """
        Compare with program I wrote previously in Mathematica. Also confirms
        that I don't accidentally mess up the program by editing.
        """
        s_data = coh_tmm("s", self.n_list, self.d_list, self.th_0, self.lam_vac)

        np.testing.assert_allclose(
            s_data["r"], -0.60331226568845775 - 0.093522181653632019j
        )
        np.testing.assert_allclose(
            s_data["t"], 0.44429533471192989 + 0.16921936169383078j
        )
        np.testing.assert_allclose(s_data["R"], 0.37273208839139516)
        np.testing.assert_allclose(s_data["T"], 0.22604491247079261)

        p_data = coh_tmm("p", self.n_list, self.d_list, self.th_0, self.lam_vac)

        np.testing.assert_allclose(
            p_data["r"], 0.60102654255772481 + 0.094489146845323682j
        )
        np.testing.assert_allclose(
            p_data["t"], 0.4461816467503148 + 0.17061408427088917j
        )
        np.testing.assert_allclose(p_data["R"], 0.37016110373044969)
        np.testing.assert_allclose(p_data["T"], 0.22824374314132009)

    def test_coh_stack_ellipso(self):
        ellips_data = ellips(self.n_list, self.d_list, self.th_0, self.lam_vac)
        np.testing.assert_allclose(ellips_data["psi"], 0.78366777347038352)
        np.testing.assert_allclose(ellips_data["Delta"], 0.0021460774404193292)


class TestPositionResolved(unittest.TestCase):
    """
    Compare with program I wrote previously in Mathematica. Also, various
    consistency checks.
    """

    d_list = [np.inf, 100, 300, np.inf]  # in nm
    n_list = [1, 2.2 + 0.2j, 3.3 + 0.3j, 1]
    th_0 = np.pi / 4
    lam_vac = 400

    def test_poynting(self):
        layer = 1
        dist = 37

        pol = "p"
        coh_tmm_data = coh_tmm(pol, self.n_list, self.d_list, self.th_0, self.lam_vac)
        np.testing.assert_allclose(
            coh_tmm_data["kz_list"][1], 0.0327410685922732 + 0.003315885921866465j
        )
        data = position_resolved(layer, dist, coh_tmm_data)
        np.testing.assert_allclose(data["poyn"], 0.7094950598055798)
        np.testing.assert_allclose(data["absor"], 0.005135049118053356)
        np.testing.assert_allclose(1, sum(absorp_in_each_layer(coh_tmm_data)))

        pol = "s"
        coh_tmm_data = coh_tmm(pol, self.n_list, self.d_list, self.th_0, self.lam_vac)
        np.testing.assert_allclose(
            coh_tmm_data["kz_list"][1], 0.0327410685922732 + 0.003315885921866465j
        )
        data = position_resolved(layer, dist, coh_tmm_data)
        np.testing.assert_allclose(data["poyn"], 0.5422594735025152)
        np.testing.assert_allclose(data["absor"], 0.004041912286816303)
        np.testing.assert_allclose(1, sum(absorp_in_each_layer(coh_tmm_data)))

    def test_poynting_derivative(self):
        # Poynting vector derivative should equal absorption
        layer = 1
        dist = 37

        for pol in ["s", "p"]:
            coh_tmm_data = coh_tmm(
                pol, self.n_list, self.d_list, self.th_0, self.lam_vac
            )
            data1 = position_resolved(layer, dist, coh_tmm_data)
            data2 = position_resolved(layer, dist + 0.001, coh_tmm_data)
            np.testing.assert_allclose(
                (data1["absor"] + data2["absor"]) / 2,
                (data1["poyn"] - data2["poyn"]) / 0.001,
            )

    def test_poynting_vector_end(self):
        # Poynting vector at end should equal T
        layer = 2
        dist = 300

        for pol in ["s", "p"]:
            coh_tmm_data = coh_tmm(
                pol, self.n_list, self.d_list, self.th_0, self.lam_vac
            )
            data = position_resolved(layer, dist, coh_tmm_data)
            np.testing.assert_allclose(data["poyn"], coh_tmm_data["T"])

    def test_poynting_vector_start(self):
        # Poynting vector at start should equal power_entering
        layer = 1
        dist = 0

        for pol in ["s", "p"]:
            coh_tmm_data = coh_tmm(
                pol, self.n_list, self.d_list, self.th_0, self.lam_vac
            )
            data = position_resolved(layer, dist, coh_tmm_data)
            np.testing.assert_allclose(data["poyn"], coh_tmm_data["power_entering"])

    def test_poynting_vector_continuity(self):
        # Poynting vector should be continuous
        for pol in ["s", "p"]:
            layer = 1
            dist = 100
            coh_tmm_data = coh_tmm(
                pol, self.n_list, self.d_list, self.th_0, self.lam_vac
            )
            data = position_resolved(layer, dist, coh_tmm_data)
            poyn1 = data["poyn"]
            layer = 2
            dist = 0
            coh_tmm_data = coh_tmm(
                pol, self.n_list, self.d_list, self.th_0, self.lam_vac
            )
            data = position_resolved(layer, dist, coh_tmm_data)
            poyn2 = data["poyn"]
            np.testing.assert_allclose(poyn1, poyn2)


class TestPositionResolved2(unittest.TestCase):
    """
    Similar to TestPositionResolved, but with initial and final medium
    having a complex refractive index.
    """

    d_list = [np.inf, 100, 300, np.inf]  # in nm
    # "00" is before the 0'th layer. This is easy way to generate th0, ensuring
    # that n0*sin(th0) is real.
    n00 = 1
    th00 = np.pi / 4
    n0 = 1 + 0.1j
    th_0 = snell(n00, n0, th00)
    n_list = [n0, 2.2 + 0.2j, 3.3 + 0.3j, 1 + 0.4j]
    lam_vac = 400

    def test_poynting(self):
        layer = 1
        dist = 37

        for pol in ["s", "p"]:
            coh_tmm_data = coh_tmm(
                pol, self.n_list, self.d_list, self.th_0, self.lam_vac
            )
            data = position_resolved(layer, dist, coh_tmm_data)
            np.testing.assert_allclose(1, sum(absorp_in_each_layer(coh_tmm_data)))

    def test_poynting_derivative(self):
        # Poynting vector derivative should equal absorption
        layer = 1
        dist = 37

        for pol in ["s", "p"]:
            coh_tmm_data = coh_tmm(
                pol, self.n_list, self.d_list, self.th_0, self.lam_vac
            )
            data1 = position_resolved(layer, dist, coh_tmm_data)
            data2 = position_resolved(layer, dist + 0.001, coh_tmm_data)
            np.testing.assert_allclose(
                (data1["absor"] + data2["absor"]) / 2,
                (data1["poyn"] - data2["poyn"]) / 0.001,
            )

    def test_poynting_vector_end(self):
        # Poynting vector at end should equal T
        layer = 2
        dist = 300

        for pol in ["s", "p"]:
            coh_tmm_data = coh_tmm(
                pol, self.n_list, self.d_list, self.th_0, self.lam_vac
            )
            data = position_resolved(layer, dist, coh_tmm_data)
            np.testing.assert_allclose(data["poyn"], coh_tmm_data["T"])

    def test_poynting_vector_start(self):
        # Poynting vector at start should equal power_entering
        layer = 1
        dist = 0

        for pol in ["s", "p"]:
            coh_tmm_data = coh_tmm(
                pol, self.n_list, self.d_list, self.th_0, self.lam_vac
            )
            data = position_resolved(layer, dist, coh_tmm_data)
            np.testing.assert_allclose(data["poyn"], coh_tmm_data["power_entering"])

    def test_poynting_vector_continuity(self):
        # Poynting vector should be continuous

        for pol in ["s", "p"]:
            layer = 1
            dist = 100
            coh_tmm_data = coh_tmm(
                pol, self.n_list, self.d_list, self.th_0, self.lam_vac
            )
            data = position_resolved(layer, dist, coh_tmm_data)
            poyn1 = data["poyn"]
            layer = 2
            dist = 0
            coh_tmm_data = coh_tmm(
                pol, self.n_list, self.d_list, self.th_0, self.lam_vac
            )
            data = position_resolved(layer, dist, coh_tmm_data)
            poyn2 = data["poyn"]
            np.testing.assert_allclose(poyn1, poyn2)


class TestAbsorbAnalytics(unittest.TestCase):
    """
    Test absorp_analytic_fn functions
    """

    d_list = [np.inf, 100, 300, np.inf]  # in nm
    n_list = [1, 2.2 + 0.2j, 3.3 + 0.3j, 1]
    th_0 = np.pi / 4
    lam_vac = 400
    layer = 1
    d = d_list[layer]
    dist = 37

    def test1(self):
        layer = 1
        d = self.d_list[layer]
        dist = 37
        for pol in ["s", "p"]:
            coh_tmm_data = coh_tmm(
                pol, self.n_list, self.d_list, self.th_0, self.lam_vac
            )
            expected_absorp = position_resolved(layer, dist, coh_tmm_data)["absor"]
            absorp_fn = absorp_analytic_fn()
            absorp_fn.fill_in(coh_tmm_data, layer)
            np.testing.assert_allclose(absorp_fn.run(dist), expected_absorp)
            absorp_fn2 = absorp_fn.copy().flip()
            dist_from_other_side = d - dist
            np.testing.assert_allclose(
                absorp_fn2.run(dist_from_other_side), expected_absorp
            )


class TestIncoherentStack(unittest.TestCase):
    def test_incoherent1(self):
        """
        test inc_tmm(). To do: Add more tests.
        """

        # 3-incoherent-layer test, real refractive indices (so that R and T are the
        # same in both directions)
        n0 = 1
        n1 = 2
        n2 = 3
        n_list = [n0, n1, n2]
        d_list = [np.inf, 567, np.inf]
        c_list = ["i", "i", "i"]
        th0 = np.pi / 3
        th1 = snell(n0, n1, th0)
        th2 = snell(n0, n2, th0)
        lam_vac = 400

        for pol in ["s", "p"]:
            inc_data = inc_tmm(pol, n_list, d_list, c_list, th0, lam_vac)
            R0 = abs(interface_r(pol, n0, n1, th0, th1) ** 2)
            R1 = abs(interface_r(pol, n1, n2, th1, th2) ** 2)
            T0 = 1 - R0
            RR = R0 + R1 * T0**2 / (1 - R0 * R1)
            np.testing.assert_allclose(inc_data["R"], RR)
            np.testing.assert_allclose(inc_data["R"] + inc_data["T"], 1)

    def test_incoherent2(self):
        # One finite layer with incoherent layers on both sides. Should agree with
        # coherent program
        n0 = 1 + 0.1j
        n1 = 2 + 0.2j
        n2 = 3 + 0.4j
        n_list = [n0, n1, n2]
        d_list = [np.inf, 100, np.inf]
        c_list = ["i", "c", "i"]
        n00 = 1
        th00 = np.pi / 3
        th0 = snell(n00, n0, th00)
        lam_vac = 400
        for pol in ["s", "p"]:
            inc_data = inc_tmm(pol, n_list, d_list, c_list, th0, lam_vac)
            coh_data = coh_tmm(pol, n_list, d_list, th0, lam_vac)
            np.testing.assert_allclose(inc_data["R"], coh_data["R"])
            np.testing.assert_allclose(inc_data["T"], coh_data["T"])
            np.testing.assert_allclose(1, sum(inc_absorp_in_each_layer(inc_data)))

    def test_incoherent3(self):
        # One finite layer with three incoherent layers. Should agree with
        # manual calculation + coherent program
        n0 = 1 + 0.1j
        n1 = 2 + 0.2j
        n2 = 3 + 0.004j
        n3 = 4 + 0.2j
        d1 = 100
        d2 = 10000
        n_list = [n0, n1, n2, n3]
        d_list = [np.inf, d1, d2, np.inf]
        c_list = ["i", "c", "i", "i"]
        n00 = 1
        th00 = np.pi / 3
        th0 = snell(n00, n0, th00)
        lam_vac = 400
        for pol in ["s", "p"]:
            inc_data = inc_tmm(pol, n_list, d_list, c_list, th0, lam_vac)
            coh_data = coh_tmm(pol, [n0, n1, n2], [np.inf, d1, np.inf], th0, lam_vac)
            th2 = snell(n0, n2, th0)
            th3 = snell(n0, n3, th0)
            coh_bdata = coh_tmm(pol, [n2, n1, n0], [np.inf, d1, np.inf], th2, lam_vac)
            R02 = coh_data["R"]
            R20 = coh_bdata["R"]
            T02 = coh_data["T"]
            T20 = coh_bdata["T"]
            P2 = np.exp(
                -4 * np.pi * d2 * (n2 * np.cos(th2)).imag / lam_vac
            )  # fraction passing through
            R23 = interface_R(pol, n2, n3, th2, th3)
            T23 = interface_T(pol, n2, n3, th2, th3)
            # T = T02 * P2 * T23 + T02 * P2 * R23 * P2 * R20 * P2 * T23 + ...
            T = T02 * P2 * T23 / (1 - R23 * P2 * R20 * P2)
            # R = R02
            #    + T02 * P2 * R23 * P2 * T20
            #    + T02 * P2 * R23 * P2 * R20 * P2 * R23 * P2 * T20 + ...
            R = R02 + T02 * P2 * R23 * P2 * T20 / (1 - R20 * P2 * R23 * P2)
            np.testing.assert_allclose(inc_data["T"], T)
            np.testing.assert_allclose(inc_data["R"], R)

    def test_incoherent4(self):
        # The coherent program with a thick but randomly-varying-thickness substrate
        # should agree with the incoherent program.
        nair = 1 + 0.1j
        nfilm = 2 + 0.2j
        nsub = 3
        nf = 3 + 0.4j
        n_list = [nair, nfilm, nsub, nf]
        n00 = 1
        th00 = np.pi / 3
        th0 = snell(n00, nair, th00)
        lam_vac = 400
        for pol in ["s", "p"]:
            d_list_inc = [np.inf, 100, 1, np.inf]  # sub thickness doesn't matter here
            c_list = ["i", "c", "i", "i"]
            inc_data = inc_tmm(pol, n_list, d_list_inc, c_list, th0, lam_vac)
            coh_Rs = []
            coh_Ts = []
            for dsub in np.linspace(10000, 30000, 357):
                d_list = [np.inf, 100, dsub, np.inf]
                coh_data = coh_tmm(pol, n_list, d_list, th0, lam_vac)
                coh_Rs.append(coh_data["R"])
                coh_Ts.append(coh_data["T"])
            np.testing.assert_allclose(np.average(coh_Rs), inc_data["R"])
            np.testing.assert_allclose(np.average(coh_Ts), inc_data["T"])

    def test_incoherent5(self):
        # The coherent program with a thick substrate and randomly-varying wavelength
        # should agree with the incoherent program.
        n0 = 1 + 0.0j
        n_list = [n0, 2 + 0.0002j, 3 + 0.0001j, 3 + 0.4j]
        n00 = 1
        th00 = np.pi / 3
        th0 = snell(n00, n0, th00)
        d_list = [np.inf, 10000, 10200, np.inf]
        c_list = ["i", "i", "i", "i"]
        for pol in ["s", "p"]:
            inc_absorp = np.array([0.0, 0.0, 0.0, 0.0])
            coh_absorp = np.array([0.0, 0.0, 0.0, 0.0])
            num_pts = 234
            for lam_vac in np.linspace(40, 50, num_pts):
                inc_data = inc_tmm(pol, n_list, d_list, c_list, th0, lam_vac)
                inc_absorp += np.array(inc_absorp_in_each_layer(inc_data))
                coh_data = coh_tmm(pol, n_list, d_list, th0, lam_vac)
                coh_absorp += np.array(absorp_in_each_layer(coh_data))
            inc_absorp /= num_pts
            coh_absorp /= num_pts
            print(
                "Coherent with random wavelength should agree with incoherent. "
                + "The two rows of this array should be the same:"
            )
            np.testing.assert_allclose(inc_absorp, coh_absorp)


class TestRT(unittest.TestCase):
    def test_RT1(self):
        """
        Tests of formulas for R and T
        """
        # When ni is real [see manual], R+T should equal 1
        ni = 2
        nf = 3.0 + 0.2j
        thi = np.pi / 5
        thf = snell(ni, nf, thi)
        for pol in ["s", "p"]:
            T = interface_T(pol, ni, nf, thi, thf)
            R = interface_R(pol, ni, nf, thi, thf)
            np.testing.assert_allclose(1, R + T)

    def test_RT2(self):
        # For a single interface, power_entering should equal T
        ni = 2 + 0.1j
        n00 = 1
        th00 = np.pi / 5
        thi = snell(n00, ni, th00)
        nf = 3.0 + 0.2j
        thf = snell(ni, nf, thi)
        for pol in ["s", "p"]:
            r = interface_r(pol, ni, nf, thi, thf)
            pe = power_entering_from_r(pol, r, ni, thi)
            T = interface_T(pol, ni, nf, thi, thf)
            np.testing.assert_allclose(pe, T)


class TestOverflow(unittest.TestCase):
    def test_coh_overflow(self):
        """
        Test whether very very opaque layers will break the coherent program
        """
        n_list = [1.0, 2 + 0.1j, 1 + 3j, 4.0, 5.0]
        d_list = [np.inf, 50, 1e5, 50, np.inf]
        lam = 200
        # alpha_d = np.imag(n_list[2]) * 4 * np.pi * d_list[2] / lam
        # Very opaque layer: Calculation should involve e^(-alpha_d)!)
        data = coh_tmm("s", n_list, d_list, 0, lam)
        n_list2 = n_list[0:3]
        d_list2 = d_list[0:3]
        d_list2[-1] = np.inf
        data2 = coh_tmm("s", n_list2, d_list2, 0, lam)

        np.testing.assert_allclose(
            data["vw_list"][:3], data2["vw_list"][:3], atol=1e-16
        )

    def test_inc_overflow(self):
        """
        Test whether very very opaque layers will break the incoherent program
        """
        n_list = [1.0, 2.0, 1 + 3j, 4.0, 5.0]
        d_list = [np.inf, 50, 1e5, 50, np.inf]
        c_list = ["i", "i", "i", "i", "i"]
        lam = 200
        # alpha_d = np.imag(n_list[2]) * 4 * np.pi * d_list[2] / lam
        # "Very opaque layer: Calculation should involve e^(-alpha_d)!"
        data = inc_tmm("s", n_list, d_list, c_list, 0, lam)
        n_list2 = n_list[0:3]
        d_list2 = d_list[0:3]
        d_list2[-1] = np.inf
        c_list2 = c_list[0:3]
        data2 = inc_tmm("s", n_list2, d_list2, c_list2, 0, lam)

        np.testing.assert_allclose(
            data["power_entering_list"][0:3], data2["power_entering_list"][0:3]
        )


class TestVectorized(unittest.TestCase):
    """Run some of the previous tests with vectorized data"""

    n_list = [1, 2 + 4j, 3 + 0.3j, 1 + 0.1j]
    n_list = np.stack([n_list] * 2, axis=-1)
    d_list = [np.inf, 2, 3, np.inf]
    th_0 = 0.1
    lam_vac = [100, 100]

    def test_coh_stack_rt(self):
        s_data = coh_tmm("s", self.n_list, self.d_list, self.th_0, self.lam_vac)

        np.testing.assert_allclose(
            s_data["r"], [-0.60331226568845775 - 0.093522181653632019j] * 2
        )
        np.testing.assert_allclose(
            s_data["t"], [0.44429533471192989 + 0.16921936169383078j] * 2
        )
        np.testing.assert_allclose(s_data["R"], [0.37273208839139516] * 2)
        np.testing.assert_allclose(s_data["T"], [0.22604491247079261] * 2)

        p_data = coh_tmm("p", self.n_list, self.d_list, self.th_0, self.lam_vac)

        np.testing.assert_allclose(
            p_data["r"], [0.60102654255772481 + 0.094489146845323682j] * 2
        )
        np.testing.assert_allclose(
            p_data["t"], [0.4461816467503148 + 0.17061408427088917j] * 2
        )
        np.testing.assert_allclose(p_data["R"], [0.37016110373044969] * 2)
        np.testing.assert_allclose(p_data["T"], [0.22824374314132009] * 2)

    def test_coh_stack_ellipso(self):
        ellips_data = ellips(self.n_list, self.d_list, self.th_0, self.lam_vac)
        np.testing.assert_allclose(ellips_data["psi"], [0.78366777347038352] * 2)
        np.testing.assert_allclose(ellips_data["Delta"], [0.0021460774404193292] * 2)

    def test_incoherent(self):
        # 3-incoherent-layer test, real refractive indices (so that R and T are the
        # same in both directions)
        n0 = 1
        n1 = 2
        n2 = 3
        n_list = [n0, n1, n2]
        n_list = np.stack([n_list] * 2, axis=-1)
        d_list = [np.inf, 567, np.inf]
        c_list = ["i", "i", "i"]
        th0 = np.pi / 3
        th1 = snell(n0, n1, th0)
        th2 = snell(n0, n2, th0)
        lam_vac = [400, 400]

        for pol in ["s", "p"]:
            inc_data = inc_tmm(pol, n_list, d_list, c_list, th0, lam_vac)
            R0 = abs(interface_r(pol, n0, n1, th0, th1) ** 2)
            R1 = abs(interface_r(pol, n1, n2, th1, th2) ** 2)
            T0 = 1 - R0
            RR = R0 + R1 * T0**2 / (1 - R0 * R1)
            np.testing.assert_allclose(inc_data["R"], RR)
            np.testing.assert_allclose(inc_data["R"] + inc_data["T"], [1, 1])

    def test_poynting(self):
        d_list = [np.inf, 100, 300, np.inf]  # in nm
        n00 = 1
        th00 = np.pi / 4
        n0 = 1 + 0.1j
        th_0 = snell(n00, n0, th00)
        n_list = [n0, 2.2 + 0.2j, 3.3 + 0.3j, 1 + 0.4j]
        n_list = np.stack([n_list] * 2, axis=-1)
        lam_vac = [400, 400]
        layer = 1
        dist = 37

        for pol in ["s", "p"]:
            coh_tmm_data = coh_tmm(pol, n_list, d_list, th_0, lam_vac)
            data = position_resolved(layer, dist, coh_tmm_data)
            np.testing.assert_allclose([1, 1], sum(absorp_in_each_layer(coh_tmm_data)))

    def test_absorb_analytics(self):
        d_list = [np.inf, 100, 300, np.inf]  # in nm
        n_list = [1, 2.2 + 0.2j, 3.3 + 0.3j, 1]
        n_list = np.stack([n_list] * 2, axis=-1)
        th_0 = np.pi / 4
        lam_vac = [400, 400]
        layer = 1
        d = d_list[layer]
        dist = 37

        layer = 1
        d = d_list[layer]
        dist = 37
        for pol in ["s", "p"]:
            coh_tmm_data = coh_tmm(pol, n_list, d_list, th_0, lam_vac)
            expected_absorp = position_resolved(layer, dist, coh_tmm_data)["absor"]
            absorp_fn = absorp_analytic_fn()
            absorp_fn.fill_in(coh_tmm_data, layer)
            np.testing.assert_allclose(absorp_fn.run(dist), expected_absorp)
            absorp_fn2 = absorp_fn.copy().flip()
            dist_from_other_side = d - dist
            np.testing.assert_allclose(
                absorp_fn2.run(dist_from_other_side), expected_absorp
            )
