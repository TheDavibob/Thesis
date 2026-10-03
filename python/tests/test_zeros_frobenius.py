"""Root finding (``zeros``) and series solutions about regular singular points
(``frobenius``)."""
import numpy as np

from shear_flow_scattering import frobenius, zeros


class TestNewtonAndArgumentPrinciple:
    f = staticmethod(lambda z: (z - 0.3) * (z + 0.2j))
    df = staticmethod(lambda z: 2 * z - 0.3 + 0.2j)

    def test_newton_raphson(self):
        assert np.isclose(zeros.newton_raphson(self.f, self.df, 0.5 + 0.1j), 0.3)

    def test_zero_count_and_simple_finder(self):
        assert np.isclose(zeros.zeros_poles(self.f, self.df, [-1, 1], [-1, 1]), 2)
        zz = np.sort_complex(zeros.simple_zero_finder(self.f, self.df, [-1, 1], [-1, 1]))
        assert np.allclose(zz, np.sort_complex([0.3, -0.2j]))

    def test_complete_finder_many_zeros(self):
        # box edges chosen so no zero sits on a bisection line
        found = np.sort(zeros.complete_zero_finder(lambda z: np.sin(3 * z), lambda z: 3 * np.cos(3 * z),
                                                   [-3.3, 3.2], [-1, 0.9]).real)
        assert np.allclose(found, np.pi / 3 * np.arange(-3, 4), atol=1e-6)


class TestRecursiveZeroFinder:
    def test_finds_lower_half_plane_zeros(self):
        roots = np.array([1.0, 2 + 1j, 1 - 1j])
        f = lambda z: np.prod([z - r for r in roots], axis=0)
        out = zeros.zero_finder_recursive(f, [1.1, 1.9 + 0.9j, 1.1 - 0.9j])
        key = lambda a: sorted(a, key=lambda z: (round(z.real, 6), round(z.imag, 6)))
        assert np.allclose(key(out.z), key(roots), atol=1e-8)


class TestZerosBox:
    def test_counts_by_half_plane(self):
        f = lambda z: (z - 0.5) * (z - 0.2j) * (z + 0.3 + 0.4j)
        Z = zeros.zeros_box(f, [-1, 1], [-1, 1])
        assert (Z.z, Z.r, Z.U, Z.L) == (3, 1, 1, 1)


class TestFrobeniusTwoSolutions:
    @staticmethod
    def _residual(fn, dfn, z, P, Q):
        h = 1e-4
        d2 = (dfn(z + h) - dfn(z - h)) / (2 * h)
        return abs(z ** 2 * d2 + P(z) * z * dfn(z) + Q(z) * fn(z))

    def test_bessel_j0_repeated_root(self):
        y = frobenius.two_solutions_complex(lambda z: 1 + 0 * z, lambda z: z ** 2)
        assert y.I == 2
        assert np.isclose(y.y1(0.5 + 0j), 0.9384698072408129)
        for fn, dfn in ((y.y1, y.dy1), (y.y2, y.dy2)):
            assert self._residual(fn, dfn, 0.5 + 0.1j, lambda z: 1, lambda z: z ** 2) < 1e-6

    def test_roots_differing_by_an_integer(self):
        # z^2 y'' - 2 y + z^2 y = 0 : indicial roots 2 and -1
        y = frobenius.two_solutions_complex(lambda z: 0 * z, lambda z: -2 + z ** 2)
        assert y.I == 3 and np.isclose(y.s1, 2) and np.isclose(y.s2, -1)
        for fn, dfn in ((y.y1, y.dy1), (y.y2, y.dy2)):
            assert self._residual(fn, dfn, 0.4 + 0.2j, lambda z: 0, lambda z: -2 + z ** 2) < 1e-5
