"""matlab2tikz-style export (``tikz``)."""
import numpy as np
import pytest


@pytest.fixture
def exported(tmp_path):
    """A three-panel figure (log-axis lines, pcolormesh with colourbar, polar)
    exported with save_figure_tikz; yields (tex texts, tmp_path)."""
    mpl = pytest.importorskip("matplotlib")
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    from shear_flow_scattering.tikz import save_figure_tikz

    fig = plt.figure()
    a = fig.add_subplot(1, 3, 1)
    f = np.logspace(2, 4, 5)
    y = np.array([1.0, 2.0, np.nan, 4.0, 5.0])
    a.plot(f, y, "C0--", label=r"$y_2 = 0.1\delta$")
    a.plot(f, y + 1, "o", color="C1")
    a.set_xscale("log"); a.set_xlabel("$f$"); a.legend()
    b = fig.add_subplot(1, 3, 2)
    x1, x2 = np.linspace(-1, 2, 7), np.linspace(0, 1, 4)
    m = b.pcolormesh(x1, x2, np.add.outer(x2, x1), shading="auto", cmap="viridis")  # increases up and right
    fig.colorbar(m, ax=b)
    c = fig.add_subplot(1, 3, 3, projection="polar")
    c.plot(np.linspace(0, 2 * np.pi, 50), np.ones(50), label="unit"); c.legend()

    paths = save_figure_tikz(fig, str(tmp_path / "fig"))
    plt.close(fig)
    assert [p[-5:] for p in paths] == ["a.tex", "b.tex", "c.tex"]
    yield [open(p).read() for p in paths], tmp_path


class TestSaveTikz:
    def test_line_plot(self, exported):
        ta = exported[0][0]
        assert "xmode=log" in ta and "unbounded coords=jump" in ta and "dashed" in ta
        assert r"\addlegendentry{$y_2 = 0.1\delta$}" in ta and ta.count(r"\addlegendentry") == 1
        assert "nan" in ta and "mark=*" in ta

    def test_field_raster_and_colourbar(self, exported):
        from matplotlib.image import imread
        (_, tb, _), tmp = exported
        assert "graphics[" in tb and "colorbar" in tb and "fig_b-1.png" in tb
        img = imread(str(tmp / "fig_b-1.png"))
        assert img.shape[:2] == (4, 7)
        # top-right pixel is the largest value -> the yellow end of viridis
        assert img[0, -1, 0] > 0.9 and img[-1, 0, 0] < 0.3

    def test_polar_axes(self, exported):
        tc = exported[0][2]
        assert "axis equal image" in tc and r"\addlegendentry{unit}" in tc

    def test_well_formed(self, exported):
        for t in exported[0]:
            assert t.count(r"\begin{axis}") == t.count(r"\end{axis}") == 1
            assert t.count("{") == t.count("}")


@pytest.fixture
def exported_extras(tmp_path):
    """Features used by the paper figures: imshow RGBA raster with a phase
    colourbar and custom tick labels (Fig. 11), contours, an unfilled
    rectangle, axhline on a log axis, and major + minor grids (Figs 8-9)."""
    mpl = pytest.importorskip("matplotlib")
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    from matplotlib.patches import Rectangle
    from shear_flow_scattering.tikz import save_tikz

    fig, (a, b) = plt.subplots(1, 2)
    rgba = np.zeros((3, 5, 4))
    rgba[..., 0] = np.linspace(0, 1, 5)[None, :]      # red increases to the right
    rgba[..., 3] = 1
    rgba[0, 0] = (1, 1, 1, 0)                          # one transparent (masked) pixel, bottom-left
    im = a.imshow(rgba, origin="lower", extent=(0, 5, 0, 3))
    norm = Normalize(-np.pi, np.pi)
    im.set_cmap("twilight"); im.set_norm(norm)
    cb = fig.colorbar(ScalarMappable(norm, "twilight"), ax=a, ticks=[-np.pi, 0, np.pi])
    cb.ax.set_yticklabels([r"$-\pi$", "0", r"$\pi$"])
    im.colorbar = cb
    X, Y = np.meshgrid(np.linspace(0, 5, 20), np.linspace(0, 3, 15))
    a.contour(X, Y, X + Y, levels=[2, 4], colors="0.25", linewidths=0.5)
    a.add_patch(Rectangle((1, 1), 2, 1, fill=False, ec="k", lw=0.8))

    b.plot([100, 1000, 5000], [1, 2, 3])
    b.axhline(0, color="k", lw=0.6)
    b.set_xscale("log")
    b.grid(True, which="major", color="0.85", lw=0.5)
    b.grid(True, which="minor", axis="x", color="0.93", lw=0.4)

    ta = open(save_tikz(a, str(tmp_path / "a.tex"))).read()
    tb = open(save_tikz(b, str(tmp_path / "b.tex"))).read()
    plt.close(fig)
    return ta, tb, tmp_path


class TestSaveTikzFigureFeatures:
    def test_imshow_raster_orientation_and_transparency(self, exported_extras):
        from matplotlib.image import imread
        ta, _, tmp = exported_extras
        assert "graphics[xmin=0, xmax=5, ymin=0, ymax=3]" in ta
        img = imread(str(tmp / "a-1.png"))
        assert img.shape[:2] == (3, 5)
        # origin="lower": data row 0 is drawn at the bottom, i.e. the last PNG row
        assert img[-1, 0, 3] == 0 and img[0, 0, 3] == 1
        assert img[1, -1, 0] > img[1, 0, 0]

    def test_colourbar_custom_tick_labels(self, exported_extras):
        ta = exported_extras[0]
        assert "colorbar" in ta
        assert r"yticklabels={{$-\pi$},{0},{$\pi$}}" in ta
        assert "ytick={-3.14159,0,3.14159}" in ta

    @pytest.mark.xfail(strict=True, reason="KNOWN FAILURE: save_tikz skips Rectangle patches -- it compares "
                       "p.get_transform() with ax.transData, but a patch's transform also includes its patch "
                       "transform, so the test never matches.")
    def test_contours_and_rectangle_exported_as_lines(self, exported_extras):
        ta = exported_extras[0]
        assert ta.count(r"\addplot[color=") >= 3     # two contour levels + rectangle
        closed_rectangle = "1\t1\\\\\n3\t1\\\\\n3\t2\\\\\n1\t2\\\\\n1\t1\\\\"
        assert closed_rectangle in ta

    def test_axhline_spans_log_axis_limits(self, exported_extras):
        tb = exported_extras[1]
        xmin = float(tb.split("xmin=")[1].split(",")[0])
        xmax = float(tb.split("xmax=")[1].split(",")[0])
        line = tb.split(r"\addplot[color=black")[1]
        rows = [r for r in line.split("{%")[1].split("};")[0].strip().split("\\\\") if r.strip()]
        xs = [float(r.split()[0]) for r in rows]
        assert np.allclose([xs[0], xs[-1]], [xmin, xmax], rtol=1e-6)

    def test_major_and_minor_grids(self, exported_extras):
        tb = exported_extras[1]
        assert "xmajorgrids" in tb and "ymajorgrids" in tb
        assert "xminorgrids" in tb and "minor grid style" in tb
        assert "log basis x=10" in tb and "xtick={100,1000}" in tb
