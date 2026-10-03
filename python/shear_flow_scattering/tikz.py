"""Export matplotlib axes to pgfplots/TikZ, in the style of matlab2tikz.

The thesis and paper figures were made by plotting in MATLAB and then running
matlab2tikz by hand (no ``.m`` file calls it; only its output survives, e.g.
``code/GenericInversion/Scattering/WallPressure_Compact.tex``).  This module
reproduces that output for the Python port so that the LaTeX side is unchanged:

* one ``tikzpicture`` per axes, sized by ``\\fwidth`` / ``\\fheight`` (define
  them before ``\\input``, as for matlab2tikz);
* line data written inline as ``\\addplot table[row sep=crcr]``, with NaNs
  breaking the line (``unbounded coords=jump``);
* field plots (``pcolormesh`` / ``imshow``) rasterised to ``<stem>-1.png`` and
  placed with ``\\addplot graphics`` inside the axis, with the colormap and a
  pgfplots ``colorbar`` (so ticks and labels stay vector and use the
  document's fonts);
* polar axes drawn on a Cartesian axis with a hand-drawn grid, as matlab2tikz
  does for ``polarplot``.

Only what the example scripts use is supported: lines (with markers),
``pcolormesh`` on a rectilinear grid, ``imshow``, colorbars, log axes, titles,
labels and legends.  Text annotations, patches and twin axes are ignored.

Usage::

    from shear_flow_scattering.tikz import save_tikz
    save_tikz(ax, "WallPressure_base.tex")

LaTeX side::

    \\usepackage{pgfplots}\\pgfplotsset{compat=1.16}
    \\newlength\\fheight \\newlength\\fwidth
    \\setlength\\fwidth{0.8\\linewidth}\\setlength\\fheight{0.5\\fwidth}
    \\input{WallPressure_base.tex}
"""
import os

import numpy as np

_DASH = {"-": None, "--": "dashed", ":": "dotted", "-.": "dashdotted"}
_MARK = {"o": "*", ".": "*", "s": "square*", "^": "triangle*", "v": "triangle*, mark options={rotate=180}",
         "d": "diamond*", "D": "diamond*", "x": "x", "+": "+", "*": "asterisk"}
_GRID = "white!83!darkgray"


def _num(v):
    return "nan" if not np.isfinite(v) else f"{v:.6g}"


def _table(x, y):
    rows = "\n".join(f"{_num(a)}\t{_num(b)}\\\\" for a, b in zip(x, y))
    return "table[row sep=crcr]{%\n" + rows + "\n};\n"


class _Colours:
    """Collect colours as ``mycolorN`` definitions, as matlab2tikz does."""

    def __init__(self):
        self.names = {}

    def __call__(self, rgba):
        from matplotlib.colors import to_rgb
        rgb = tuple(round(float(c), 5) for c in to_rgb(rgba))
        if rgb == (0.0, 0.0, 0.0):
            return "black"
        if rgb == (1.0, 1.0, 1.0):
            return "white"
        if rgb not in self.names:
            self.names[rgb] = f"mycolor{len(self.names) + 1}"
        return self.names[rgb]

    def preamble(self):
        return "".join(f"\\definecolor{{{n}}}{{rgb}}{{{r:.5f},{g:.5f},{b:.5f}}}%\n"
                       for (r, g, b), n in self.names.items())


def _line_options(line, colours, lw_scale):
    opts = [f"color={colours(line.get_color())}"]
    ls = line.get_linestyle()
    if ls in ("None", "none", " ", ""):
        opts.append("draw=none")
    elif _DASH.get(ls):
        opts.append(_DASH[ls])
    opts.append(f"line width={line.get_linewidth() * lw_scale:.2g}pt")
    m = line.get_marker()
    if m not in (None, "None", "none", "", " "):
        opts.append(f"mark={_MARK.get(m, '*')}")
        opts.append(f"mark size={line.get_markersize() / 2.5:.2g}pt")
        opts.append(f"mark options={{solid, fill={colours(line.get_markerfacecolor())}}}")
    else:
        opts.append("mark=none")
    return opts


def _legend_label(line, ax):
    label = line.get_label()
    shown = ax.get_legend() is not None and not label.startswith("_")
    return label if shown else None


def _colormap(cmap, n=64):
    cols = cmap(np.linspace(0, 1, n))
    stops = "; ".join(f"rgb({i}pt)=({c[0]:.6g},{c[1]:.6g},{c[2]:.6g})" for i, c in enumerate(cols))
    return f"colormap={{mymap}}{{[1pt] {stops}}}"


def _raster(mappable, png_path):
    """Write the mappable's colours to a PNG oriented for ``\\addplot graphics``;
    return its extent (xmin, xmax, ymin, ymax)."""
    import matplotlib.image as mpimg
    from matplotlib.collections import QuadMesh

    if isinstance(mappable, QuadMesh):
        xy = mappable.get_coordinates()
        X, Y = xy[..., 0], xy[..., 1]
        if not (np.allclose(X, X[:1, :]) and np.allclose(Y, Y[:, :1])):
            raise ValueError("save_tikz only rasterises pcolormesh on a rectilinear grid")
        rgba = mappable.to_rgba(np.ma.asarray(mappable.get_array()).reshape(X.shape[0] - 1, X.shape[1] - 1))
        x, y = X[0], Y[:, 0]
        if x[-1] < x[0]:
            rgba = rgba[:, ::-1]
        if y[-1] > y[0]:            # image row 0 is drawn at the top (ymax)
            rgba = rgba[::-1]
        extent = (x.min(), x.max(), y.min(), y.max())
    else:                           # AxesImage
        rgba = mappable.to_rgba(mappable.get_array())
        l, r, b, t = mappable.get_extent()
        if r < l:
            rgba = rgba[:, ::-1]
        if mappable.origin == "lower":
            rgba = rgba[::-1]
        extent = (min(l, r), max(l, r), min(b, t), max(b, t))
    mpimg.imsave(png_path, rgba)
    return extent


def _polar(ax, colours, lw_scale, legend):
    """Body of a polar axes drawn on a Cartesian axis."""
    rmax = ax.get_rmax()
    offset, direction = ax.get_theta_offset(), ax.get_theta_direction()
    to_xy = lambda th, r: (r * np.cos(direction * np.asarray(th) + offset),
                           r * np.sin(direction * np.asarray(th) + offset))
    body = []
    t = np.linspace(0, 2 * np.pi, 181)
    body.append(f"\\addplot[area legend, draw={_GRID}, fill=white, forget plot]\n" + _table(*to_xy(t, rmax)))
    for r in [r for r in ax.get_yticks() if 0 < r < rmax * (1 - 1e-9)]:
        body.append(f"\\addplot[color={_GRID}, line width=0.5pt, forget plot]\n" + _table(*to_xy(t, r)))
    for a in np.deg2rad(np.arange(0, 360, 45)):
        body.append(f"\\addplot[color={_GRID}, line width=0.5pt, forget plot]\n" + _table(*to_xy([a, a], [0, rmax])))
        lx, ly = to_xy(a, 1.12 * rmax)
        body.append(f"\\node at (axis cs:{lx:.4g},{ly:.4g}) {{${np.rad2deg(a):.0f}^\\circ$}};\n")
    lab = np.deg2rad(22.5)
    for r in [r for r in ax.get_yticks() if 0 < r <= rmax * (1 + 1e-9)]:
        lx, ly = to_xy(lab, r)
        body.append(f"\\node[anchor=south west, font=\\scriptsize] at (axis cs:{lx:.4g},{ly:.4g}) {{{r:g}}};\n")
    for line in ax.get_lines():
        opts = _line_options(line, colours, lw_scale)
        label = _legend_label(line, ax) if legend else None
        if label is None:
            opts.append("forget plot")
        th, r = line.get_data()
        body.append(f"\\addplot[{', '.join(opts)}]\n" + _table(*to_xy(th, r)))
        if label is not None:
            body.append(f"\\addlegendentry{{{label}}}\n")
    lim = 1.25 * rmax
    axis = [f"xmin={-lim:.4g}", f"xmax={lim:.4g}", f"ymin={-lim:.4g}", f"ymax={lim:.4g}",
            "axis equal image", "axis line style={draw=none}", "ticks=none"]
    return axis, body


def save_tikz(ax, path, width=r"0.951\fwidth", height=r"\fheight", graphics_dir=None,
              lw_scale=1.0, legend=True, title=True):
    """Write axes ``ax`` to ``path`` as a matlab2tikz-style ``tikzpicture``.

    ``width``/``height`` are pgfplots lengths (``height`` is ignored for polar
    axes, which keep a 1:1 aspect).  Rasters are written next to ``path`` as
    ``<stem>-1.png``; ``graphics_dir`` is the directory prefix LaTeX should use
    to find them (default: the bare file name, i.e. relative to the main
    document, as matlab2tikz does).  ``lw_scale`` scales matplotlib line widths
    (points) to pgfplots points.
    """
    colours = _Colours()
    stem = os.path.splitext(path)[0]
    opts = [f"width={width}"]
    body = []
    if ax.name == "polar":
        axis, body = _polar(ax, colours, lw_scale, legend)
        opts += axis
    else:
        opts.append(f"height={height}")
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        opts += [f"xmin={_num(min(x0, x1))}", f"xmax={_num(max(x0, x1))}",
                 f"ymin={_num(min(y0, y1))}", f"ymax={_num(max(y0, y1))}"]
        if x1 < x0:
            opts.append("x dir=reverse")
        if y1 < y0:
            opts.append("y dir=reverse")
        for a in "xy":
            if getattr(ax, f"get_{a}scale")() == "log":
                opts += [f"{a}mode=log", f"{a}minorticks=true"]
            label = getattr(ax, f"get_{a}label")()
            if label:
                opts.append(f"{a}label={{{label}}}")
        mappables = list(ax.collections) + list(ax.images)
        mappables = [m for m in mappables if hasattr(m, "get_array") and m.get_array() is not None]
        for i, m in enumerate(mappables, 1):
            png = f"{stem}-{i}.png"
            xa, xb, ya, yb = _raster(m, png)
            name = os.path.basename(png) if graphics_dir is None else graphics_dir.rstrip("/") + "/" + os.path.basename(png)
            body.append(f"\\addplot[forget plot] graphics[xmin={_num(xa)}, xmax={_num(xb)}, "
                        f"ymin={_num(ya)}, ymax={_num(yb)}] {{{name}}};\n")
            if i == 1:
                opts += ["axis on top", _colormap(m.cmap),
                         f"point meta min={_num(m.norm.vmin)}", f"point meta max={_num(m.norm.vmax)}"]
                if getattr(m, "colorbar", None) is not None:
                    opts.append("colorbar")
                    clabel = m.colorbar.ax.get_ylabel()
                    if clabel:
                        opts.append(f"colorbar style={{ylabel={{{clabel}}}}}")
        jump = False
        for line in ax.get_lines():
            lo = _line_options(line, colours, lw_scale)
            label = _legend_label(line, ax) if legend else None
            if label is None:
                lo.append("forget plot")
            x, y = (np.asarray(v, dtype=float) for v in line.get_data())
            jump |= not np.all(np.isfinite(y))
            body.append(f"\\addplot[{', '.join(lo)}]\n" + _table(x, y))
            if label is not None:
                body.append(f"\\addlegendentry{{{label}}}\n")
        if jump:
            opts.append("unbounded coords=jump")
    opts.append("scale only axis")
    if title and ax.get_title():
        opts.append(f"title={{{ax.get_title()}}}")
    if legend and ax.get_legend() is not None:
        opts.append("legend style={legend cell align=left, align=left, draw=white!15!black, font=\\footnotesize}")
    text = ("% This file was created by shear_flow_scattering.tikz (matlab2tikz-compatible).\n%\n"
            + colours.preamble() + "%\n\\begin{tikzpicture}\n\n\\begin{axis}[%\n"
            + ",\n".join(opts) + "\n]\n" + "".join(body) + "\\end{axis}\n\\end{tikzpicture}%\n")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)
    return path


def save_figure_tikz(fig, stem, names=None, **kw):
    """Write every axes of ``fig`` (colorbar axes excluded) to ``<stem>_<name>.tex``.

    ``names`` defaults to ``a, b, c, ...`` in axes order.
    """
    axes = [a for a in fig.axes if not hasattr(a, "_colorbar") and a.get_label() != "<colorbar>"]
    names = names or [chr(ord("a") + i) for i in range(len(axes))]
    return [save_tikz(a, f"{stem}_{n}.tex", **kw) for a, n in zip(axes, names)]
