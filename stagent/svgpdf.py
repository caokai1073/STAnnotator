"""The page's figures (dot plots, bubble plots, circle plots: SVG it draws itself) as PDF, drawn again with matplotlib.

Only what those figures use: svg, g, text (text-anchor, font-size, font-weight, rotate(a x y)), circle, rect, line
(stroke-dasharray), path (M L H V Q C Z), fill / stroke / their opacities, inherited from the parents; <title>
(tooltips) is left out. Text is real text in the PDF (TrueType, editable), in Arial or Liberation Sans where
installed, DejaVu Sans for the glyphs they lack (→, ≤, …).
"""
import io, re
import xml.etree.ElementTree as ET
from matplotlib.figure import Figure
from matplotlib.patches import Circle, PathPatch, Rectangle
from matplotlib.lines import Line2D
from matplotlib.path import Path
from matplotlib import font_manager, rc_context

NS = "{http://www.w3.org/2000/svg}"
TEXT = "#222222"  # currentColor: the page's text color in light mode
INHERIT = ("fill", "stroke", "stroke-width", "fill-opacity", "stroke-opacity", "font-size", "font-weight", "text-anchor")


def _fonts():
    have = {f.name for f in font_manager.fontManager.ttflist}
    return [f for f in ("Arial", "Liberation Sans") if f in have] + ["DejaVu Sans"]


def _color(v, opacity=1.0):
    if v is None or v == "none":
        return None
    if v == "currentColor":
        v = TEXT
    m = re.fullmatch(r"rgba?\(([^)]*)\)", v.strip())
    if m:
        p = [float(x) for x in m.group(1).split(",")]
        return (p[0] / 255, p[1] / 255, p[2] / 255, (p[3] if len(p) > 3 else 1.0) * opacity)
    h = v.strip().lstrip("#")
    if len(h) == 3:
        h = "".join(c * 2 for c in h)
    if re.fullmatch(r"[0-9a-fA-F]{6}", h):
        return (int(h[0:2], 16) / 255, int(h[2:4], 16) / 255, int(h[4:6], 16) / 255, opacity)
    return (0, 0, 0, opacity)


def _path(d):
    toks = re.findall(r"[MLHVQCZmlhvqcz]|-?[\d.]+(?:e-?\d+)?", d)
    verts, codes, i, cmd, x, y = [], [], 0, None, 0.0, 0.0
    num = lambda: float(toks[i])
    while i < len(toks):
        if toks[i].isalpha():
            cmd = toks[i].upper(); i += 1
            if cmd == "Z":
                verts.append(verts[0] if verts else (x, y)); codes.append(Path.CLOSEPOLY)
                continue
        n = {"M": 2, "L": 2, "H": 1, "V": 1, "Q": 4, "C": 6}[cmd]
        a = [float(t) for t in toks[i:i + n]]; i += n
        if cmd == "H":
            x = a[0]; verts.append((x, y)); codes.append(Path.LINETO)
        elif cmd == "V":
            y = a[0]; verts.append((x, y)); codes.append(Path.LINETO)
        elif cmd in ("M", "L"):
            x, y = a; verts.append((x, y)); codes.append(Path.MOVETO if cmd == "M" else Path.LINETO)
            cmd = "L" if cmd == "M" else cmd  # numbers after M continue as lines
        else:
            pts = list(zip(a[::2], a[1::2]))
            verts += pts; codes += [Path.CURVE3 if cmd == "Q" else Path.CURVE4] * len(pts)
            x, y = pts[-1]
    return Path(verts, codes)


def to_pdf(svg):
    root = ET.fromstring(svg)
    W, H = float(root.get("width")), float(root.get("height"))
    fig = Figure(figsize=(W / 72, H / 72))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.axis("off")
    fonts, z = _fonts(), [0]

    def style(el, up):
        s = dict(up)
        s.update({k: el.get(k) for k in INHERIT if el.get(k) is not None})
        return s

    def paint(s):
        fill = _color(s.get("fill", "#000"), float(s.get("fill-opacity", 1)))
        stroke = _color(s.get("stroke"), float(s.get("stroke-opacity", 1)))
        return dict(facecolor=fill or "none", edgecolor=stroke or "none",
                    linewidth=float(s.get("stroke-width", 1)) if stroke else 0, clip_on=False)

    def walk(el, up):
        tag = el.tag.replace(NS, "")
        if tag == "title":
            return
        s = style(el, up)
        z[0] += 1
        f = lambda k: float(el.get(k, 0))
        if tag == "circle":
            ax.add_patch(Circle((f("cx"), f("cy")), f("r"), zorder=z[0], **paint(s)))
        elif tag == "rect":
            ax.add_patch(Rectangle((f("x"), f("y")), f("width"), f("height"), zorder=z[0], **paint(s)))
        elif tag == "path":
            ax.add_patch(PathPatch(_path(el.get("d")), zorder=z[0], **paint(s)))
        elif tag == "line":
            ln = Line2D([f("x1"), f("x2")], [f("y1"), f("y2")], color=_color(s.get("stroke"), float(s.get("stroke-opacity", 1))),
                        linewidth=float(s.get("stroke-width", 1)), zorder=z[0], clip_on=False)
            if el.get("stroke-dasharray"):
                ln.set_dashes([float(v) for v in re.split(r"[ ,]+", el.get("stroke-dasharray").strip())])
            ax.add_line(ln)
        elif tag == "text":
            text = (el.text or "") + "".join(c.tail or "" for c in el)
            m = re.search(r"rotate\(\s*(-?[\d.]+)", el.get("transform", ""))
            ax.text(f("x"), f("y"), text.strip(), fontsize=float(s.get("font-size", 11)), family=fonts,
                    fontweight="bold" if s.get("font-weight") in ("600", "700", "bold") else "normal",
                    ha={"middle": "center", "end": "right"}.get(s.get("text-anchor"), "left"), va="baseline",
                    rotation=-float(m.group(1)) if m else 0, rotation_mode="anchor",
                    color=_color(s.get("fill", TEXT)), zorder=z[0], clip_on=False)
        for c in el:
            walk(c, s)

    walk(root, {"fill": root.get("fill", "#000")})
    buf = io.BytesIO()
    with rc_context({"pdf.fonttype": 42}):
        fig.savefig(buf, format="pdf")
    return buf.getvalue()
