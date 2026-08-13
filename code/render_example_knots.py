"""
Render every `*.pd.txt` under preliminary_results/example_knots/ as a PNG of
the standard orthogonal link projection produced by plink.

Pipeline per file:
  PD code -> snappy.Link -> OrthogonalLinkDiagram -> plink.LinkEditor
          -> save_as_pdf -> pdftoppm -> PNG

The PDF is kept alongside the PNG. Requires plink (Tk-capable, so a display
is needed) and `pdftoppm` on PATH.

Run from the project root (in any environment where snappy + plink import):

    python code/render_example_knots.py
"""

import ast
import glob
import os
import shutil
import subprocess
import sys

import plink
import snappy
from spherogram.links.orthogonal import OrthogonalLinkDiagram

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
EXAMPLES_DIR = os.path.join(_ROOT, "preliminary_results", "example_knots")
DPI = 300


def _require_pdftoppm():
    if shutil.which("pdftoppm") is None:
        sys.exit(
            "error: pdftoppm not found on PATH. Install poppler "
            "(`brew install poppler`) and retry."
        )


def _load_pd(path):
    with open(path) as f:
        return ast.literal_eval(f.read().strip())


def _force_black(editor):
    """Recolor every arrow black. Used for single-component knots so they
    don't collide visually with the red R component of the RBG link."""
    for arrow in editor.Arrows:
        arrow.color = "#000000"
        try:
            editor.canvas.itemconfig(arrow.lines[0], fill="#000000") if arrow.lines else None
        except Exception:
            pass


def _render_link_to_pdf(link, pdf_path, force_black=False):
    editor = plink.LinkEditor()
    diagram = OrthogonalLinkDiagram(link)
    editor.unpickle(*diagram.plink_data())
    if force_black:
        _force_black(editor)
    try:
        editor.zoom_to_fit()
    except Exception:
        pass
    editor.save_as_pdf(pdf_path)
    try:
        editor.window.destroy()
    except Exception:
        pass


def _pdf_to_png(pdf_path, png_path):
    stem = png_path[:-4] if png_path.endswith(".png") else png_path
    subprocess.run(
        ["pdftoppm", "-r", str(DPI), "-png", "-singlefile", pdf_path, stem],
        check=True,
    )


def render_pd_file(pd_path):
    pd = _load_pd(pd_path)
    link = snappy.Link(pd)
    base = pd_path[: -len(".pd.txt")] if pd_path.endswith(".pd.txt") else os.path.splitext(pd_path)[0]
    pdf_path = base + ".pdf"
    png_path = base + ".png"
    force_black = len(link.link_components) == 1
    _render_link_to_pdf(link, pdf_path, force_black=force_black)
    _pdf_to_png(pdf_path, png_path)
    return png_path


def main():
    _require_pdftoppm()
    pd_files = sorted(glob.glob(os.path.join(EXAMPLES_DIR, "*", "*.pd.txt")))
    if not pd_files:
        print(f"No *.pd.txt files under {EXAMPLES_DIR}")
        return
    for pd_path in pd_files:
        rel = os.path.relpath(pd_path, _ROOT)
        try:
            out = render_pd_file(pd_path)
            print(f"  ✓ {rel} -> {os.path.relpath(out, _ROOT)}")
        except Exception as e:
            print(f"  ✗ {rel}: {type(e).__name__}: {e}")


if __name__ == "__main__":
    main()
