"""Rebuild atds_brief.html from body.html, the style/script parts and the result files.

Usage (from the repository root):  python atds_brief/page_source/build_page.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RES = HERE.parent / "results"


def data(name):
    return json.dumps(json.load(open(RES / name)), separators=(",", ":"))


page = (HERE / "body.html").read_text()
page = page.replace("__ANIM__", (HERE / "anim.js").read_text())
page = page.replace("__STYLE__", (HERE / "style.css").read_text()).replace("__LIB__", (HERE / "lib.js").read_text())
for key, fname in [("__DATA__", "theory.json"), ("__TESTS2__", "tests2.json"), ("__EXP3__", "exp3.json"),
                   ("__EXP4__", "exp4.json"), ("__EXP5__", "exp5.json"), ("__TESTS__", "tests.json")]:
    page = page.replace(key, data(fname))
assert "__" + "EXP" not in page and "__DATA__" not in page, "unfilled placeholder"
(HERE.parent / "atds_brief.html").write_text(page)
print("wrote", HERE.parent / "atds_brief.html", len(page), "bytes")
