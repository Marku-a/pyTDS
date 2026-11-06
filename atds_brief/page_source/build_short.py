"""Build the short brief atds_brief_short.html from short_body.html and the result files.

Usage (from the repository root):  python atds_brief/page_source/build_short.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RES = HERE.parent / "results"
E4 = json.load(open(RES / "exp4.json"))
dump = lambda o: json.dumps(o, separators=(",", ":"))

page = (HERE / "short_body.html").read_text()
page = page.replace("__STYLE__", (HERE / "style.css").read_text()).replace("__LIB__", (HERE / "lib.js").read_text())
page = page.replace("__EX__", dump(E4["U12"][0]["example"])).replace("__U4__", dump(E4["U4"])).replace("__U1__", dump(E4["U1"]))
assert "__EX__" not in page and "__U1__" not in page
(HERE.parent / "atds_brief_short.html").write_text(page)
print("wrote", HERE.parent / "atds_brief_short.html", len(page), "chars")
