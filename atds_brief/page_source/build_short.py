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
import re
body = (HERE / "body.html").read_text()
blocks = {m.group(1): m.group(0) for m in re.finditer(r'<div class="anim">\s*<canvas id="(a\d)".*?</button></div>\s*</div>', body, re.S)}
anim_css = re.search(r"\.anim \{.*?\.anim-bar button:focus-visible \{[^}]*\}", body, re.S).group(0)
for k in ("a1", "a2", "a3"):
    page = page.replace("__" + k.upper() + "__", blocks[k])
page = page.replace("__ANIMCSS__", anim_css).replace("__ANIM__", (HERE / "anim.js").read_text())
assert not re.search(r"__[A-Z0-9]+__", page), re.findall(r"__[A-Z0-9]+__", page)
(HERE.parent / "atds_brief_short.html").write_text(page)
print("wrote", HERE.parent / "atds_brief_short.html", len(page), "chars")
