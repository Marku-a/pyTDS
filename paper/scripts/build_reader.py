"""Assemble reader.html (all aTDS reports in one page) from the paper's style and chart code + reader/main.html + reader/extra.js."""
import os
import re

P = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
paper = open(os.path.join(P, "atds_paper.html"), encoding="utf-8").read()
main = open(os.path.join(P, "reader", "main.html"), encoding="utf-8").read()
extra = open(os.path.join(P, "reader", "extra.js"), encoding="utf-8").read()

style = re.search(r"<style>.*?</style>", paper, re.S).group(0)
style = style.replace("</style>", """.sw { display: inline-block; width: 16px; height: 4px; border-radius: 2px; background: var(--c); vertical-align: middle; margin: 0 5px 0 2px; }
.sw.dash { background: repeating-linear-gradient(90deg, var(--c) 0 5px, transparent 5px 8px); }
.part { display: inline-block; font: 600 12px/1 var(--f-ui); letter-spacing: .08em; text-transform: uppercase; color: var(--accent); border: 1px solid var(--accent); border-radius: 3px; padding: 4px 7px; margin-right: 10px; vertical-align: 4px; }
.verdict { background: var(--tint); border-radius: 6px; padding: 14px 18px; margin: 10px 0 20px; }
.verdict p { margin: 4px 0 0; }
table.score td { font-size: 13.5px; }
</style>""")
js = paper[paper.index("<script>"):paper.index("function renderAll")]
hidden = '<div hidden><div id="fig1b"></div><div id="leg1b"></div><div id="fig4a"></div><div id="leg4a"></div></div>'
loader = """let R;
""" + extra + """
function renderAll() {
  [fig1, fig12, fig4, fig5, fig7, fig8, fig21, fig22, tab23, fig24, fig25, fig26, fig31, tabScore].forEach(fn => { try { fn(); } catch (e) { console.error(fn.name, e); } });
  if (E) [fig2, fig9].forEach(fn => { try { fn(); } catch (e) { console.error(fn.name, e); } });
}
Promise.all([fetch("figures/figdata.json").then(r => r.json()), fetch("figures/reader.json").then(r => r.json()),
             fetch("examples.json").then(r => r.ok ? r.json() : null).catch(() => null)])
  .then(([d, rr, e]) => { D = d; R = rr; E = e;
    if (!E) ["fig2a", "fig9a"].forEach(id => { const n = document.getElementById(id); if (n) n.innerHTML = '<div class="loading">Worked-example traces are not in the public repository; run paper/scripts/examples.py to generate examples.json.</div>'; });
    renderAll();
    let rt; addEventListener("resize", () => { clearTimeout(rt); rt = setTimeout(renderAll, 150); });
    const mq = matchMedia("(prefers-color-scheme: dark)"); (mq.addEventListener ? mq.addEventListener("change", renderAll) : mq.addListener(renderAll));
    new MutationObserver(renderAll).observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] }); })
  .catch(err => { document.querySelectorAll(".chart").forEach(n => n.innerHTML = '<div class="loading">Figure data could not be loaded.</div>'); console.error(err); });
})();
</script>"""
head = '<meta charset="utf-8">\n<title>aTDS Evidence Reader</title>\n' + paper[paper.index("<link rel=\"preconnect\""):paper.index("<style>")]
out = head + style + "\n" + main + "\n" + hidden + '\n<div class="tip" id="tip" hidden></div>\n' + js + loader
open(os.path.join(P, "reader.html"), "w", encoding="utf-8").write(out)
print("reader.html", len(out) // 1024, "KB")
