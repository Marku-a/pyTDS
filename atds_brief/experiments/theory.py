import sys, json
from dataclasses import replace
import numpy as np
W = sys.argv[1]
sys.path.insert(0, W + "/experiments/adaptive"); sys.path.insert(0, W)
from adaptive import adaptive_params, bartlett_factor
from generators import ar1, coupled_pair
from pyTDS.params import TDSParams
from pyTDS.core import time_delay_interaction, stable_label, tds_score
out = {}
# Fig 1: memory vs independent samples in a 60-sample window
rng = np.random.default_rng(1)
rows = []
for tc in [1, 2, 4, 8, 16, 32, 64]:
    Bs = [bartlett_factor(ar1(30000, tc, rng)) for _ in range(5)]
    B = float(np.median(Bs))
    rows.append(dict(tc=tc, B=round(B, 2), indep_fixed=round(60 / B, 2),
                     atds_window=int(round(30 * B)), indep_atds=round(round(30 * B) / B, 1)))
out["memory"] = rows
# Fig 1/5: worked example: tc=8, true delay 20, coupling c=0.6 (superseded by exp3.py T9, which sets parameters on a separate training pair)
rng = np.random.default_rng(7)
x, y, _ = coupled_pair(30000, 8, 8, 20, 0.6, rng)
fixed = TDSParams()
ad, info = adaptive_params([x, y], window_rule="bartlett", tolerance="calibrated")
def run(p):
    tau, t, _ = time_delay_interaction(x, y, p)
    lab = stable_label(tau, p)
    return tau, t, lab
ex = {}
for name, p in [("fixed", fixed), ("atds", ad)]:
    tau, t, lab = run(p)
    m = t <= 12000
    ex[name] = dict(params=dict(window=p.window, step=p.overlap, max_lag=p.max_lag, tol=p.tolerance),
                    score=round(float(tds_score(lab)), 1),
                    near_true=round(float(np.mean(np.abs(tau + 20) <= 2)) * 100, 1),
                    t=t[m].astype(int).tolist(), tau=tau[m].astype(int).tolist(), stable=lab[m].astype(int).tolist())
out["example"] = ex
out["example_info"] = dict(bartlett=[round(b, 2) for b in info["bartlett"]], window=info["window"],
                           calibration_failed=bool(info["calibration_failed"]))
# full tolerance curve: null (circular shift) vs real
T = len(x); r2 = np.random.default_rng(3)
lo, hi = max(ad.window, T // 10), T - max(ad.window, T // 10)
null_taus = [time_delay_interaction(x, np.roll(y, int(r2.integers(lo, hi))), ad)[0] for _ in range(20)]
real_tau = time_delay_interaction(x, y, ad)[0]
curve = []
for tol in range(0, 31):
    p = replace(ad, tolerance=tol)
    curve.append(dict(tol=tol, null=round(float(np.mean([tds_score(stable_label(t, p)) for t in null_taus])), 2),
                      real=round(float(tds_score(stable_label(real_tau, p))), 2)))
out["tol_curve"] = curve
json.dump(out, open(sys.argv[2], "w"))
print(json.dumps(out["memory"])); print({k: (v["params"], v["score"], v["near_true"]) for k, v in ex.items()}); print(out["example_info"])
print([(c["tol"], c["null"], c["real"]) for c in curve])
