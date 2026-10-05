"""R1.1: thesis port on all 35 D1 subjects vs stored MATLAB real_links. Writes results/r1_real_links_check.json."""
import json, os, sys
from concurrent.futures import ProcessPoolExecutor
import numpy as np, scipy.io as sio
sys.path.insert(0, os.path.dirname(__file__))
from thesis_port import links_strength

DATA, REPO = sys.argv[1], sys.argv[2]
OUT = os.path.join(os.path.dirname(__file__), "..", "results")


def one(i):
    z = np.load(DATA); n = z["names"][i]
    return links_strength(z["X_" + n], z["S_" + n]).flatten(order="F")


if __name__ == "__main__":
    ref = sio.loadmat(f"{REPO}/01 repreduce NCOM/surrogate test/links_strength2000Sur.mat")["real_links"]
    with ProcessPoolExecutor() as ex:
        mine = np.array(list(ex.map(one, range(35))))
    d = np.abs(mine - ref)
    res = {"subjects": 35, "values": int(d.size), "max_abs_diff": float(np.nanmax(d)),
           "n_diff_gt_1e-9": int((d > 1e-9).sum()), "nan_mine": int(np.isnan(mine).sum()), "nan_ref": int(np.isnan(ref).sum())}
    os.makedirs(OUT, exist_ok=True)
    json.dump(res, open(os.path.join(OUT, "r1_real_links_check.json"), "w"), indent=1)
    np.save(os.path.join(os.environ.get("SCRATCH", "/tmp"), "real_links_port.npy"), mine)
    print(res)
