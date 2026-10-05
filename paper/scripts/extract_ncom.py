"""Decode the thesis's D1 data (35 subjects x 10 signals, 1 Hz, sleep stages) from MATLAB table objects.

Input : <msc-project>/01 repreduce NCOM/surrogate test/All_data_clean.mat (MAT v5; struct p_data_stucrt
        with fields p_name, data (MATLAB `table`), time_add, sleep_stages).
Output: an .npz OUTSIDE git with X_<name> (T x 10 float) and S_<name> (T str stages).

MATLAB `table` objects live in the file's __function_workspace__ (MCOS subsystem), which scipy loads as
raw bytes. We parse that byte stream as a MAT v5 variable; the table payloads are the (10,) cells of
column vectors, in subject order, alternating with (10,) cells of column names. The subject order is
verified by matching every subject's data length to its sleep-stage length.

usage: python extract_ncom.py <msc-project dir> <out.npz>
"""
import io
import sys

import numpy as np
import scipy.io as sio
from scipy.io.matlab._mio5 import MatFile5Reader

SIGNALS = ["HR", "Resp", "Chin", "Leg", "Eye", "delta", "theta", "alpha", "sigma", "beta"]


def main(repo: str, out: str) -> None:
    fn = f"{repo}/01 repreduce NCOM/surrogate test/All_data_clean.mat"
    raw = sio.loadmat(fn)["__function_workspace__"].tobytes()
    f = io.BytesIO(b"MATLAB 5.0 MAT-file".ljust(116, b" ") + b"\x00" * 8 + raw[:4] + raw[8:])
    r = MatFile5Reader(f, struct_as_record=False, squeeze_me=True)
    r.mat_stream.seek(0)
    r.initialize_read()
    r.mat_stream.seek(128)
    hdr, _ = r.read_var_header()
    arr = r.read_var_array(hdr, process=True).MCOS[0]["arr"]

    cells = [x for x in arr.flat if getattr(x, "shape", None) == (10,)]
    data = [np.column_stack([np.asarray(c, float) for c in x]) for x in cells if isinstance(x[0], np.ndarray)]
    names = [list(x) for x in cells if isinstance(x[0], str)]
    assert all(n == SIGNALS for n in names), names[:1]

    d = sio.loadmat(fn, squeeze_me=True, struct_as_record=False)["p_data_stucrt"]
    stages = [np.array([e if isinstance(e, str) else "none" for e in s]) for s in d.sleep_stages]
    assert len(data) == len(stages) == 35
    assert all(len(x) == len(s) for x, s in zip(data, stages)), "subject order mismatch"
    assert not any(np.isnan(x).any() for x in data)

    out_d = {f"X_{n}": x for n, x in zip(d.p_name, data)}
    out_d.update({f"S_{n}": s for n, s in zip(d.p_name, stages)})
    np.savez_compressed(out, names=np.array(list(d.p_name)), signals=np.array(SIGNALS), **out_d)
    print(f"wrote {out}: {len(data)} subjects, {sum(len(x) for x in data) / 3600:.1f} h")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
