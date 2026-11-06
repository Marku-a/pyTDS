# Adaptive TDS (aTDS) — idea brief

A self-contained brief on **adaptive Time Delay Stability**: TDS that sets its
own window, lag range and tolerance from the signals. The only mandatory input
is the signals; the longest plausible delay `d_max` is optional.

- **Technical report:** `atds_brief.html` (open in any browser; interactive charts and
  animations, no internet needed except for fonts). Published privately at
  https://claude.ai/artifact/Ji7ju8xZzQQTRVyRsDFyP8 — share it from the page's
  Share menu.
- **Short brief (send this):** `atds_brief_short.html`, published at
  https://claude.ai/artifact/T2wbNs2czosxrEtzzbx3H3 — the idea, the step-by-step
  rules with animations, why each number was chosen, the key results and limits.
- **Review record:** `REVIEW.md` — five review rounds by an AI model (Claude)
  acting as a network-physiology referee, word for word, with the changes made
  after each. Final decision: approved to show as an idea after corrections.

## The rule (aTDS v2.1)

| Parameter | Rule |
|---|---|
| Window L | `30 · max_i B_i` (B = Bartlett factor, `1 + 2 Σ r(k)²` up to the first zero crossing), at least 61 samples (or `2·d_max + 1` if larger), at most 10 % of the record; the minimum wins over the cap |
| Step | `L // 2` |
| Lag range | `(L − 1) // 2`, capped at half the shortest detected rhythm period |
| Tolerance | largest tol with mean + 2 SE of the chance score ≤ 5 %, chance score from circularly rotated copies pooled over all available recordings |
| Stability rule | unchanged: 4 of 5 windows |

The base code is `experiments/atds_code/adaptive.py` (copied from branch
`claude/practical-mccarthy-n2mec6`, where it was written). The v2.1 changes
(minimum 61, pooled upper-bound calibration, optional `d_max`) are implemented
in `experiments/exp4.py` (`calib`, `atds2`) and `experiments/exp6.py`
(`atds_opt`).

## Folder layout

| Path | Contents |
|---|---|
| `atds_brief.html` | The full technical report |
| `atds_brief_short.html` | The short brief (build: `python atds_brief/page_source/build_short.py`) |
| `page_source/` | Page source: `body.html` (text, charts code), `style.css`, `lib.js` (chart helpers), `anim.js` (step animations), `build_page.py` |
| `experiments/` | Every script behind a figure or table (see below) |
| `experiments/atds_code/` | Copies of `adaptive.py`, `generators.py`, `common.py` used by the scripts |
| `results/` | The raw outputs (JSON) the page is built from |
| `REVIEW.md` | Referee reports |

## Which script produced what

| Script / tasks | Result file | Used in the brief |
|---|---|---|
| `theory.py` | `theory.json` | Fig. 2 (independent samples per window) |
| `tests.py` (test A) | `tests.json` | first version of Fig. 3 (kept for reference) |
| `tests2.py` | `tests2.json` | first versions of Figs 4 and 6 (kept for reference) |
| `exp3.py` T1–T9 | `exp3.json` | round-1 revision (kept for reference) |
| `exp4.py` U1 | `exp4.json` | Figs 6 and 7 (100 + 100 pairs, intervals, baselines, oracle) |
| `exp4.py` U2, U2m | `exp4.json` | Fig. 3 |
| `exp4.py` U3 | `exp4.json` | Fig. 3b |
| `exp4.py` U4 | `exp4.json` | Fig. 4 |
| `exp4.py` U5–U10 | `exp4.json` | Table 2 (long memory, rhythm, EEG, mixed speeds, changing dynamics) |
| `exp4.py` U12 | `exp4.json` | Figs 1 and 5 |
| `exp5.py` V1 | `exp5.json` | Table 2, Rulkov row |
| `exp5.py` V2 | `exp5.json` | Table 3 (networks, 999 rotations, FDR, common driver) |
| `exp6.py` W1 | `exp6.json` | Table 2, fast signal with long delay and the optional `d_max` |

## How to rebuild the page

```bash
python atds_brief/page_source/build_page.py
```

This fills `page_source/body.html` with the styles, scripts and the files in
`results/`, and writes `atds_brief.html`.

## How to rerun the experiments

The scripts need a checkout that contains both `pyTDS/` and
`experiments/adaptive/` (branch `claude/practical-mccarthy-n2mec6`), and they
import each other from the `experiments/` folder. Example:

```bash
git worktree add ../atds-code origin/claude/practical-mccarthy-n2mec6
cd atds_brief/experiments
python -I exp4.py ../../../atds-code . ../results/exp4.json U1,U2,U2m   # tasks to run, comma-separated
python -I exp5.py ../../../atds-code . ../results/exp5.json V1,V2
python -I exp6.py ../../../atds-code . ../results/exp6.json
```

Arguments: `<checkout with pyTDS and experiments/adaptive> <folder of these scripts> <output json> <tasks>`.
All random seeds are fixed, so reruns reproduce the saved results. The scripts
use 4 worker processes; `exp4.py` U1 and `exp5.py` V2 are the slowest tasks.
The EEG test reads `data/data-delta.txt` and `data/data-sigma.txt` from this
repository.

## Status

Approved as an idea by the AI referee (see `REVIEW.md`); not yet reviewed by a human expert. Not yet tested: staged multi-subject sleep
cohorts with cross-subject, stage-matched surrogates; coupling that switches
with the regime; stronger common drivers. See "Not tested here" and "Next
steps" in the brief.
