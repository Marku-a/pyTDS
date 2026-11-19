# Network-physiology journal article guidelines

Target: an Original Research article for a journal such as *Frontiers in Network Physiology*, *Physical Review E* or *PLoS ONE*.

## General setup
* Document class `article`, single column, 11 pt; keep the preamble of `template.tex` unchanged.
* Citations: numeric, `natbib` with `unsrtnat` (`\cite{key}`); use only keys from the reference library.
* Figures: single column, `\includegraphics[width=\textwidth]{figures/<file>.png}` (or `0.55\textwidth` for square plots).
* Tables: `booktabs`; wide tables in `tabularx` or `\resizebox{\textwidth}{!}{...}`.

## Length
* Main text about 8,000–10,000 words (roughly 15–19 pages in this format), excluding references and appendix.
* Abstract: one paragraph, at most 300 words.

## Content rules
* This is a methods paper in network physiology / time-series analysis, not a machine-learning conference paper: write for physicists and physiologists. Define every symbol.
* Report effect sizes with 95 % intervals as given in the experimental log; never round in a way that changes a conclusion.
* A Discussion section with explicit Limitations and "Not tested" paragraphs is mandatory; every limitation in the experimental log must appear.
* Mandatory structure: (1) Introduction, (2) Related Work, (3) Idea and algorithm, (4) Results on simulated data, (5) Results on real data, (6) Discussion, (7) Conclusion. Keep simulated and real results in separate sections.
* No overclaiming: real-data results are mixed for sleep networks (no significant gain on the main stage contrast; weaker wake contrast); say so in the abstract. All TDS is computed with the author's pyTDS.
* Author: Asaf Markuza (single author, not anonymous). Affiliation to be added by the author.
