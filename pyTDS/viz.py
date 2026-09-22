"""
Visualization helpers for TDS results.

All functions return (fig, ax) so they can be embedded in larger figures,
or called standalone (they create their own figure if ax is not provided).
"""

from __future__ import annotations
import numpy as np


def plot_tau_series(
    tau: np.ndarray,
    t_vec: np.ndarray,
    stbl_lbl: np.ndarray = None,
    title: str = "Time Delay Interaction",
    ax=None,
):
    """
    Plot the time delay series τ₀(t) with stable regions highlighted in red.

    Reproduces Figure 1b style from Bashan et al. 2012 (Nature Communications).

    Parameters
    ----------
    tau : np.ndarray
        Time delay series from time_delay_interaction().
    t_vec : np.ndarray
        Timestamps (x-axis).
    stbl_lbl : np.ndarray, optional
        Binary stability labels. Stable points plotted in red, unstable in blue.
    title : str
    ax : matplotlib.axes.Axes, optional

    Returns
    -------
    (fig, ax)
    """
    import matplotlib.pyplot as plt

    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 3))

    if stbl_lbl is not None:
        # Handle case where stbl_lbl has different length (middle_ref method)
        if len(stbl_lbl) == len(tau):
            mask_stable = stbl_lbl == 1
            mask_unstable = ~mask_stable
            ax.plot(t_vec[mask_unstable], tau[mask_unstable],
                    'o', color='steelblue', markersize=4, alpha=0.6, label='Unstable')
            ax.plot(t_vec[mask_stable], tau[mask_stable],
                    'o', color='crimson', markersize=5, alpha=0.9, label='Stable')
        else:
            ax.plot(t_vec, tau, 'o--', color='steelblue', markersize=4, alpha=0.7)
    else:
        ax.plot(t_vec, tau, 'o--', color='steelblue', markersize=4, alpha=0.7)

    ax.axhline(0, color='gray', linewidth=0.5, linestyle='--')
    ax.set_xlabel("Time (samples)")
    ax.set_ylabel("Time Delay (samples)")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    if stbl_lbl is not None and len(stbl_lbl) == len(tau):
        ax.legend(loc='upper right', fontsize=8)

    if fig is not None:
        fig.tight_layout()

    return fig, ax


def plot_tds_matrix(
    mat: np.ndarray,
    labels: list[str] = None,
    threshold: float = None,
    title: str = "TDS Matrix",
    ax=None,
    cmap: str = "hot",
):
    """
    Plot a TDS matrix as a heatmap.

    Reproduces Figure 2 style from Bashan et al. 2012.

    Parameters
    ----------
    mat : np.ndarray
        Shape (N, N) — TDS score matrix.
    labels : list[str], optional
        Axis tick labels for each signal.
    threshold : float, optional
        If provided, draws a colorbar marker at the threshold level.
    title : str
    ax : matplotlib.axes.Axes, optional
    cmap : str
        Colormap. Default 'hot' matches the NCOM paper style.

    Returns
    -------
    (fig, ax)
    """
    import matplotlib.pyplot as plt

    N = mat.shape[0]
    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 6))

    im = ax.imshow(mat, cmap=cmap, vmin=0, vmax=100, aspect='auto')
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("TDS (%)")
    if threshold is not None:
        cbar.ax.axhline(threshold, color='cyan', linewidth=2, label=f'Threshold ({threshold:.1f}%)')

    if labels is not None:
        ax.set_xticks(range(N))
        ax.set_yticks(range(N))
        ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
        ax.set_yticklabels(labels, fontsize=8)

    ax.set_title(title)

    if fig is not None:
        fig.tight_layout()

    return fig, ax


def plot_network(
    graph,
    tds_mat: np.ndarray = None,
    labels: list[str] = None,
    threshold: float = None,
    title: str = "Physiological Network",
    ax=None,
):
    """
    Plot a physiological network as a circular diagram.

    Reproduces the network diagrams in Figure 2 of Bashan et al. 2012.
    Edge thickness and color represent link strength (TDS score).

    Parameters
    ----------
    graph : networkx.Graph
        From to_networkx().
    tds_mat : np.ndarray, optional
        TDS matrix for edge coloring. If None, uses graph edge weights.
    labels : list[str], optional
        Displayed node labels.
    threshold : float, optional
        For title annotation.
    title : str
    ax : matplotlib.axes.Axes, optional

    Returns
    -------
    (fig, ax)
    """
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors

    try:
        import networkx as nx
    except ImportError:
        raise ImportError("networkx is required: pip install networkx")

    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 7))

    nodes = list(graph.nodes())
    N = len(nodes)
    pos = nx.circular_layout(graph)

    # Node drawing
    nx.draw_networkx_nodes(graph, pos, ax=ax, node_size=600,
                           node_color='white', edgecolors='black', linewidths=1.5)
    nx.draw_networkx_labels(graph, pos, ax=ax, font_size=9, font_weight='bold')

    # Edge drawing — width and color by TDS weight
    edges = list(graph.edges(data=True))
    if edges:
        weights = np.array([d.get('weight', 1.0) for _, _, d in edges])
        norm = mcolors.Normalize(vmin=0, vmax=100)
        cmap_edges = cm.get_cmap('RdYlGn')
        edge_colors = [cmap_edges(norm(w)) for w in weights]
        edge_widths = 1.0 + 4.0 * (weights / 100.0)

        nx.draw_networkx_edges(graph, pos, ax=ax,
                               edge_color=edge_colors,
                               width=edge_widths,
                               alpha=0.85)

    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.axis('off')

    if fig is not None:
        fig.tight_layout()

    return fig, ax
