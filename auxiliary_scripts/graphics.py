#!/usr/bin/env python3.7

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns

from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, optimal_leaf_ordering
from scipy.cluster.hierarchy import leaves_list

import matplotlib.patches as mpatches

def dotplot(
    size_df,
    color_df,
    mask=None,
    size_legend_exp=False,
    size_scale=1.0,
    cmap='viridis',
    size_legend_title=None,
    size_legend_vmin=None, size_legend_vmax=None,
    padding=0.5,
    max_dot_fraction=0.8,
    figsize_xscale=0.5,
    figsize_yscale=0.5,
    center=None,
    title=None,
    xlabel=None,
    ylabel=None,
    clabel=None,
    slabel=None,
    xrotation=0,
    yrotation=0,
    xgroup=None,
    xgroup_rotation=0,
    ygroup=None,
    ygroup_rotation=0,
    edgecolor='none',
    edgewidth=1,
    show=True
):
    """
    Create a dotplot where circle size and color are determined by two DataFrames.
    Automatically rescales dots to avoid overlap.

    Parameters
    ----------
    size_df : pd.DataFrame
        DataFrame of numeric values controlling circle size.
    color_df : pd.DataFrame
        DataFrame of numeric values controlling circle color.
        Must have same shape and index/columns as size_df.
    mask : pd.DataFrame
        Remove values from plot.
        Must have same shape and index/columns as size_df.
    size_scale : float
        Additional scaling multiplier for dot size (applied after normalization).
    cmap : str
        Matplotlib colormap name for color encoding.
    size_legend_title : str
        Title for the size legend.
    padding : float
        Extra spacing around the plot edges (in data units).
    max_dot_fraction : float
        Fraction of the grid cell (0–1) that the largest dot should occupy.
        1.0 means circles just touch; <1 leaves spacing.
    """
    if not size_df.shape == color_df.shape:
        raise ValueError("size_df and color_df must have the same shape")
    if not (size_df.index.equals(color_df.index) and size_df.columns.equals(color_df.columns)):
        raise ValueError("size_df and color_df must have same index and columns")

    x_labels = size_df.columns
    y_labels = size_df.index
    x, y = np.meshgrid(np.arange(len(x_labels)), np.arange(len(y_labels)))
    x = x.flatten()
    y = y.flatten()

    sizes = size_df.values.flatten()
    colors = color_df.values.flatten()
    
    if mask is not None:
        idx = np.where(mask.values.flatten())
        sizes = sizes[idx]
        colors = colors[idx]
        x = x[idx]
        y = y[idx]

    fig, ax = plt.subplots(figsize=(len(x_labels)*figsize_xscale, len(y_labels)*figsize_yscale))

    # === Normalize sizes so they don't overlap ===
    # Circle area in matplotlib scatter = s (points^2)
    # To ensure largest circle fits inside cell width (1.0), convert to area accordingly.
    max_area = (72 * max_dot_fraction / fig.dpi) ** 2 * (fig.dpi / 72) ** 2 * 100  # base scaling
    norm_sizes = (sizes / sizes.max()) * max_area * size_scale
    if center is None:
        norm = None
    else:
        norm = mpl.colors.CenteredNorm(vcenter=center)

    sc = ax.scatter(x, y, s=norm_sizes, c=colors, cmap=cmap, edgecolor=edgecolor, lw=edgewidth, norm=norm)

    # Axis limits with padding
    ax.set_xlim(-padding, len(x_labels) - 1 + padding)
    ax.set_ylim(-padding, len(y_labels) - 1 + padding)

    # Axis labels
    ax.set_xticks(np.arange(len(x_labels)))
    ax.set_xticklabels(x_labels, rotation=xrotation)
    ax.set_yticks(np.arange(len(y_labels)))
    ax.set_yticklabels(y_labels, rotation=yrotation)
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    # Colorbar
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label(clabel)

    # --- Horizontal size legend ---
    if size_legend_vmin is None:
        size_legend_vmin = sizes[sizes != 0].min()
    if size_legend_vmax is None:
        size_legend_vmax = sizes.max()
    size_vals = np.linspace(size_legend_vmin, size_legend_vmax, 3)
    legend_handles = [
        plt.scatter([], [], s=(val / sizes.max()) * max_area * size_scale, color='gray', 
                    edgecolor=edgecolor, lw=edgewidth)
        for val in size_vals
    ]
    
    legend_labels = [f"{np.exp(-val):.2f}" if size_legend_exp else f"{val:.2f}" for val in size_vals]

    ax.legend(
        legend_handles,
        legend_labels,
        title=size_legend_title,
        loc='upper center',
        bbox_to_anchor=(0.5, -0.3),
        ncol=len(size_vals),
        frameon=True, borderpad=1,
    )

    
    # --- Add grouping lines/brackets ---
    def add_group_lines(axis, groups, is_x=True):
        """Draws lines and labels separating groups on the specified axis."""
        if not groups:
            return
        ticks = np.arange(len(x_labels) if is_x else len(y_labels))
        labels = x_labels if is_x else y_labels

        # Map each label to its tick index
        label_to_idx = {label: i for i, label in enumerate(labels)}

        for i, (group_name, members) in enumerate(groups.items()):
            idxs = sorted(label_to_idx[m] for m in members if m in label_to_idx)
            if not idxs:
                continue
            start, end = min(idxs), max(idxs)

            if is_x:
                # vertical lines and label above
                if i != len(groups)-1:
                    ax.axvline(end + 0.5, color='black', lw=1)
                ax.text((start + end) / 2, -padding - 0.5, group_name,
                        ha='left' if xgroup_rotation > 0 else 'center', 
                        va='bottom' if xgroup_rotation > 0 else 'top', 
                        fontsize=10, fontweight='bold', 
                        rotation=xgroup_rotation)
                ax.set_title('')
            else:
                # horizontal lines and label to the left
                if i != len(groups)-1:
                    ax.axhline(end + 0.5, color='black', lw=1)
                ax.text(-padding*15 - 0.5, (start + end) / 2, group_name,
                        ha='center' if ygroup_rotation > 0 else 'right', va='center', rotation=ygroup_rotation, fontsize=10, fontweight='bold')
    
    if xgroup is not None:
        add_group_lines(ax, xgroup, is_x=True)
    if ygroup is not None:
        add_group_lines(ax, ygroup, is_x=False)
    
    plt.tight_layout()
    
    if show:
        plt.show()
    else:
        return fig

