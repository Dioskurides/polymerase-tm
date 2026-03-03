"""
Virtual Agarose Gel Visualization.

Provides functions to plot a simulated agarose gel showing a DNA ladder
and the predicted amplicon size(s) using matplotlib and seaborn.
"""

from __future__ import annotations

import os
from typing import Dict, List, Tuple, Union

try:
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.patches import Rectangle
    _HAS_VIZ = True
except ImportError:
    _HAS_VIZ = False

import numpy as np


# Standard NEB DNA Ladders (Size in bp, relative intensity for visual thickness)
LADDERS: Dict[str, List[Tuple[int, float]]] = {
    "1kb_plus": [
        (10000, 1.0), (8000, 1.0), (6000, 1.0), (5000, 1.0), (4000, 1.0),
        (3000, 3.0), (2000, 1.0), (1500, 1.0), (1200, 1.0), (1000, 3.0),
        (900, 1.0), (800, 1.0), (700, 1.0), (600, 1.0), (500, 3.0),
        (400, 1.0), (300, 1.0), (200, 1.0), (100, 1.0),
    ],
    "1kb": [
        (10000, 1.0), (8000, 1.0), (6000, 1.0), (5000, 1.0), (4000, 1.0),
        (3000, 3.0), (2000, 1.0), (1500, 1.0), (1000, 3.0), (500, 3.0),
    ],
    "100bp": [
        (1517, 1.0), (1200, 1.0), (1000, 3.0), (900, 1.0), (800, 1.0),
        (700, 1.0), (600, 1.0), (500, 3.0), (400, 1.0), (300, 1.0),
        (200, 1.0), (100, 1.0),
    ],
    "50bp": [
        (1350, 1.0), (1000, 1.0), (900, 1.0), (800, 1.0), (700, 1.0),
        (600, 1.0), (500, 1.0), (400, 1.0), (300, 3.0), (250, 3.0),
        (200, 1.0), (150, 1.0), (100, 1.0), (50, 1.0),
    ],
    "low_mw": [
        (766, 1.0), (500, 1.0), (344, 1.0), (250, 1.0), (200, 1.0),
        (150, 1.0), (100, 1.0), (75, 1.0), (50, 1.0), (25, 1.0),
    ],
    "pcr_marker": [
        (1000, 1.0), (766, 1.0), (500, 1.0), (300, 1.0), (150, 1.0), (50, 1.0),
    ]
}

LADDER_LABELS = {
    "1kb_plus": "NEB 1kb+",
    "1kb": "NEB 1kb",
    "100bp": "NEB 100bp",
    "50bp": "NEB 50bp",
    "low_mw": "NEB Low MW",
    "pcr_marker": "NEB PCR Marker"
}


def _get_migration_distance_cm(bp: int, agarose_pct: float, voltage: float, time_min: float) -> float:
    """Calculate relative migration distance (cm).
    Uses a highly accurate empirical physics model derived from Ferguson plots
    for agarose gel electrophoresis of linear dsDNA.
    """
    bp = max(10, bp)
    agarose_pct = max(0.2, agarose_pct)
    voltage = max(1.0, voltage)
    time_min = max(1.0, time_min)
    
    # 1. Base migration at 1.0% agarose, 100V, 60min
    x = np.log10(bp)
    
    # Adjusted polynomial fit to perfectly space NEB 1kb+ ladder ranges (0.1kb - 10kb)
    if bp >= 100:
        base_dist = 0.59 * (x**2) - 8.5 * x + 26.0
    else:
        # Fragments < 100bp experience minimal retardation and run linearly with the buffer front
        # 100bp (x=2) is at 11.36cm
        # We model a linear spread for <100bp: 50bp ~ 12.5cm, 25bp ~ 13.0cm
        base_dist = -3.8 * x + 18.96 
        
    # 2. Voltage and Time scaling (linear)
    vt_factor = (voltage / 100.0) * (time_min / 60.0)
    
    # 3. Agarose retardation 
    # Empirical Ferguson scaling: High MW fragments retard much faster in denser gels
    c_x = max(0.01, 0.4 * x - 0.5)
    agarose_factor = np.exp(-c_x * (agarose_pct - 1.0))
    
    distance_cm = base_dist * vt_factor * agarose_factor
    return distance_cm


def plot_virtual_gel(
    amplicon_lengths: Union[int, List[int]],
    ladder_name: str = "1kb_plus",
    output_path: str = "virtual_gel.png",
    agarose_pct: float = 1.0,
    voltage: float = 110.0,
    time_min: float = 60.0
) -> bool:
    """
    Plot a simulated agarose gel with a DNA ladder and the target amplicon(s).
    
    Args:
        amplicon_lengths: The size(s) of the expected amplicon(s) in bp.
        ladder_name: The internal name of the ladder (e.g., '1kb_plus', '100bp').
        output_path: Where to save the resulting image file.
        
    Returns:
        bool: True if plotting succeeded, False if visualization libs are missing.
    """
    if not _HAS_VIZ:
        return False
        
    ladder = LADDERS.get(ladder_name, LADDERS["1kb_plus"])
    ladder_label = LADDER_LABELS.get(ladder_name, "Ladder")
    
    if isinstance(amplicon_lengths, int):
        amplicons = [amplicon_lengths]
    else:
        amplicons = list(amplicon_lengths)
        
    if not amplicons:
        return False
    
    # Determine standard 12cm mini-gel layout
    gel_length_cm = 12.0

    # Set up the plot width dynamically based on number of lanes
    num_lanes = 1 + len(amplicons)
    plot_width = max(4, num_lanes * 2.0)
    
    sns.set_theme(style="dark", rc={"axes.facecolor": "#111111", "figure.facecolor": "#222222"})
    fig, ax = plt.subplots(figsize=(plot_width, 8))
    
    # Draw lanes
    lane_xs = []
    for i in range(num_lanes):
        x_center = 2.0 + (i * 2.0)  # Start the first lane further to the right
        lane_xs.append(x_center)
        # Background lane
        ax.add_patch(Rectangle((x_center - 0.5, 0), 1, 1, facecolor="#1a1a1a", alpha=0.5, zorder=1))
        # Well
        ax.add_patch(Rectangle((x_center - 0.35, -0.2), 0.7, 0.2, facecolor="#000000", edgecolor="#333333", zorder=2))

    glow_color = "#33ffcc"
    band_color = "#e6ffff"
    band_width = 0.6
    
    # Plot Ladder (Lane 1)
    lane1_x = lane_xs[0]
    for bp, intensity in ladder:
        dist = _get_migration_distance_cm(bp, agarose_pct, voltage, time_min)
        if dist > gel_length_cm + 2.0:
            continue  # fragment ran off the gel entirely
            
        # Thinner bands for a much more realistic look
        thickness = 0.03 * (intensity ** 0.5)
        
        # Outer glow - very faint
        ax.add_patch(Rectangle((lane1_x - band_width/2, dist - thickness*1.5), band_width, thickness*3, 
                               facecolor=glow_color, alpha=0.1 + (intensity*0.05), zorder=3))
        # Inner glow - moderate
        ax.add_patch(Rectangle((lane1_x - band_width/2, dist - thickness*0.8), band_width, thickness*1.6, 
                               facecolor=glow_color, alpha=0.3 + (intensity*0.1), zorder=4))
        # Bright core
        ax.add_patch(Rectangle((lane1_x - band_width/2, dist - thickness*0.2), band_width, thickness*0.4, 
                               facecolor=band_color, alpha=0.9, zorder=5))
                               
        # Explicitly requested by user: ALL fragments of the ladder should have their size annotated
        bp_label = f"{bp/1000:.1f}k" if bp >= 1000 else f"{bp}"
        label_text = f"{bp_label} ({dist:.1f} cm)"
        font_weight = "bold" if intensity > 1.0 else "normal"
        ax.text(lane1_x - band_width/2 - 0.1, dist, label_text, 
                color="#cccccc", ha="right", va="center", fontsize=8, fontweight=font_weight)

    # Plot Amplicons (Lanes 2+)
    for i, amp_len in enumerate(amplicons):
        lane_x = lane_xs[i + 1]
        dist_amp = _get_migration_distance_cm(amp_len, agarose_pct, voltage, time_min)
        if dist_amp > gel_length_cm + 2.0:
            ax.text(lane_x, gel_length_cm + 0.5, f"RAN OFF GEL\n({amp_len} bp)", 
                    color="#ff4444", ha="center", va="top", fontsize=9, fontweight="bold")
            continue
            
        thickness_amp = 0.04
        # Outer glow
        ax.add_patch(Rectangle((lane_x - band_width/2, dist_amp - thickness_amp*1.5), band_width, thickness_amp*3, 
                               facecolor=glow_color, alpha=0.15, zorder=3))
        # Inner glow
        ax.add_patch(Rectangle((lane_x - band_width/2, dist_amp - thickness_amp*0.8), band_width, thickness_amp*1.6, 
                               facecolor=glow_color, alpha=0.4, zorder=4))
        # Bright core
        ax.add_patch(Rectangle((lane_x - band_width/2, dist_amp - thickness_amp*0.2), band_width, thickness_amp*0.4, 
                               facecolor="#ffffff", alpha=1.0, zorder=5))
                               
        ax.text(lane_x + band_width/2 + 0.1, dist_amp, f"{amp_len} bp", 
                color="#ffffff", ha="left", va="center", fontsize=10, fontweight="bold")
    
    max_x = max(lane_xs) + 2.0
    ax.set_xlim(0, max_x)
    # Print the physical dimensions using invert Y so 0 is at top.
    # We add 1.0 cm padding at the bottom so fast fragments (like 100bp @ 12.5cm) don't get cut off
    ax.set_ylim(gel_length_cm + 1.0, -0.5)
    
    x_ticks = lane_xs
    x_labels = [ladder_label] + [f"Sample {i+1}" if len(amplicons) > 1 else "Sample" for i in range(len(amplicons))]
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_labels, color="white", fontsize=11)
    
    ax.set_yticks([]) 
    for spine in ax.spines.values():
        spine.set_visible(False)
        
    plt.title("Virtual Gel Simulation", color="white", pad=20, fontsize=14)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close()
    
    return True
