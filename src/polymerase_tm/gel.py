"""
Virtual Agarose Gel Visualization.

Provides functions to plot a simulated agarose gel showing a DNA ladder
and the predicted amplicon size using matplotlib and seaborn.
"""

from __future__ import annotations

import os
from typing import Dict, List, Tuple

try:
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.patches import Rectangle
    import numpy as np
    _HAS_VIZ = True
except ImportError:
    _HAS_VIZ = False


# Standard NEB DNA Ladders (Size in bp, relative intensity for visual thickness)
LADDERS: Dict[str, List[Tuple[int, float]]] = {
    "1kb_plus": [
        (10000, 1.0),
        (8000, 1.0),
        (6000, 1.0),
        (5000, 1.0),
        (4000, 1.0),
        (3000, 3.0),  # Increased intensity reference band
        (2000, 1.0),
        (1500, 1.0),
        (1200, 1.0),
        (1000, 3.0),  # Increased intensity reference band
        (900, 1.0),
        (800, 1.0),
        (700, 1.0),
        (600, 1.0),
        (500, 3.0),   # Increased intensity reference band
        (400, 1.0),
        (300, 1.0),
        (200, 1.0),
        (100, 1.0),
    ],
    "100bp": [
        (1517, 1.0),
        (1200, 1.0),
        (1000, 3.0),  # Increased intensity reference band
        (900, 1.0),
        (800, 1.0),
        (700, 1.0),
        (600, 1.0),
        (500, 3.0),   # Increased intensity reference band
        (400, 1.0),
        (300, 1.0),
        (200, 1.0),
        (100, 1.0),
    ]
}


def _get_migration_distance(bp: int, max_bp: int, min_bp: int) -> float:
    """Calculate relative migration distance (logarithmic scale).
    Lower bp migrate further (higher distance).
    """
    # Safeguard against log(0)
    bp = max(1, bp)
    min_bp = max(1, min_bp)
    
    log_max = np.log10(max_bp)
    log_min = np.log10(min_bp)
    log_bp = np.log10(bp)
    
    # Invert so smaller fragments have higher values (migrate further down)
    # Scale from 0 (top/wells) to 1 (bottom)
    distance = (log_max - log_bp) / (log_max - log_min)
    return distance


def plot_virtual_gel(
    amplicon_length: int,
    ladder_name: str = "1kb_plus",
    output_path: str = "virtual_gel.png",
) -> bool:
    """
    Plot a simulated agarose gel with a DNA ladder and the target amplicon.
    
    Args:
        amplicon_length: The size of the expected amplicon in bp.
        ladder_name: The internal name of the ladder (e.g., '1kb_plus', '100bp').
        output_path: Where to save the resulting image file.
        
    Returns:
        bool: True if plotting succeeded, False if visualization libs are missing.
    """
    if not _HAS_VIZ:
        return False
        
    ladder = LADDERS.get(ladder_name, LADDERS["1kb_plus"])
    
    # Determine bounds for the gel visualization
    all_sizes = [b[0] for b in ladder] + [amplicon_length]
    max_bp = max(all_sizes) * 1.2  # Add 20% headroom at the top
    min_bp = min(all_sizes) * 0.8  # Add 20% space at the bottom (or fixed min)
    if ladder_name == "100bp":
        min_bp = 50
    else:
        min_bp = min(50, min_bp)
        

    # Set up the plot with a dark aesthetic resembling a UV-illuminated EtBr/GelRed gel
    sns.set_theme(style="dark", rc={"axes.facecolor": "#111111", "figure.facecolor": "#222222"})
    fig, ax = plt.subplots(figsize=(4, 8))
    
    # Draw lanes (subtle background rectangles)
    ax.add_patch(Rectangle((0.5, 0), 1, 1, facecolor="#1a1a1a", alpha=0.5, zorder=1))
    ax.add_patch(Rectangle((2.5, 0), 1, 1, facecolor="#1a1a1a", alpha=0.5, zorder=1))
    
    # Draw Wells at the top (distance 0)
    ax.add_patch(Rectangle((0.6, -0.05), 0.8, 0.04, facecolor="#000000", edgecolor="#333333", zorder=2))
    ax.add_patch(Rectangle((2.6, -0.05), 0.8, 0.04, facecolor="#000000", edgecolor="#333333", zorder=2))


    # Setup color modifiers for glowing effect
    glow_color = "#33ffcc"  # Cyan-ish glow
    band_color = "#e6ffff"  # Bright center
    
    # Plot Ladder (Lane 1)
    lane1_x = 1.0
    band_width = 0.6
    
    for bp, intensity in ladder:
        dist = _get_migration_distance(bp, max_bp, min_bp)
        
        # Base thickness adjusted by intensity relative to distance scaling
        thickness = 0.005 * intensity
        
        # Draw glowing halo
        ax.add_patch(Rectangle((lane1_x - band_width/2, dist - thickness*1.5), band_width, thickness*3, 
                               facecolor=glow_color, alpha=0.3 + (intensity*0.1), zorder=3))
        # Draw bright core
        ax.add_patch(Rectangle((lane1_x - band_width/2, dist - thickness/2), band_width, thickness, 
                               facecolor=band_color, alpha=0.9, zorder=4))
                               
        # Add text labels for reference bands or every alternating band
        if intensity > 1.0 or bp in [10000, 100]:
            ax.text(lane1_x - band_width/2 - 0.1, dist, f"{bp/1000:.1f}k" if bp >= 1000 else f"{bp}", 
                    color="#cccccc", ha="right", va="center", fontsize=9, fontweight="bold" if intensity > 1.0 else "normal")


    # Plot Amplicon (Lane 2)
    lane2_x = 3.0
    dist_amp = _get_migration_distance(amplicon_length, max_bp, min_bp)
    
    # Strong bright band for the amplicon
    thickness_amp = 0.01
    ax.add_patch(Rectangle((lane2_x - band_width/2, dist_amp - thickness_amp*1.5), band_width, thickness_amp*3, 
                           facecolor=glow_color, alpha=0.6, zorder=3))
    ax.add_patch(Rectangle((lane2_x - band_width/2, dist_amp - thickness_amp/2), band_width, thickness_amp, 
                           facecolor="#ffffff", alpha=1.0, zorder=4))
                           
    # Label the amplicon
    ax.text(lane2_x + band_width/2 + 0.1, dist_amp, f"Expected:\n{amplicon_length} bp", 
            color="#ffffff", ha="left", va="center", fontsize=10, fontweight="bold")

    
    # Configure axes
    ax.set_xlim(0, 4)
    # Y-axis goes from -0.1 (top, above wells) to 1.05 (bottom)
    ax.set_ylim(1.05, -0.1)
    
    ax.set_xticks([1.0, 3.0])
    ax.set_xticklabels(["Ladder\n(NEB 1kb+)", "Sample"], color="white", fontsize=11)
    
    ax.set_yticks([]) # Hide y-ticks, we rely on labels
    
    # Cleanup spines
    for spine in ax.spines.values():
        spine.set_visible(False)
        
    plt.title(f"Virtual Gel Simulation", color="white", pad=20, fontsize=14)
    
    # Save the plot
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close()
    
    return True
