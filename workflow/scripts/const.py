# TODO:
# •Fontsize, font (Ariel)
# •	Axes labels
# •	Colors for activity, diff activity etc.
# •	Top and right box boundaries off
# •	No grid
# •	White background
# •	Each scatterplot with density (choose new cmap but with lowest density as brightest color, closest to white)
# •	Legend, title and clear labels
# •	legend outside of the plot
# •	no figure title
# •	Savefig function: Save file as png and eps, 300 dpi
# •	Input paths

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.figure import Figure

# Standard font size and family
FONT_SIZE_small: int = 18  # changed from 16 NM 06-02-2025
FONT_SIZE_big: int = 20  # changed from 18 NM 06-02-2025

FONT_FAMILY: str = "arial"

# Define color maps for different data
pos_active_ctrl_color: str = "g"
neg_active_ctrl_color: str = "r"
highlight_color: str = "y"

DIFF_ACTIVITY_COLOR: str = "green"

# Custom colormap for scatterplots

colors: list[str] = ["#EBF4FF", "#E1ECFA", "#D0E0F5", "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026"]
custom_cmap: LinearSegmentedColormap = LinearSegmentedColormap.from_list("custom_cmap", colors, N=256)

colors: list[str] = ["#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026"]
custom_cmap_bolder: LinearSegmentedColormap = LinearSegmentedColormap.from_list("custom_cmap", colors, N=256)


# Set default figure style
def set_plot_style() -> None:
    """
    Set standardized figure settings for matplotlib.
    """
    plt.rcParams.update(
        {
            "axes.titlesize": FONT_SIZE_big,
            "axes.labelsize": FONT_SIZE_big,
            "xtick.labelsize": FONT_SIZE_small,
            "ytick.labelsize": FONT_SIZE_small,
            "legend.fontsize": FONT_SIZE_small,
            "legend.title_fontsize": FONT_SIZE_small,
            "axes.labelweight": "bold",
            #        'font.family': FONT_FAMILY, # TODO: change font
            "axes.linewidth": 1.0,
            "figure.figsize": (8, 8),
            "axes.grid": False,  # No grid
            "axes.spines.top": False,  # Top border off
            "axes.spines.right": False,  # Right border off
            "figure.facecolor": "none",
            "axes.facecolor": "none",
            #        'legend.frameon': False  # No frame for legend
        }
    )


#     _ = plt.figure(figsize=(9,9)) # TODO: adjust figure size to comply with fontsize


def save_fig(fig: Figure, name: str, path: str) -> None:
    """
    Save the figure to the specified path in PNG and EPS formats at 500 dpi.
    """
    fig.savefig(f"{path}/{name}.png", dpi=500, bbox_inches="tight", transparent=True)
    fig.savefig(f"{path}/{name}.svg", bbox_inches="tight", transparent=True)
    fig.savefig(f"{path}/{name}.eps", dpi=500, bbox_inches="tight")  # increasd DPI to 500 NM 19/09
    fig.savefig(f"{path}/{name}.pdf", dpi=500, bbox_inches="tight", transparent=True)


def set_equal_plot_limits(x: np.ndarray, y: np.ndarray) -> None:
    """
    Sets the x and y axis limits to the same range based on the min and max values of x and y.

    Parameters:
    x (array-like): Data for the x-axis.
    y (array-like): Data for the y-axis.
    """
    min_limit = min(np.min(x), np.min(y))
    max_limit = max(np.max(x), np.max(y))

    plt.xlim([min_limit, max_limit])
    plt.ylim([min_limit, max_limit])


plot_color_pallete: dict[str, str] = {
    "default_color": "#AEAEAE",
    "cCRE": "#3D9F95",  # orange
    "barcode": "#227C9D",  # turquoise
    "read": "#FFC25F",  # purple
    "cCRE-barcode-pair": "#9383B8",  # magenta
}
