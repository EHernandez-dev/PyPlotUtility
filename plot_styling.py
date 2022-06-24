import matplotlib.pyplot as plt
import math
"""
    get_figure(aspect_ratio::Float64=1.1; x_pixels::Int64=600)

Get a `Figure` object with a given aspect ratio for a fixed number of pixels in x-diretion.
"""


def get_figure(aspect_ratio=1.1, x_pixels=600):
    plt.rc("figure", dpi=100)
    plt.figure(figsize=(aspect_ratio*1.e-2*x_pixels, 1.e-2*x_pixels))#, constrained_layout=true)


"""
    plot_styling!(;axis_label_font_size::Int64=20,
                        legend_font_size::Int64=15,
                        tick_label_size::Int64=15)

LMB default plot styling
"""

def plot_styling( x_pixels = 600,\
                        axis_label_font_size = 16,\
                        legend_font_size = 15,\
                        color="k"):

    axis_label_font_size = math.floor(x_pixels / 600 * axis_label_font_size)
    legend_font_size = math.floor(x_pixels / 600 * legend_font_size)

    plt.rc("font", size = axis_label_font_size, family = "stixgeneral")
    plt.rc("legend", fontsize = legend_font_size)
    plt.rc("mathtext", fontset = "stix")
    plt.rc("text", color=color)
    plt.rc("axes", labelcolor=color)
    plt.rc("xtick", color=color)
    plt.rc("ytick", color=color)


"""
    axis_ticks_styling!(ax::PyCall.PyObject; size_minor_ticks::Int64=6, 
                        width_minor_ticks::Int64=1, major_tick_width_scale::Int64=1,
                        tick_label_size::Int64=15, color::String="k")

LMB default axis tick styling.
"""
def axis_ticks_styling(ax, size_minor_ticks=3, 
                             width_minor_ticks=1, major_tick_width_scale=1,
                             tick_label_size=15, color="k"):

    ax.tick_params(reset=True, direction="in", axis="both", labelsize=tick_label_size,
                    which="major", size=size_minor_ticks<<1, 
                    width=major_tick_width_scale*width_minor_ticks, color=color)

    ax.tick_params(reset=True, direction="in", axis="both", labelsize=tick_label_size,
                    which="minor", size=size_minor_ticks, width=width_minor_ticks, color=color)

    ax.minorticks_on()

    for spine in ax.spines:
        spine.set_edgecolor(color)



    return ax



"""
    cb_ticks_styling!(ax::PyCall.PyObject; size_minor_ticks::Int64=6, 
                        width_minor_ticks::Int64=1, major_tick_width_scale::Int64=1,
                        tick_label_size::Int64=15, color::String="k")

LMB default colorbar tick styling.
"""
def cb_ticks_styling(ax, size_minor_ticks=3, 
                             width_minor_ticks=1, major_tick_width_scale=1,
                             tick_label_size=15, color="k"):

    ax.tick_params( direction="in", labelsize=tick_label_size,
                    which="major", size=size_minor_ticks<<1, 
                    width=major_tick_width_scale*width_minor_ticks, color=color)

    ax.tick_params( direction="in", labelsize=tick_label_size,
                    which="minor", size=size_minor_ticks, width=width_minor_ticks, color=color)

    ax.minorticks_on()

    for spine in ax.spines:
        spine.set_edgecolor(color)


    return ax



"""
    pixel_size(fig::Figure)

Helper function to get markers with size 1 pixel.
"""
def pixel_size(fig):
    return (72.0/fig.dpi)*(72.0/fig.dpi)
