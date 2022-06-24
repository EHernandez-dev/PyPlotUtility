import numpy as np
import matplotlib

"""
    get_imshow(ax::PyCall.PyObject, image::Array{<:Real}, 
                    x_lim::Array{<:Real}=zeros(2), y_lim::Array{<:Real}=zeros(2), 
                    vmin::Real=0.0, vmax::Real=0.0; 
                    cmap::String="viridis", cnorm=matplotlib.colors.NoNorm(),
                    ticks_color::String="white",
                    interpolation::String="none")

Helper function to plot an `imshow` with linear colorbar.
"""
def get_imshow(ax, image, x_lim=np.zeros(2), y_lim=np.zeros(2),  vmin=0.0, vmax=0.0, cmap="viridis", cnorm=matplotlib.colors.NoNorm(),ticks_color="white",interpolation="none"):

    if vmin == 0.0 and vmax == 0.0:
        vmin = min(image)
        vmax = max(image)


    if x_lim == np.zeros(2) and y_lim == np.zeros(2):
        x_lim = [1, len(image[:,1])]
        y_lim = [1, len(image[1,:])]

    
    im = ax.imshow(image, norm=cnorm,\
                    vmin=vmin, vmax=vmax, \
                    cmap = cmap,\
                    origin="lower",\
                    extent= [x_lim[1],\
                             x_lim[2],\
                             y_lim[1],\
                             y_lim[2]],\
                    interpolation=interpolation)

    for spine in ax.spines:
        spine.set_edgecolor(ticks_color)


    return im


"""
    get_imshow_log(ax::PyCall.PyObject, image::Array{<:Real}, 
                   x_lim::Array{<:Real}=zeros(2), y_lim::Array{<:Real}=zeros(2), 
                   vmin::Real=0.0, vmax::Real=0.0; 
                   cmap::String="viridis",
                    ticks_color::String="white",
                    interpolation::String="none")

Helper function to plot an `imshow` with logarithmic colorbar.
"""

def get_imshow_log(ax, image, \
                    x_lim=np.zeros(2), y_lim=np.zeros(2), \
                    vmin=0.0, vmax=0.0, \
                    cmap="viridis", ticks_color="white", \
                    interpolation="none"):\

    return get_imshow(ax, image, x_lim, y_lim, vmin, vmax, cmap=cmap, cnorm=matplotlib.colors.LogNorm(), 
                      ticks_color=ticks_color, interpolation=interpolation)
