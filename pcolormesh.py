import matplotlib
import matplotlib.pyplot as plt
from numpy import linspace
import math
import numpy as np
"""
    get_pcolormesh_log( map, x_lim, y_lim, c_lim=nothing; 
                        image_cmap="plasma", cnorm=matplotlib.colors.LogNorm(),
                        set_bad::Bool=true, bad_color::String="white" )

Plot a pcolormesh of `map` with logarithmic xy-axis given by `x_lim`, `y_lim`.
If `c_lim` is not defined the colorbar is limited by minimum and maximum of the `map`.

## Keyword Arguments:
- `image_cmap="plasma"`: Name of the colormap to use.
- `cnorm=matplotlib.colors.LogNorm()`: Norm for `pcolormesh`. 
- `set_bad::Bool=true`: Whether to set invalid pixels to a different color.
- `bad_color::String="white"`: Color to use if `set_bad=true`.
"""
def get_pcolormesh_log(map, x_lim, y_lim, c_lim=None,image_cmap="plasma", cnorm=matplotlib.colors.LogNorm(), set_bad=True, bad_color="white"):

    # reconstruct corner point for pcolormesh
    X = 10.0**linspace(math.log10(x_lim[1]), math.log10(x_lim[2]), np.size(map,2))
    Y = 10.0**linspace(math.log10(y_lim[1]), math.log10(y_lim[2]), np.size(map,1))

    if set_bad == True:
        # get colormap object
        cmap = plt.get_cmap(image_cmap)
        # set invalid pixels to white
        cmap.set_bad(bad_color)


    # if no colorbar limits are set -> use minmax
    if c_lim == None:
        c_lim = [min(map), max(map)]


    # plot pcolormesh
    im1 = plt.pcolormesh(X, Y, map, \
                     cmap=image_cmap, vmin=c_lim[1], vmax=c_lim[2], \
                     norm=cnorm)

    return im1



"""
    get_pcolormesh( map, x_lim, y_lim, c_lim=nothing; 
                    image_cmap="plasma", cnorm=matplotlib.colors.LogNorm(),
                    set_bad::Bool=true, bad_color::String="white" )

Plot a pcolormesh of `map` with linear xy-axis given by `x_lim`, `y_lim`.
If `c_lim` is not defined the colorbar is limited by minimum and maximum of the `map`.

## Keyword Arguments:
- `image_cmap="plasma"`: Name of the colormap to use.
- `cnorm=matplotlib.colors.LogNorm()`: Norm for `pcolormesh`. 
- `set_bad::Bool=true`: Whether to set invalid pixels to a different color.
- `bad_color::String="white"`: Color to use if `set_bad=true`.
"""
def get_pcolormesh(map, x_lim, y_lim, c_lim=None,image_cmap="plasma", cnorm=matplotlib.colors.LogNorm(),set_bad=True, bad_color="white"):

    # reconstruct corner point for pcolormesh
    X = linspace(x_lim[1], x_lim[2], np.size(map,2))
    Y = linspace(y_lim[1], y_lim[2], np.size(map,1))

    if set_bad == True:
        # get colormap object
        cmap = plt.get_cmap(image_cmap)
        # set invalid pixels to white
        cmap.set_bad(bad_color)


    # if no colorbar limits are set -> use minmax
    if c_lim == None:
        c_lim = [min(map), max(map)]


    # plot pcolormesh
    im1 = plt.pcolormesh(X, Y, map, \
                     cmap=image_cmap, vmin=c_lim[1], vmax=c_lim[2], \
                     norm=cnorm)

    return im1
