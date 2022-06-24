import matplotlib.pyplot as PyPlot
import mpl_toolkits.axes_grid1.axes_divider as axes_divider
import colorbar


"""
    slice(rng::UnitRange)
Helper function to convert `Julia` `UnitRange` to `Python` `slice`.
"""
def slice(rng:UnitRange):
    pycall(pybuiltin("slice"), PyObject, first(rng), last(rng))


"""
    get_gs(gs::PyObject, row, col)
Helper function to get the correct `gridspec` panel for a given (range of) `row` and `col`.
"""
def get_gs(gs, row, col):

    if (typeof(row) <: UnitRange  typeof(col) <: UnitRange ):
        return get(gs, (slice(row), slice(col)))
    elif typeof(row) <: UnitRange:
        return get(gs, (slice(row), col))
    elif typeof(col) <: UnitRange:
        return get(gs, (row, slice(col)))
    else:
        return get(gs, (row, col))