import matplotlib.pyplot as PyPlot
import numpy as np

"""
    linestyle(sequence::String="";
                   dash::Real=3.7, dot::Real=1.0, space::Real=3.7,
                   buffer::Bool=false, offset::Int=0)

Get custom `IDL`-like linestyles. `sequence` has to be an even number of entries.

## Keyword Arguments
- `dash::Real=3.7`: Size of dash.
- `dot::Real=1.0`: Size of dot.
- `space::Real=3.7`: Size of space.
- `buffer::Bool=false`: Whether to add a buffer to the end of the sequence.
- `offset::Int=0`: No clue.
"""
def linestyle(sequence = "",\
    dash = 3.7, dot= 1.0, space = 3.7,\
    buffer = False, offset = 0):


    if sequence == "-" :
        return sequence

    if (sequence.count('.') > 0) or ( sequence.count(':')> 0):
        dash = 6.4

    if (sequence == "" or sequence == "-" or sequence == "_"):
        return (offset, [dash])

    elif sequence == ":" :
        return (offset, [dot, space, dot, space, dot, space])
    else:
        reftype = Dict('-' => [dash, space],\
            '_' => [2dash, space],\
            '.' => [dot, space],\
            ':' => [dot, space, dot, space, dot, space])

        onoffseq = Vector{Float64}(undef, 0)

        for c ∈ sequence :
            onoffseq = [onoffseq; reftype[c]]


        if buffer==True :
            onoffseq[-1] = 1.0

        return (offset, onoffseq)




"""
    ls(sequence::String="";
       dash::Real=3.7, dot::Real=1.0, space::Real=3.7,
       buffer::Bool=false, offset::Int=0)

Get custom `IDL`-like linestyles. `sequence` has to be an even number of entries.
See also [`linestyle`](@ref)
"""

def  ls(sequence = "",\
    dash = 3.7, dot = 1.0, space = 3.7,\
    buffer = False, offset= 0):

    return linestyle(sequence, dash = dash, dot = dot, space = space, buffer = buffer, offset = offset)
