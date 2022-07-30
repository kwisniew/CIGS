#=
# CIGS pn junction: stationary with trap and Schottky contacts.
Simulating stationary charge transport in a pn junction with hole traps and a Schottky boundary condition.
=#


module Example_CIGS_Schottky_with_traps_2D_unstr_grid

using VoronoiFVM
using ChargeTransport
using ExtendableGrids
using GridVisualize
using PyPlot
using DelimitedFiles
using SimplexGridFactory
using Triangulate

#############################################################################################################
## function to initialize the grid for a possible extension to other p-i-n devices.
function initialize_pin_grid(refinementfactor, h_pdoping_left, h_pdoping_right)
    coord_pdoping_left    = collect(range(0.0, stop = h_pdoping_left, length = 3 * refinementfactor))
    coord_pdoping_right  = collect(range(h_pdoping_left, stop = (h_pdoping_left + h_pdoping_right), length = 3 * refinementfactor))
                                 
    coord            = glue(coord_pdoping_left, coord_pdoping_right)

    return coord
end
#############################################################################################################


end #module