# See https://wias-berlin.github.io/PDELib.jl/dev/intro/#Grid-generation-and-visualization for more details and examples

using ChargeTransport
using ExtendableGrids
using GridVisualize
using SimplexGridFactory
using Triangulate
using PyPlot
 
 ## region numbers
 regionCIGSLeft       = 1  # CIGS
 regionCIGSRight      = 2  # grain boundary region
 regions              = [regionCIGSLeft, regionCIGSRight]
 numberOfRegions      = length(regions)

 ## boundary region numbers
 bregionCIGSLeft  = 1
 bregionCIGSRight = 2
 bregionInterface = 3
 bregionNoFlux    = 4
 bregions                = [bregionCIGSRight, bregionCIGSLeft, bregionNoFlux]
 numberOfBoundaryRegions = length(bregions)   

 ## grid
 total_width = 2 * μm
 width_CIGS_Left_down  = 1.0 * μm
 width_CIGS_Left_up    = 1.0 * μm  
 width_CIGS_Right_down = total_width - width_CIGS_Left_down 
 width_CIGS_Right_up   = total_width - width_CIGS_Left_up 
 height                = 1.0 * μm

b                = SimplexGridBuilder(Generator=Triangulate)

## specify boundary nodes
length_0         = point!(b, 0.0, 0.0)
length_CIGS_left = point!(b, width_CIGS_Left_down, 0.0)
length_L         = point!(b, width_CIGS_Left_down + width_CIGS_Right_down, 0.0)

height_0         = point!(b, 0.0, height)
height_CIGS_left = point!(b, width_CIGS_Left_up, height)
height_L         = point!(b, width_CIGS_Left_up + width_CIGS_Right_up, height)

## specify boundary regions
## metal interface
facetregion!(b, bregionCIGSLeft)
facet!(b, length_0, height_0)
facetregion!(b, bregionCIGSRight)
facet!(b, length_L, height_L)

## no flux
facetregion!(b, bregionNoFlux)
facet!(b, length_0, length_L)
facetregion!(b, bregionNoFlux)
facet!(b, height_0, height_L)

#  ## inner interface
facetregion!(b, bregionInterface)
facet!(b, length_CIGS_left, height_CIGS_left)

## cell regions
cellregion!(b, regionCIGSLeft)
regionpoint!(b, width_CIGS_Left_down/2, height/2)
cellregion!(b,regionCIGSRight)
regionpoint!(b, width_CIGS_Left_down + width_CIGS_Right_down/2, height/2)

options!(b,maxvolume=(0.1*μm)^2)

#  grid           = simplexgrid(b)

# GridVisualize.gridplot(grid, Plotter= PyPlot, resolution=(600,400),linewidth=0.5, legend=:lt)
# PyPlot.title("Grid")
builderplot(b,Plotter=PyPlot,resolution=(750,700))