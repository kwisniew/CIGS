#File for grids creation in external files


module grids

using VoronoiFVM
using ChargeTransport
using ExtendableGrids
using GridVisualize
using PyPlot
using DelimitedFiles
using SimplexGridFactory
using Triangulate


# in onther files => using grids (and when calling functions grids.func)

function schottky()
    ## region numbers
    regionCIGS      = 1  # CIGS
    regionGrainBoundary  = 2  # grain boundary region
    #regionCIGS      = 3  # CIGS
    regions              = [regionCIGS, regionGrainBoundary]
    numberOfRegions      = length(regions)

    ## boundary region numbers
    bregionCIGSLeft  = 1
    bregionCIGSRight = 2
    bregionGB        = 3
    bregionNoFlux    = 4   

    ## grid
    total_width = 2 * μm
    width_CIGS_Left_down  = 1.0 * μm
    width_CIGS_Left_up    = 1.5 * μm  
    width_grain_boundary  = 0.1 * μm
    width_CIGS_Right_down = total_width - width_CIGS_Left_down - width_grain_boundary
    width_CIGS_Right_up   = total_width - width_CIGS_Left_up - width_grain_boundary 
    height                = 1.0 * μm

    b                = SimplexGridBuilder(Generator=Triangulate)

    ## specify boundary nodes
    length_0         = point!(b, 0.0, 0.0)
    length_CIGS_left = point!(b, width_CIGS_Left_down, 0.0)
    length_gb        = point!(b, width_CIGS_Left_down + width_grain_boundary, 0.0)
    length_L         = point!(b, width_CIGS_Left_down + width_grain_boundary + width_CIGS_Right_down, 0.0)

    height_0         = point!(b, 0.0, height)
    height_CIGS_left = point!(b, width_CIGS_Left_up, height)
    height_gb        = point!(b, width_CIGS_Left_up + width_grain_boundary, height)
    height_L         = point!(b, width_CIGS_Left_up + width_grain_boundary + width_CIGS_Right_up, height)

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
    facetregion!(b, bregionGB)
    facet!(b, length_CIGS_left, height_CIGS_left)
    facetregion!(b, bregionGB)
    facet!(b, length_gb, height_gb)

    ## cell regions
    cellregion!(b, regionCIGS)
    regionpoint!(b, width_CIGS_Left_down/2, height/2)
    cellregion!(b,regionGrainBoundary)
    regionpoint!(b, width_CIGS_Left_down + width_grain_boundary/2, height/100)
    cellregion!(b,regionCIGS)
    regionpoint!(b, width_CIGS_Left_down + width_grain_boundary + width_CIGS_Right_down/2, height/2)

    options!(b,maxvolume=(0.1*μm)^2)

    # grid           = simplexgrid(b)

    # GridVisualize.gridplot(grid, Plotter= PyPlot, resolution=(600,400),linewidth=0.5, legend=:lt)
    # PyPlot.title("Grid")
    
    builderplot(b,Plotter=PyPlot,resolution=(750,700), legend=:lt)
    
end #Schottky

######################################################################################################

function pn()

     ## region numbers
 regionZnO            = 1  # ZnO
 regionCIGS           = 2  # CIGS
 regionGrainBoundary  = 3  # grain boundary region
 regions              = [regionZnO, regionCIGS, regionGrainBoundary]
 numberOfRegions      = length(regions)

 ## boundary region numbers
 bregionZnO       = 1
 bregionCIGSRight = 2
 bregionInner     = 3
 bregionNoFlux    = 4   

 ## grid
 total_width = 2.5 * μm
 width_ZnO  = 0.5 * μm 
 width_CIGS_Left_down  = 1.0 * μm
 width_CIGS_Left_up    = 1.0 * μm #bigger or smaller will bend the grain boundary e.g. 1.5 * μm  
 width_grain_boundary  = 0.1 * μm
 width_CIGS_Right_down = total_width - width_CIGS_Left_down - width_grain_boundary
 width_CIGS_Right_up   = total_width - width_CIGS_Left_up - width_grain_boundary 
 height                = 1.0 * μm

b                = SimplexGridBuilder(Generator=Triangulate)

## specify boundary nodes
length_0         = point!(b, 0.0, 0.0)
length_ZnO       = point!(b, width_ZnO, 0.0)
length_CIGS_left = point!(b, width_ZnO + width_CIGS_Left_down, 0.0)
length_gb        = point!(b, width_ZnO + width_CIGS_Left_down + width_grain_boundary, 0.0)
length_L         = point!(b, width_ZnO + width_CIGS_Left_down + width_grain_boundary + width_CIGS_Right_down, 0.0)

height_0         = point!(b, 0.0, height)
height_ZnO       = point!(b, width_ZnO, height)
height_CIGS_left = point!(b, width_ZnO + width_CIGS_Left_up, height)
height_gb        = point!(b, width_ZnO + width_CIGS_Left_up + width_grain_boundary, height)
height_L         = point!(b, width_ZnO + width_CIGS_Left_up + width_grain_boundary + width_CIGS_Right_up, height)

## specify boundary regions
## metal interface
facetregion!(b, bregionZnO)
facet!(b, length_0, height_0)
facetregion!(b, bregionCIGSRight)
facet!(b, length_L, height_L)

## no flux
facetregion!(b, bregionNoFlux)
facet!(b, length_0, length_L)
facetregion!(b, bregionNoFlux)
facet!(b, height_0, height_L)

#  ## inner interface
facetregion!(b, bregionInner)
facet!(b, length_ZnO, height_ZnO)
facetregion!(b, bregionInner)
facet!(b, length_CIGS_left, height_CIGS_left)
facetregion!(b, bregionInner)
facet!(b, length_gb, height_gb)


## cell regions
cellregion!( b, regionZnO)
regionpoint!(b, width_ZnO/2, height/2)
cellregion!( b, regionCIGS)
regionpoint!(b, width_ZnO+width_CIGS_Left_down/2, height/2)
cellregion!( b, regionGrainBoundary)
regionpoint!(b, width_ZnO+width_CIGS_Left_down + width_grain_boundary/2, height/100)
cellregion!( b, regionCIGS)
regionpoint!(b, width_ZnO+width_CIGS_Left_down + width_grain_boundary + width_CIGS_Right_down/2, height/2)

options!(b,maxvolume=(0.1*μm)^2)

grid           = simplexgrid(b)

GridVisualize.gridplot(grid, Plotter= PyPlot, resolution=(600,400),linewidth=0.5, legend=:lt)
PyPlot.title("Grid")
# builderplot(b,Plotter=PyPlot,resolution=(750,700), legend=:lt)

end #PN

######################################################################################################

function etc()

end # 


end #module