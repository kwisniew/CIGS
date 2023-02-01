#=
# CIGS pn junction: stationary with trap and Schottky contacts.
Simulating stationary charge transport in a pn junction with hole traps and a Schottky boundary condition.
=#

#set LC_NUMERIC="C"

module Example_CIGS_Schottky_with_traps_2D_unstr_grid

using ChargeTransport
using ExtendableGrids
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

function main(;n = 3, voltageMin=-0.5, voltageMax=0.1, Plotter = PyPlot, plotting = false, unknown_storage=:sparse, verbose = false)
 
## region numbers
 regionCIGS             = 1  # CIGS
 regionGrainBoundary    = 2  # grain boundary region
 #regionCIGS            = 3  # CIGS
#  regionInterface        = 3
 regions              = [regionCIGS, regionGrainBoundary]
 numberOfRegions      = length(regions)



 ## boundary region numbers
 regionCIGSLeft  = 1
 regionCIGSRight = 2
 bregionGB       = 3
 bregionNoFlux   = 4   
 #bregions                = [regionCIGSRight, regionCIGSLeft, bregionInterface, bregionNoFlux]
 numberOfBoundaryRegions = length(regions)   

 ## grid
 total_width = 2 * μm
 width_CIGS_Left_down  = 1.0 * μm #Adam: to jest pochylenie 
 width_CIGS_Left_up    = 1.0 * μm  
 width_grain_boundary  = 0.01 * μm #Adam: To jest szerokość 
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
facetregion!(b, regionCIGSLeft)
facet!(b, length_0, height_0)
facetregion!(b, regionCIGSRight)
facet!(b, length_L, height_L)

## no flux
facetregion!(b, bregionNoFlux)
facet!(b, length_0, length_L)
facetregion!(b, bregionNoFlux)
facet!(b, height_0, height_L)

##   inner interface
facetregion!(b, bregionGB)
facet!(b, length_CIGS_left, height_CIGS_left)
facetregion!(b, bregionGB)
facet!(b, length_gb, height_gb)

## cell regions
cellregion!(b, regionCIGS)
regionpoint!(b, width_CIGS_Left_down/2, height/2)
cellregion!(b,regionGrainBoundary)
# maxvolume!(b,(0.01*μm)^2)
regionpoint!(b, width_CIGS_Left_down + width_grain_boundary/2, height/100)
cellregion!(b,regionCIGS)
regionpoint!(b, width_CIGS_Left_down + width_grain_boundary + width_CIGS_Right_down/2, height/2)

options!(b,maxvolume=(0.1*μm)^2)

grid  = simplexgrid(b)

# GridVisualize.gridplot(grid, Plotter= PyPlot, resolution=(600,400),linewidth=0.5, legend=:lt)
# PyPlot.title("Grid")
#builderplot(b,Plotter=PyPlot,resolution=(750,700), legend=:lt)

if plotting
    
    Plotter.figure()
    gridplot(grid, Plotter = Plotter, legend=:lt)
    Plotter.title("Grid")

end


############################################################################################################
#physical parameters and model

    iphin             = 1 # index electron quasi Fermi potential
    iphip             = 2 # index hole quasi Fermi potential
    numberOfCarriers  = 2 # electrons, holes and no traps
    
    Ec_CIGS           = 1.1                  *  eV
    Ev_CIGS           = 0.0                  *  eV
    Et                = Ev_CIGS + 0.5        *  eV

    Nc                = 4.351959895879690e17 / (cm^3)
    Nv                = 9.139615903601645e18 / (cm^3)
    Nt                = 1e18                 / (cm^3)               
    mun_CIGS          = 100.0                * (cm^2) / (V * s)
    mup_CIGS          = 25                   * (cm^2) / (V * s)
    mut               = 0                    * (cm^2) / (V * s)  # no flux for traps
    εr_CIGS           = 13.6                 *  1.0            
    T                 = 300.0                *  K

    An                = 4 * pi * q * mₑ * kB^2 / Planck_constant^3
    Ap                = 4 * pi * q * mₑ * kB^2 / Planck_constant^3
    vn                = An * T^2 / (q*Nc)
    vp                = Ap * T^2 / (q*Nv)
    # barrier_right     = Ev_CIGS + 0.4 * eV
    barrier_left      = 0.0 * eV

    ## recombination parameters
    ni_CIGS           = sqrt(Nc * Nv) * exp(-(Ec_CIGS - Ev_CIGS) / (2 * kB * T)) # intrinsic concentration
    n0_CIGS           = Nc * Boltzmann( (Et-Ec_CIGS) / (kB*T) )                  # Boltzmann equilibrium concentration
    p0_CIGS           = ni_CIGS^2 / n0_CIGS                                      # Boltzmann equilibrium concentration

    
    Auger             = 1.0e-29  * cm^6 / s          # 1.0e-41 m^6 / s
    SRH_LifeTime      = 1.0e-3   * ns               
    Radiative         = 1.0e-10  * cm^3 / s          # 1.0e-16 m^3 / s
    G                 = 1.0e20   / (cm^3 * s)
    ## ???
    A_CIGS            = 1.0e5    / cm
    N0                = 1e17     / cm^2/s

    ## doping -- trap doping will not be set and thus automatically zero
    Na                = 5.0e15 / (cm^3)   

    ## we will impose this applied voltage on one boundary
    voltageMin   = voltageMin * V 
    voltageMax   = voltageMax * V

    ################################################################################

    ## initialize Data instance and fill in data
    data                                = Data(grid, numberOfCarriers)

    ## possible choices: model_stationary, model_transient
    data.modelType                      = Stationary

    ## possible choices: Boltzmann, FermiDiracOneHalfBednarczyk, FermiDiracOneHalfTeSCA FermiDiracMinusOne, Blakemore
    data.F                             .= [FermiDiracOneHalfTeSCA, FermiDiracOneHalfTeSCA]

    data.bulkRecombination              = set_bulk_recombination(;iphin = iphin, iphip = iphip,
                                                                 bulk_recomb_Auger = true,
                                                                 bulk_recomb_radiative = true,
                                                                 bulk_recomb_SRH = true)

   
    zt = 1
    Nt_vector = similar(data.params.densityOfStates[zt,:])
    Nt_vector[regionGrainBoundary] = Nt
    #enable_traps!(data = data, traps = iphip, regions = regions)
    # enable_traps!(data,z=zt,Nt=Nt_vector)   #traps
    NT = [0, Nt]
    add_trap_density!(data=data, zt = zt, Nt = NT)
    
    ## Possible choices: GenerationNone, GenerationUniform, GenerationBeerLambert
    data.generationModel                = GenerationNone#GenerationBeerLambert

    ## Possible choices: OhmicContact, SchottkyContact (outer boundary) and InterfaceModelNone,
    ## InterfaceModelSurfaceReco (inner boundary).
    #data.boundary_type[bregionAcceptorLeft ]    = ohmic_contact#schottky_contact                       
    #data.boundary_type[bregionAcceptorRight]    = ohmic_contact#schottky_contact   
    data.boundaryType[regionCIGSLeft ]    = SchottkyContact                       
    data.boundaryType[regionCIGSRight]    = OhmicContact
    # data.boundaryType[regionInterface]    = InterfaceModelNone    
    
    ## Choose flux discretization scheme: ScharfetterGummel, ScharfetterGummelGraded,
    ## ExcessChemicalPotential, ExcessChemicalPotentialGraded, DiffusionEnhanced, GeneralizedSG
    data.fluxApproximation              .= ExcessChemicalPotential
   

############################################################################################################
    #Define System and fill in information about model

     ## physical parameters
     params                                              = Params(grid, numberOfCarriers)
     params.temperature                                  = T
     params.UT                                           = (kB * params.temperature) / q
     params.chargeNumbers[iphin]                         =  -1
     params.chargeNumbers[iphip]                         =  1
     # params.chargeNumbers[iphit]                         =  1  # +1: hole trap is used
 
     for ibreg in 1:numberOfBoundaryRegions   # boundary region data
         params.bDensityOfStates[iphin, ibreg]           = Nc
         params.bDensityOfStates[iphip, ibreg]           = Nv
         # params.bDensityOfStates[iphit, ibreg]           = Nt_low
         # params.bBandEdgeEnergy[iphit, ibreg]            = Et
     end
 
     params.bBandEdgeEnergy[iphin, regionCIGSRight]     = Ec_CIGS
     params.bBandEdgeEnergy[iphip, regionCIGSRight]     = Ev_CIGS
     params.bBandEdgeEnergy[iphin, regionCIGSLeft]      = Ec_CIGS
     params.bBandEdgeEnergy[iphip, regionCIGSLeft]      = Ev_CIGS
 
     for ireg in 1:numberOfRegions           # interior region data
 
         params.dielectricConstant[ireg]                 = εr_CIGS*ε0       
 
         ## effective DOS, band-edge energy and mobilities
         params.densityOfStates[iphin, ireg]             = Nc
         params.densityOfStates[iphip, ireg]             = Nv
         params.bandEdgeEnergy[iphin, ireg]              = Ec_CIGS
         params.bandEdgeEnergy[iphip, ireg]              = Ev_CIGS
         params.mobility[iphin, ireg]                    = mun_CIGS
         params.mobility[iphip, ireg]                    = mup_CIGS
 
         ## recombination parameters
         params.recombinationRadiative[ireg]             = Radiative
         params.recombinationSRHLifetime[iphin, ireg]    = SRH_LifeTime
         params.recombinationSRHLifetime[iphip, ireg]    = SRH_LifeTime
         params.recombinationSRHTrapDensity[iphin, ireg] = n0_CIGS
         params.recombinationSRHTrapDensity[iphip, ireg] = p0_CIGS
         params.recombinationAuger[iphin, ireg]          = Auger
         params.recombinationAuger[iphip, ireg]          = Auger
 
         ## generation parameters
         params.generationAbsorption[ireg]               = A_CIGS
         params.generationIncidentPhotonFlux[ireg]       = N0
         params.generationUniform[ireg]                  = G
         
     end

############################################################################################################

    ## doping -- since we do not set any doping for the traps it is automatically zero
    params.doping[iphip, regionCIGSLeft]             = Na        
    params.doping[iphip, regionCIGSRight]            = Na        

    ## boundary doping
    params.bDoping[iphip, regionCIGSRight]          = Na        
    params.bDoping[iphip, regionCIGSLeft]           = Na   
    
    ## values for the schottky contacts
    params.SchottkyBarrier[regionCIGSLeft]             = barrier_left
    params.bVelocity[iphin,regionCIGSLeft]             = vn 
    params.bVelocity[iphip,regionCIGSLeft]             = vp 

    #params.SchottkyBarrier[bregionAcceptorRight]                = barrier_right
    #params.bVelocity[iphin,bregionAcceptorRight]                = vn 
    #params.bVelocity[iphip,bregionAcceptorRight]                = vp 

    data.params                                         = params
    ctsys                                               = System(grid, data, unknown_storage=unknown_storage)

    show_params(ctsys)
    println("*** done\n")

    ################################################################################
    println("Define outerior boundary conditions and enabled layers")
    ################################################################################

    ## set ohmic contact in bregionAcceptorRight and schottky contact in bregionAcceptorLeft
    #set_schottky_contact!(ctsys, bregionAcceptorRight, appliedVoltage = 0.0)
    #set_schottky_contact!(ctsys, bregionAcceptorLeft , appliedVoltage = 0.0)
    #set_ohmic_contact!(ctsys, bregionAcceptorRight   , 0.0)
    #set_ohmic_contact!(ctsys, bregionAcceptorLeft, 0.0)
    set_contact!(ctsys, regionCIGSRight, Δu = 0.0)
    set_contact!(ctsys, regionCIGSLeft,  Δu = 0.0)


    
    ################################################################################
    println("Define control parameters for Newton solver")

    control                   = NewtonControl()
    control.verbose           = verbose
    control.damp_initial      = 0.5
    control.damp_growth       = 1.21    #>= 1
    control.max_iterations    = 250
    control.tol_absolute      = 1.0e-14
    control.tol_relative      = 1.0e-14
    control.handle_exceptions = true
    control.tol_round         = 1.0e-8
    control.max_round         = 5

    ################################################################################
    println("Compute solution in thermodynamic equilibrium")
    ################################################################################

    ## initialize solution and starting vectors
    initialGuess          = unknowns(ctsys)
    equilibriumSolution   = unknowns(ctsys)
    solution              = unknowns(ctsys)

    #data.calculationType = inEquilibrium 

    ## solve thermodynamic equilibrium and update initial guess
    solution              = equilibrium_solve!(ctsys, control = control, nonlinear_steps = 20)
    equilibriumSolution  .= solution  
    initialGuess         .= solution

    ################################################################################


    if plotting 

        ipsi = data.index_psi
        X = grid[Coordinates][1,:]
        Y = grid[Coordinates][2,:]

        #Plot energies
        Plotter.figure()
        #TO DO: make the labels visible 
        Plotter.surf(X[:], Y[:], Ev_CIGS/q .- solution[ipsi, :], label="Ev")
        Plotter.surf(X[:], Y[:], Ec_CIGS/q .- solution[ipsi, :], label="Ec")
        Plotter.surf(X[:], Y[:], -q*solution[iphin, :], label="\$ \\varphi_n \$")
        Plotter.surf(X[:], Y[:], -q*solution[iphip, :], label="\$ \\varphi_p \$")
        
        Plotter.title("Band Edge Energies and qFermi levels in Equilibrium")
        Plotter.xlabel("length [m]")
        Plotter.ylabel("height [m]")
        Plotter.zlabel("Energy [eV]")
        Plotter.tight_layout()

        # Plot densities
        Plotter.figure()
        #TO DO: make the labels visible 
        Plotter.surf(X[:], Y[:], log10.(1.0e-6 .*Nc.*exp.((-(Ec_CIGS/q .- solution[ipsi, :])-solution[iphin, :])/(params.UT)))) #electron density
        Plotter.surf(X[:], Y[:], log10.(1.0e-6 .*Nv.*exp.(( (Ev_CIGS/q .- solution[ipsi, :])+solution[iphip, :])/(params.UT)))) #hole density
        
        Plotter.title("electrons and holes densities")
        Plotter.xlabel("length [m]")
        Plotter.ylabel("height [m]")
        Plotter.zlabel("Densities [log(cm-3)]")
        Plotter.tight_layout()


    end

    ################################################################################
    println("Stationary bias loop")
    ################################################################################
    
    ## set calculationType to OutOfEquilibrium for starting with respective simulation.
    data.calculationType = OutOfEquilibrium      # Rn = Rp = R, since the model type is stationary

    IV         = zeros(0)   
    biasSteps  = 101
    biasValues = collect(range(voltageMin, stop = voltageMax, length = biasSteps))
    if(!(0.0 in biasValues))
        append!(biasValues, 0.0)
        sort!(biasValues)
    end

    chargeDensities = zeros(0)
    

    w_device = 1.0    * cm  # width of device
    z_device = 1.0#    * cm  # depth of device

    ## adjust Newton parameters
    control.tol_absolute      = 1.0e-10
    control.tol_relative      = 1.0e-10
    control.tol_round         = 1.0e-7
    control.damp_initial      = 0.5
    control.damp_growth       = 1.2
    control.max_iterations    = 30
    control.max_round         = 3

    
    indexOfZero = indexin(0.0, biasValues)[1]
    i = indexOfZero

    #for i in eachindex(biasValues)
    while (i>0)
        Δu = biasValues[i] # bias

        ## set non equilibrium boundary condition
        #set_schottky_contact!(ctsys, bregionAcceptorLeft, appliedVoltage = Δu)

        set_contact!(ctsys, regionCIGSRight, Δu = Δu)

        ## increase generation rate with bias
        ctsys.data.λ2 = 10.0^(-biasSteps + i)
        #println("bias: Δu = $(Δu)")
     
        ## solve time step problems with timestep Δt
        solve!(solution, initialGuess, ctsys, control  = control, tstep = Inf)

        ## save IV data
        current = get_current_val(ctsys, solution)
        push!(IV, w_device * z_device * current)

        ## store CHARGE DENSITY in CIGS
        #push!(chargeDensities,chargeDensity(ctsys,solution)[regionAcceptorLeft])
        push!(chargeDensities,w_device * z_device *(charge_density(ctsys,solution)[regionCIGSLeft]+charge_density(ctsys,solution)[regionCIGSRight]))

        #set_contact!(ctsys, bregionAcceptorRight, Δu = jezuchryste)
        #push!(chargeDensities2,w_device * z_device *(charge_density(ctsys,solution)[regionAcceptorLeft]+charge_density(ctsys,solution)[regionAcceptorRight]))

        initialGuess .= solution
        i=i-1
    end # bias loop 1

    reverse!(IV)
    reverse!(chargeDensities)

    #plot energies and qFermi levels for voltageMin
    if plotting
        ipsi = data.index_psi
        X = grid[Coordinates][1,:]
        Y = grid[Coordinates][2,:]

        Plotter.figure()
        #TO DO: make the labels visible 
        Plotter.surf(X[:], Y[:], Ev_CIGS/q .- solution[ipsi, :], label="Ev")
        Plotter.surf(X[:], Y[:], Ec_CIGS/q .- solution[ipsi, :], label="Ec")
        Plotter.surf(X[:], Y[:], -solution[iphin, :], label="\$ \\varphi_n \$")
        Plotter.surf(X[:], Y[:], -solution[iphip, :], label="\$ \\varphi_p \$")
        
        Plotter.title("Band Edge Energies and qFermi levels at \$ $(voltageMin) \$V")
        Plotter.xlabel("length [m]")
        Plotter.ylabel("height [m]")
        Plotter.zlabel("Energy [eV]")
        Plotter.tight_layout() 
               
    end
    
    i = indexOfZero+1
    initialGuess .= equilibriumSolution

    while (i<=length(biasValues))
        Δu = biasValues[i] # bias

        ## set non equilibrium boundary condition
        #set_schottky_contact!(ctsys, bregionAcceptorLeft, appliedVoltage = Δu)
        set_contact!(ctsys, regionCIGSRight, Δu = Δu)
    
        ## solve time step problems with timestep Δt
        solve!(solution, initialGuess, ctsys, control  = control, tstep = Inf)

        ## save IV data
        current = get_current_val(ctsys, solution)
        push!(IV, w_device * z_device * current)

        ## store uncompensated CHARGE DENSITY in CIGS
        push!(chargeDensities,w_device * z_device *(charge_density(ctsys,solution)[regionCIGSLeft]+charge_density(ctsys,solution)[regionCIGSRight]))

        initialGuess .= solution
        i=i+1
    end # bias loop 2

   

    println("*** done\n")

    ## compute static capacitance: check this is correctly computed

    staticCapacitance = diff(chargeDensities) ./ diff(biasValues)
    writedlm( "staticCapacitance.csv",  staticCapacitance, ',')
    writedlm( "chargeDensities.csv"  ,  chargeDensities  , ',')
    writedlm( "biasValues.csv"       ,  biasValues       , ',')

    ## plot solution at voltageMax, IV, QV and CV curves
    if plotting

        ipsi = data.index_psi
        X = grid[Coordinates][1,:]
        Y = grid[Coordinates][2,:]

        #Plot energies
        Plotter.figure()
        #TO DO: make the labels visible 
        Plotter.surf(X[:], Y[:], Ev_CIGS/q .- solution[ipsi, :], label="Ev")
        Plotter.surf(X[:], Y[:], Ec_CIGS/q .- solution[ipsi, :], label="Ec")
        Plotter.surf(X[:], Y[:], -solution[iphin, :], label="\$ \\varphi_n \$")
        Plotter.surf(X[:], Y[:], -solution[iphip, :], label="\$ \\varphi_p \$")
        
        Plotter.title("Band Edge Energies and qFermi levels at \$ $(voltageMax) \$V")
        Plotter.xlabel("length [m]")
        Plotter.ylabel("height [m]")
        Plotter.zlabel("Energy [eV]")
        Plotter.tight_layout()

        #Plot densities
        Plotter.figure()
        #TO DO: make the labels visible 
        Plotter.surf(X[:], Y[:], log10.(1.0e-6 .*Nc.*exp.((-(Ec_CIGS/q .- solution[ipsi, :])-solution[iphin, :])/(params.UT)))) # electron density
        Plotter.surf(X[:], Y[:], log10.(1.0e-6 .*Nv.*exp.(( (Ev_CIGS/q .- solution[ipsi, :])+solution[iphip, :])/(params.UT)))) #hole density
        
        Plotter.title("electrons and holes densities at \$ $(voltageMax) \$V")
        Plotter.xlabel("length [m]")
        Plotter.ylabel("height [m]")
        Plotter.zlabel("Densities [log(cm-3)]")
        Plotter.tight_layout()

        #Plot IV, QV and CV
        Plotter.figure()
        plot_IV(Plotter, biasValues,IV, biasValues[end], plotGridpoints = true)
        Plotter.figure()
        plot_IV(Plotter, biasValues,chargeDensities, biasValues[end], plotGridpoints = true)
        Plotter.title("Charge density in donor region")
        Plotter.ylabel("Charge density [C]")
        Plotter.figure()
        plot_IV(Plotter, biasValues,abs.(staticCapacitance), biasValues[end-1], plotGridpoints = true)
        Plotter.title("Static capacitance in donor region")
        Plotter.ylabel("Static capacitance [F]")
               
    end
 
    testval = sum(filter(!isnan, solution))/length(solution)
    return testval

    println("*** done\n")

    

end #main 

############################################################################################################
end #module