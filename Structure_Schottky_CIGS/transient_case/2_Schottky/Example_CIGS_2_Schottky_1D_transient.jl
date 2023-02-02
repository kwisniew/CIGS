#=
# CIGS pn junction: stationary with trap and Schottky contacts.
Simulating stationary charge transport in a pn junction with hole traps and a Schottky boundary condition.
=#

module Example_CIGS_2_Schottky_1D_transient

using ChargeTransport
using ExtendableGrids
using PyPlot
using DelimitedFiles
using FileIO, JLD2

## function to initialize the grid for a possible extension to other p-i-n devices.
function initialize_pin_grid(refinementfactor, h_pdoping, h_pdoping_right)
    coord_pdoping_left    = collect(range(0.0, stop = h_pdoping, length = 3 * refinementfactor))
    coord_pdoping_right  = collect(range(h_pdoping, stop = (h_pdoping + h_pdoping_right), length = 3 * refinementfactor))
                                 
    coord            = glue(coord_pdoping_left, coord_pdoping_right)

    return coord
end

#In this script voltage is applied to the left side of the domain, what means that "+" means reverse voltage, and "-" means forward voltage
function main(;n = 3, voltageStep=0.5, Plotter = PyPlot, plotting = false, verbose = false, test = false, unknown_storage=:sparse)

    ################################################################################
    println("Set up grid and regions")
    ################################################################################

    ## region numbers
    regionAcceptor      = 1  # p doped region
    regions             = [regionAcceptor]
    numberOfRegions     = length(regions)

    ## boundary region numbers
    bregionAcceptorLeft     = 1
    bregionAcceptorRight    = 2
    bregions                = [bregionAcceptorRight, bregionAcceptorLeft]
    numberOfBoundaryRegions = length(bregions)

    ## grid
    refinementfactor        = 2^(n-1)
    h_pdoping          = 2.5 * μm

    coord                   = initialize_pin_grid(refinementfactor,
                                                  h_pdoping/2,
                                                  h_pdoping/2,
                                                 )

    grid                    = simplexgrid(coord)

    ## set different regions in grid, doping profiles do not intersect
    ## p doped                    
    cellmask!(grid, [0.0 * μm], [h_pdoping], regionAcceptor)    


    if plotting
        Plotter.figure()
        gridplot(grid, Plotter = Plotter, legend=:lt)
        Plotter.title("Grid")
    end

    println("*** done\n")

    ################################################################################
    println("Define physical parameters and model")
    ################################################################################

    iphin             = 1 # index electron quasi Fermi potential
    iphip             = 2 # index hole quasi Fermi potential
    numberOfCarriers  = 2 # electrons, holes and no traps
    
    Ec_CIGS           = 1.1                  *  eV
    Ev_CIGS           = 0.0                  *  eV

    Nc                = 4.351959895879690e17 / (cm^3)
    Nv                = 9.139615903601645e18 / (cm^3)
    Nt                = 5e14                 / (cm^3)
    mun_CIGS          = 100.0                * (cm^2) / (V * s)
    mup_CIGS          = 25                   * (cm^2) / (V * s)
    εr_CIGS           = 13.6                 *  1.0               
    T                 = 300.0                *  K

    An                = 4 * pi * q * mₑ * kB^2 / Planck_constant^3
    Ap                = 4 * pi * q * mₑ * kB^2 / Planck_constant^3
    vn                = An * T^2 / (q*Nc)
    vp                = Ap * T^2 / (q*Nv)
    barrier_right     = 0.7 * eV #with respect to conduction band
    barrier_left      = 0.1 * eV #with respect to conduction band

    ## recombination parameters
    Auger             = 1.0e-29  * cm^6 / s          # 1.0e-41 m^6 / s
    SRH_LifeTime      = 1.0e-3   * ns               
    Radiative         = 1.0e-10  * cm^3 / s          # 1.0e-16 m^3 / s
    G                 = 1.0e20   / (cm^3 * s)
    ## ???
    A_CIGS            = 1.0e5    / cm
    N0                = 1e17     / cm^2/s

    ## doping -- trap doping will not be set and thus automatically zero
    Na                = 1.0e15 / (cm^3)   

    ## we will impose this applied voltage on one boundary
    voltageStep   = voltageStep * V

    println("*** done\n")
    ################################################################################
    println("Define System and fill in information about model")
    ################################################################################

    ## initialize Data instance and fill in data
    data                                = Data(grid, numberOfCarriers)

    ## Possible choices: Stationary, Transient
    data.modelType                     = Transient

    ## possible choices: Boltzmann, FermiDiracOneHalfBednarczyk, FermiDiracOneHalfTeSCA FermiDiracMinusOne, Blakemore
    data.F                             .= [FermiDiracOneHalfTeSCA, FermiDiracOneHalfTeSCA]

    data.bulkRecombination              = set_bulk_recombination(;iphin = iphin, iphip = iphip,
                                                                 bulk_recomb_Auger = false,
                                                                 bulk_recomb_radiative = false,
                                                                 bulk_recomb_SRH = false)

    #enable_traps!(data)
    
    ## Possible choices: GenerationNone, GenerationUniform, GenerationBeerLambert
    # NOTE: GenerationBeerLambert doesn't work in 2D case!
    data.generationModel                = GenerationNone#GenerationBeerLambert

    ## Possible choices: OhmicContact, SchottkyContact (outer boundary) and InterfaceModelNone,
    ## InterfaceModelSurfaceReco (inner boundary).
    data.boundaryType[bregionAcceptorLeft ]    = SchottkyContact                       
    data.boundaryType[bregionAcceptorRight]    = SchottkyContact#OhmicContact   
    
    ## Choose flux discretization scheme: ScharfetterGummel, ScharfetterGummelGraded,
    ## ExcessChemicalPotential, ExcessChemicalPotentialGraded, DiffusionEnhanced, GeneralizedSG
    data.fluxApproximation              .= ExcessChemicalPotential
   
    println("*** done\n")

    ################################################################################
    println("Define Params and fill in physical parameters")
    ################################################################################

    ## physical parameters
    params                                              = Params(grid, numberOfCarriers)
    params.temperature                                  = T
    params.UT                                           = (kB * params.temperature) / q
    params.chargeNumbers[iphin]                         = -1
    params.chargeNumbers[iphip]                         =  1

    for ibreg in 1:numberOfBoundaryRegions   # boundary region data
        params.bDensityOfStates[iphin, ibreg]           = Nc
        params.bDensityOfStates[iphip, ibreg]           = Nv
        # params.bDensityOfStates[iphit, ibreg]           = Nt_low
        # params.bBandEdgeEnergy[iphit, ibreg]            = Et
    end

    params.bBandEdgeEnergy[iphin, bregionAcceptorRight]     = Ec_CIGS
    params.bBandEdgeEnergy[iphip, bregionAcceptorRight]     = Ev_CIGS
    params.bBandEdgeEnergy[iphin, bregionAcceptorLeft]      = Ec_CIGS
    params.bBandEdgeEnergy[iphip, bregionAcceptorLeft]      = Ev_CIGS

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
        params.recombinationAuger[iphin, ireg]          = Auger
        params.recombinationAuger[iphip, ireg]          = Auger

        ## generation parameters
        params.generationAbsorption[ireg]               = A_CIGS
        params.generationIncidentPhotonFlux[ireg]       = N0
        params.generationUniform[ireg]                  = G
        
    end

    ## doping -- since we do not set any doping for the traps it is automatically zero
    params.doping[iphip, regionAcceptor]             = Na    

    ## boundary doping
    params.bDoping[iphip, bregionAcceptorRight]          = Na        
    params.bDoping[iphip, bregionAcceptorLeft]           = Na   
    
    ## values for the schottky contacts
    params.SchottkyBarrier[bregionAcceptorLeft]             = barrier_left
    params.bVelocity[iphin,bregionAcceptorLeft]             = vn 
    params.bVelocity[iphip,bregionAcceptorLeft]             = vp 

    params.SchottkyBarrier[bregionAcceptorRight]                = barrier_right
    params.bVelocity[iphin,bregionAcceptorRight]                = vn 
    params.bVelocity[iphip,bregionAcceptorRight]                = vp 


    data.params                                         = params
    ctsys                                               = System(grid, data, unknown_storage=unknown_storage)

    show_params(ctsys)
    println("*** done\n")

    ################################################################################
    println("Define outerior boundary conditions and enabled layers")
    ################################################################################

    ## set ohmic contact in bregionAcceptorRight and schottky contact in bregionAcceptorLeft
    set_contact!(ctsys, bregionAcceptorRight, Δu = 0.0)
    set_contact!(ctsys, bregionAcceptorLeft,  Δu = 0.0)


    println("*** done\n")

    ################################################################################
    println("Define control parameters for Newton solver")
    ################################################################################

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

    println("*** done\n")

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


    println("*** done\n")


    if plotting 
        ## ##### set legend for plotting routines #####
        label_solution, label_density, label_energy = set_plotting_labels(data)

        ## for electrons 
        label_energy[1, iphin] = "\$E_c-q\\psi\$"       
        label_energy[2, iphin] = "\$ - q \\varphi_n\$"
        label_density[iphin]   = "n";                    
        label_solution[iphin]  = "\$ \\varphi_n\$"

        ## for holes 
        label_energy[1, iphip] = "\$E_v-q\\psi\$"       
        label_energy[2, iphip] = "\$ - q \\varphi_p\$"
        label_density[iphip]   = "p";                    
        label_solution[iphip]  = "\$ \\varphi_p\$"

        ## ##### set legend for plotting routines #####
        Plotter.figure()
        plot_energies(Plotter, ctsys, solution, "Equilibrium", label_energy)
        Plotter.figure()
        plot_densities(Plotter, ctsys, solution,"Equilibrium", label_density)
    end

    ################################################################################
    println("Stationary loop to increase bias")
    ################################################################################
    
    ## set calculationType to OutOfEquilibrium for starting with respective simulation.
    data.calculationType = OutOfEquilibrium
    biasValues           = range(0.0, stop = voltageStep, length = 11)

    ## increase slowly the bias in stationary calculations

    for Δu in biasValues

        ## Apply new voltage: set non equilibrium boundary conditions
        set_contact!(ctsys, bregionAcceptorLeft, Δu = Δu)

        if test == false
            println("bias value: Δu = $(Δu)")
        end

        # the Inf as tstep corresponds to stationary calculations.
        solve!(solution, initialGuess, ctsys, control  = control, tstep = 1e-12)


        PyPlot.clf()
        plot_solution(Plotter, ctsys, solution, "Equilibrium", label_solution)
        PyPlot.pause(0.5)

        initialGuess .= solution

    end # bias values


    ################################################################################
    println("Transient bias loop")
    ################################################################################
    
    ## set calculationType to OutOfEquilibrium for starting with respective simulation.
    data.calculationType = OutOfEquilibrium      # Rn = Rp = R, since the model type is stationary
    Voltage_step         = voltageStep            # final bias value

    ## Scan rate and time steps
    number_tsteps        = 500
    tend                 = 1e-8
    tvalues              = range(0.0, stop = tend, length = number_tsteps)
    # Δt                   = tvalues[2] - tvalues[1]

    IV         = zeros(0)     
    chargeDensities = zeros(0)
    #chargeDensities2 = zeros(0)
    

    w_device = 1.0    * cm  # width of device
    z_device = 1.0    * cm  # depth of device

    ## adjust Newton parameters
    control                = NewtonControl()
    control.verbose        = verbose
    control.tol_absolute   = 1.0e-10
    control.tol_relative   = 1.0e-10
    control.tol_round      = 1.0e-4
    control.damp_initial   = 0.5
    control.damp_growth    = 1.61
    control.max_iterations = 100
    control.max_round      = 3

    
    
    for istep = 2:number_tsteps

        t  = tvalues[istep]          # Actual time
        Δu = Voltage_step           # Applied voltage
        Δt = t - tvalues[istep-1]    # Time step size

        ## Apply new voltage: set non equilibrium boundary conditions
        set_contact!(ctsys, bregionAcceptorLeft, Δu = Δu)

        if test == false
            println("time value: t = $(t)")
        end

        solve!(solution, initialGuess, ctsys, control  = control, tstep = Δt)

        ## get I-V data
        current = get_current_val(ctsys, solution, initialGuess, Δt)
        push!(IV, w_device * z_device * current)
        push!(chargeDensities,w_device * z_device *(charge_density(ctsys,solution)[regionAcceptor]))

        initialGuess .= solution

    end # bias loop

       println("*** done\n")

    ## compute static capacitance: check this is correctly computed
    # staticCapacitance = diff(chargeDensities) ./ diff(biasValues)
    # writedlm( "staticCapacitance.csv",  staticCapacitance, ',')
    writedlm( "chargeDensities.csv"  ,  chargeDensities  , ',')
    writedlm( "biasValues.csv"       ,  biasValues       , ',')

    ## plot solution and IV curve
    if plotting 
        plot_energies(Plotter, ctsys, solution, "bias \$\\Delta u\$ = $(voltageStep), \$ t=$(tend)\$", label_energy)
        Plotter.figure()
        plot_densities(Plotter, ctsys, solution,"bias \$\\Delta u\$ = $(voltageStep), \$ t=$(tend)\$", label_density)
        Plotter.figure()
        plot_IV(Plotter, tvalues,IV, "\$ t_{end}=$(tvalues[end])s\$", plotGridpoints = true)
        Plotter.ylabel("I[A]")
        Plotter.xlabel("t[s]")
        Plotter.figure()
        plot_IV(Plotter, tvalues,chargeDensities, tvalues[end], plotGridpoints = true)
        Plotter.title("Uncompensated charge density")
        Plotter.ylabel("Charge density [C]")
        # Plotter.figure()
        # plot_IV(Plotter, tvalues,abs.(staticCapacitance), tvalues[end-1], plotGridpoints = true)
        # Plotter.title("Static capacitance in donor region")
        # Plotter.ylabel("Static capacitance [F]")
               
    end


    ################################################################################
    println("Changing time step")
    ################################################################################    
    ## Scan rate and time steps
    number_tsteps             = 500
    tend_slow                 = 0.2e-6
    tvalues_slow              = range(tvalues[end], stop = tend_slow, length = number_tsteps)

    for istep = 2:number_tsteps

        t  = tvalues_slow[istep]          # Actual time
        Δu = Voltage_step           # Applied voltage
        Δt = t - tvalues_slow[istep-1]    # Time step size

        ## Apply new voltage: set non equilibrium boundary conditions
        set_contact!(ctsys, bregionAcceptorLeft, Δu = Δu)

        if test == false
            println("time value: t = $(t)")
        end

        solve!(solution, initialGuess, ctsys, control  = control, tstep = Δt)

        ## get I-V data
        current = get_current_val(ctsys, solution, initialGuess, Δt)
        push!(IV, w_device * z_device * current)
        push!(chargeDensities,w_device * z_device *(charge_density(ctsys,solution)[regionAcceptor]))

        initialGuess .= solution

    end # bias loop

    println("*** done\n")
    FileIO.save("solution.jld2","solution",solution)
    tvalues_all = [tvalues' tvalues_slow']

    ## plot solution and IV curve
    if plotting 
        plot_energies(Plotter, ctsys, solution, "bias \$\\Delta u\$ = $(Voltage_step), \$ t=$(tend_slow)\$", label_energy)
        Plotter.figure()
        plot_densities(Plotter, ctsys, solution,"bias \$\\Delta u\$ = $(Voltage_step), \$ t=$(tend_slow)\$", label_density)
        Plotter.figure()
        plot_IV(Plotter, tvalues_all,IV, "\$ t_{end}=$(tvalues_all[end])s\$", plotGridpoints = true)
        Plotter.ylabel("total current[A]")
        Plotter.xlabel("t[s]")
        Plotter.figure()
        plot_IV(Plotter, tvalues_all,chargeDensities, tvalues_all[end], plotGridpoints = true)
        Plotter.title("Uncompensated charge density \$ t_{end}=$(tvalues_all[end])s\$")
        Plotter.ylabel("Charge density [C]")
        Plotter.xlabel("t[s]")
               
    end

    testval = solution[data.index_psi, 10]
    return testval

    println("*** done\n")

end #  main

function test()
    testval = 1.3214196490674017
    main(test = true, unknown_storage=:dense) ≈ testval && main(test = true, unknown_storage=:sparse) ≈ testval
end

println("This message should show when this module has successfully recompiled.")


end # module



