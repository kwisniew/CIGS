#=
# CIGS pn junction: stationary with trap and Schottky contacts.
Simulating stationary charge transport in a pn junction with hole traps and a Schottky boundary condition.
=#

module Example_CIGS_2_Schottky_1D

using ChargeTransport
using ExtendableGrids
using PyPlot
using DelimitedFiles
using FileIO, JLD2

## function to initialize the grid for a possible extension to other p-i-n devices.
function initialize_pin_grid(refinementfactor, h_pdoping_left, h_pdoping_right)
    coord_pdoping_left    = collect(range(0.0, stop = h_pdoping_left, length = 3 * refinementfactor))
    coord_pdoping_right  = collect(range(h_pdoping_left, stop = (h_pdoping_left + h_pdoping_right), length = 3 * refinementfactor))
                                 
    coord            = glue(coord_pdoping_left, coord_pdoping_right)

    return coord
end

function main(;n = 3, voltageMin=0.5, voltageMax=-0.1, Plotter = PyPlot, plotting = false, verbose = false, test = false, unknown_storage=:sparse)

    ################################################################################
    println("Set up grid and regions")
    ################################################################################

    ## region numbers
    regionAcceptor  = 1                           # n doped region
    # regionAcceptor = 2                           # p doped region
    regions             = [regionAcceptor]#, regionAcceptor]
    numberOfRegions     = length(regions)

    ## boundary region numbers
    bregionAcceptorLeft     = 1
    bregionAcceptorRight    = 2
    bregions                = [bregionAcceptorRight, bregionAcceptorLeft]
    numberOfBoundaryRegions = length(bregions)

    ## grid
    refinementfactor        = 2^(n-1)
    #h_pdoping_left               = 0.5 * μm
    h_pdoping_left          = 1.25 * μm

    coord                   = initialize_pin_grid(refinementfactor,
                                                  h_pdoping_left,
                                                  h_pdoping_left,
                                                 )

    grid                    = simplexgrid(coord)

    ## set different regions in grid, doping profiles do not intersect
    ## n doped 
    cellmask!(grid, [0.0 * μm], [h_pdoping_left], regionAcceptor)          
    ## p doped                    
    cellmask!(grid, [h_pdoping_left], [h_pdoping_left + h_pdoping_left], regionAcceptor)    


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
    Nt                = 5e14                / (cm^3)   
    Nt_low            = Nt#/1e3                        
    mun_CIGS          = 100.0                * (cm^2) / (V * s)
    mup_CIGS          = 25                   * (cm^2) / (V * s)
    mun_ZnO           = 100                  * (cm^2) / (V * s)
    mup_ZnO           = 25                   * (cm^2) / (V * s)
    mut               = 0                    * (cm^2) / (V * s)  # no flux for traps
    εr_CIGS           = 13.6                 *  1.0              
    εr_ZnO            = 9                    *  1.0                
    T                 = 200.0                *  K

    An                = 4 * pi * q * mₑ * kB^2 / Planck_constant^3
    Ap                = 4 * pi * q * mₑ * kB^2 / Planck_constant^3
    vn                = An * T^2 / (q*Nc)
    vp                = Ap * T^2 / (q*Nv)
    barrier_right     = 0.7 * eV
    barrier_left      = 0.1 * eV

    ## recombination parameters
    #=
    ni_CIGS           = sqrt(Nc * Nv) * exp(-(Ec_CIGS - Ev_CIGS) / (2 * kB * T)) # intrinsic concentration
    n0_CIGS           = Nc * Boltzmann( (Et-Ec_CIGS) / (kB*T) )             # Boltzmann equilibrium concentration
    p0_CIGS           = ni_CIGS^2 / n0_CIGS                                      # Boltzmann equilibrium concentration
    ni_ZnO            = sqrt(Nc * Nv) * exp(-(Ec_ZnO - Ev_ZnO) / (2 * kB * T)) # intrinsic concentration
    n0_ZnO            = Nc * Boltzmann( (Et-Ec_ZnO) / (kB*T) )             # Boltzmann equilibrium concentration
    p0_ZnO            = ni_ZnO^2 / n0_ZnO                                      # Boltzmann equilibrium concentration
    =#
    Auger             = 1.0e-29  * cm^6 / s          # 1.0e-41 m^6 / s
    SRH_LifeTime      = 1.0e-3   * ns               
    Radiative         = 1.0e-10  * cm^3 / s          # 1.0e-16 m^3 / s
    G                 = 1.0e20   / (cm^3 * s)
    ## ???
    A_CIGS            = 1.0e5    / cm
    A_ZnO             = 0.0      / cm
    N0                = 1e17     / cm^2/s

    ## doping -- trap doping will not be set and thus automatically zero
    Na                = 1.0e15 / (cm^3)   

    ## we will impose this applied voltage on one boundary
    voltageMin   = voltageMin * V 
    voltageMax   = voltageMax * V

    println("*** done\n")
    ################################################################################
    println("Define System and fill in information about model")
    ################################################################################

    ## initialize Data instance and fill in data
    data                                = Data(grid, numberOfCarriers)

    ## possible choices: model_stationary, model_transient
    data.modelType                      = Stationary

    ## possible choices: Boltzmann, FermiDiracOneHalfBednarczyk, FermiDiracOneHalfTeSCA FermiDiracMinusOne, Blakemore
    data.F                             .= [FermiDiracOneHalfTeSCA, FermiDiracOneHalfTeSCA]

    data.bulkRecombination              = set_bulk_recombination(;iphin = iphin, iphip = iphip,
                                                                 bulk_recomb_Auger = false,
                                                                 bulk_recomb_radiative = false,
                                                                 bulk_recomb_SRH = false)

    #enable_traps!(data)
    
    ## Possible choices: GenerationNone, GenerationUniform, GenerationBeerLambert
    data.generationModel                = GenerationNone#GenerationBeerLambert

    ## Possible choices: OhmicContact, SchottkyContact (outer boundary) and InterfaceModelNone,
    ## InterfaceModelSurfaceReco (inner boundary).
    #data.boundary_type[bregionAcceptorLeft ]    = ohmic_contact#schottky_contact                       
    #data.boundary_type[bregionAcceptorRight]    = ohmic_contact#schottky_contact   
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
    params.chargeNumbers[iphin]                         =  -1
    params.chargeNumbers[iphip]                         =  1
    # params.chargeNumbers[iphit]                         =  1  # +1: hole trap is used

    for ibreg in 1:numberOfBoundaryRegions   # boundary region data
        params.bDensityOfStates[iphin, ibreg]           = Nc
        params.bDensityOfStates[iphip, ibreg]           = Nv
        # params.bDensityOfStates[iphit, ibreg]           = Nt_low
        # params.bBandEdgeEnergy[iphit, ibreg]            = Et
    end

    params.bBandEdgeEnergy[iphin, bregionAcceptorRight]         = Ec_CIGS
    params.bBandEdgeEnergy[iphip, bregionAcceptorRight]         = Ev_CIGS
    params.bBandEdgeEnergy[iphin, bregionAcceptorLeft]      = Ec_CIGS
    params.bBandEdgeEnergy[iphip, bregionAcceptorLeft]      = Ev_CIGS

    for ireg in 1:numberOfRegions           # interior region data

        params.dielectricConstant[ireg]                 = εr_CIGS*ε0     

        ## effective DOS, band-edge energy and mobilities
        params.densityOfStates[iphin, ireg]             = Nc
        params.densityOfStates[iphip, ireg]             = Nv
        params.bandEdgeEnergy[iphin, ireg]              = Ec_CIGS
        params.bandEdgeEnergy[iphip, ireg]              = Ev_CIGS
        params.mobility[iphin, ireg]                    = mup_CIGS
        params.mobility[iphip, ireg]                    = mup_CIGS

        ## recombination parameters
        params.recombinationRadiative[ireg]             = Radiative
        params.recombinationSRHLifetime[iphin, ireg]    = SRH_LifeTime
        params.recombinationSRHLifetime[iphip, ireg]    = SRH_LifeTime
        #params.recombinationSRHTrapDensity[iphin, ireg] = n0_CIGS
        #params.recombinationSRHTrapDensity[iphip, ireg] = p0_CIGS
        params.recombinationAuger[iphin, ireg]          = Auger
        params.recombinationAuger[iphip, ireg]          = Auger

        ## generation parameters
        params.generationAbsorption[ireg]               = A_CIGS
        params.generationIncidentPhotonFlux[ireg]       = N0
        params.generationUniform[ireg]                  = G
        
    end

    ## overwrite parameters in ZnO donor region
    #=
    params.generationUniform[regionAcceptorLeft]                  = 0.0      # only used if for "generation_uniform"
    params.generationAbsorption[regionAcceptorLeft]               = A_ZnO    # only used if for "generation_beer_lambert"
    params.generationIncidentPhotonFlux[regionAcceptorLeft]       = N0
    params.recombinationSRHTrapDensity[iphin, regionAcceptorLeft] = n0_ZnO
    params.recombinationSRHTrapDensity[iphip, regionAcceptorLeft] = p0_ZnO
    params.bandEdgeEnergy[iphin, regionAcceptorLeft]              = Ec_ZnO
    params.bandEdgeEnergy[iphip, regionAcceptorLeft]              = Ev_ZnO
    params.dielectricConstant[regionAcceptorLeft]                 = εr_ZnO 
    params.mobility[iphin, regionAcceptorLeft]                    = mun_ZnO
    params.mobility[iphip, regionAcceptorLeft]                    = mup_ZnO
    =#
    ## hole trap density only high in grain
    # params.densityOfStates[iphit, regionAcceptorLeft]          = Nt_low
    # params.densityOfStates[iphit, regionAcceptorRight]   = Nt_low
    # params.densityOfStates[iphit, regionAcceptorTrap]   = Nt
    # params.densityOfStates[iphit, regionAcceptorRight]  = Nt_low

    ## doping -- since we do not set any doping for the traps it is automatically zero
    params.doping[iphip, regionAcceptor]             = Na        
    # params.doping[iphip, regionAcceptor]            = Na        

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
    #set_schottky_contact!(ctsys, bregionAcceptorRight, appliedVoltage = 0.0)
    #set_schottky_contact!(ctsys, bregionAcceptorLeft , appliedVoltage = 0.0)
    #set_ohmic_contact!(ctsys, bregionAcceptorRight   , 0.0)
    #set_ohmic_contact!(ctsys, bregionAcceptorLeft, 0.0)
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

        ## for traps 
        # label_energy[1, iphit] = "\$E_{\\tau}-q\\psi\$"; label_energy[2, iphit] = "\$ - q \\varphi_{\\tau}\$"
        # label_density[iphit]   = "\$n_{\\tau}\$";        label_solution[iphit]  = "\$ \\varphi_{\\tau}\$"
        ## ##### set legend for plotting routines #####
        Plotter.figure()
        plot_energies(Plotter, ctsys, solution, "Equilibrium", label_energy)
        Plotter.figure()
        plot_densities(Plotter, ctsys, solution,"Equilibrium", label_density)

        # plot_solution(Plotter, ctsys, solution, "Equilibrium", label_solution)
        # Plotter.figure()
    end

    ################################################################################
    println("Stationary bias loop")
    ################################################################################
    
    ## set calculationType to OutOfEquilibrium for starting with respective simulation.
    data.calculationType = OutOfEquilibrium      # Rn = Rp = R, since the model type is stationary
    endVoltage           = voltageMin            # final bias value

    IV         = zeros(0)   
    #maxBias    = voltageMin  
    #maxBias2   = -1*voltageMin   
    biasSteps  = 101
    biasValues = collect(range(voltageMin, stop = voltageMax, length = biasSteps))
    if(!(0.0 in biasValues))
        append!(biasValues, 0.0)
        sort!(biasValues)
    end
    #print(biasValues)
    #biasValues2 = collect(range(0, stop= maxBias2, length = biasSteps))
    chargeDensities = zeros(0)
    #chargeDensities2 = zeros(0)
    

    w_device = 1.0    * cm  # width of device
    z_device = 1.0    * cm  # depth of device

    ## adjust Newton parameters
    control.tol_absolute      = 1.0e-10
    control.tol_relative      = 1.0e-10
    control.tol_round         = 1.0e-7
    control.damp_initial      = 0.5
    control.damp_growth       = 1.2
    control.max_iterations    = 30
    control.max_round         = 3

    AllSolution = Matrix{Float64}[]
    indexOfZero = indexin(0.0, biasValues)[1]
    i = indexOfZero

    #for i in eachindex(biasValues)
    while (i>0)
        Δu = biasValues[i] # bias

        ## set non equilibrium boundary condition
        #set_schottky_contact!(ctsys, bregionAcceptorLeft, appliedVoltage = Δu)
        set_contact!(ctsys, bregionAcceptorLeft, Δu = Δu)
     
        ## solve time step problems with timestep Δt
        solve!(solution, initialGuess, ctsys, control  = control, tstep = Inf)

        ## save IV data
        current = get_current_val(ctsys, solution)
        push!(IV, w_device * z_device * current)

        ## store CHARGE DENSITY in CIGS
        #push!(chargeDensities,chargeDensity(ctsys,solution)[regionAcceptorLeft])
        push!(chargeDensities,w_device * z_device *(charge_density(ctsys,solution)[regionAcceptor]))#+charge_density(ctsys,solution)[regionAcceptor]))

        push!(AllSolution, solution)

        initialGuess .= solution
        i=i-1
    end # bias loop 1

    # FileIO.save("steady_state_solution.jld2","steady_state_solution",solution)
    reverse!(IV)
    reverse!(chargeDensities)
    reverse!(AllSolution)

    ## plot energies and qFermi levels for voltageMin
    if plotting 
        plot_energies(Plotter, ctsys, solution, "bias \$\\Delta u\$ = $(endVoltage), \$ t=$(0)\$", label_energy)
        Plotter.figure()
                
    end
    
    i = indexOfZero+1
    initialGuess .= equilibriumSolution

    while (i<=length(biasValues))
        Δu = biasValues[i] # bias

        ## set non equilibrium boundary condition
        #set_schottky_contact!(ctsys, bregionAcceptorLeft, appliedVoltage = Δu)
        set_contact!(ctsys, bregionAcceptorLeft, Δu = Δu)
    
        ## solve time step problems with timestep Δt
        solve!(solution, initialGuess, ctsys, control  = control, tstep = Inf)

        ## save IV data
        current = get_current_val(ctsys, solution)
        push!(IV, w_device * z_device * current)

        ## store uncompensated CHARGE DENSITY in CIGS
        push!(chargeDensities,w_device * z_device *(charge_density(ctsys,solution)[regionAcceptor]))#+charge_density(ctsys,solution)[regionAcceptor]))
        push!(AllSolution, solution)

        initialGuess .= solution
        i=i+1
    end # bias loop 2

    FileIO.save("steady_state_solution.jld2","steady_state_solution",solution)

    println("*** done\n")

    ## compute static capacitance: check this is correctly computed

    staticCapacitance = diff(chargeDensities) ./ diff(biasValues)
    writedlm( "staticCapacitance.csv",  staticCapacitance, ',')
    writedlm( "chargeDensities.csv"  ,  chargeDensities  , ',')
    writedlm( "biasValues.csv"       ,  biasValues       , ',')

    ## plot solution and IV curve
    if plotting 
        plot_energies(Plotter, ctsys, solution, "bias \$\\Delta u\$ = $(voltageMax), \$ t=$(0)\$", label_energy)
        Plotter.figure()
        plot_densities(Plotter, ctsys, solution,"bias \$\\Delta u\$ = $(voltageMax), \$ t=$(0)\$", label_density)
        Plotter.figure()
        # plot_solution(Plotter, ctsys, solution, "bias \$\\Delta u\$ = $(endVoltage), \$ t=$(0)\$", label_solution)
        # Plotter.figure()
        plot_IV(Plotter, biasValues,IV, biasValues[end], plotGridpoints = true)
        Plotter.figure()
        plot_IV(Plotter, biasValues,chargeDensities, biasValues[end], plotGridpoints = true)
        Plotter.title("Charge density in donor region")
        Plotter.ylabel("Charge density [C]")
        Plotter.figure()
        plot_IV(Plotter, biasValues,abs.(staticCapacitance)*1e9, biasValues[end-1], plotGridpoints = true)
        Plotter.title("Static capacitance in donor region")
        Plotter.ylabel("Static capacitance [nF]")
               
    end
    Plotter.figure()
    plt.yscale("log")
    plot_IV(Plotter, biasValues,abs.(staticCapacitance)*1e9, biasValues[end-1], plotGridpoints = true)
    Plotter.title("Static capacitance in donor region")
    Plotter.ylabel("Static capacitance [nF]")

    println("*** done\n")

    ################################################################################
    println("Dynamic Capacitance")
    ################################################################################
    data.modelType = Transient
    ctsys = System(grid, data, unknown_storage=unknown_storage)
    NumCapPoints = length(biasValues)
    #1/τC serves as "frequency" in e.g. Impedance calculations. In more complicated structures results can strongly depend on τC! 
    τC = 1e-6 * s
    ΔV = 1e-3 * V
    Capacitance = zeros(NumCapPoints)

    for iCapMeas in eachindex(Capacitance)
            initialGuess .= AllSolution[iCapMeas]
            # charge_den_before = w_device * z_device * (charge_density(ctsys, initialGuess)[regionAcceptor])
            set_contact!(ctsys, bregionAcceptorLeft, Δu=biasValues[iCapMeas] + ΔV)
            solve!(solution, initialGuess, ctsys, control=control, tstep=τC)
            charge_den_after_plus   = w_device * z_device * (charge_density(ctsys, solution)[regionAcceptor])
            set_contact!(ctsys, bregionAcceptorLeft, Δu=biasValues[iCapMeas] - ΔV)
            solve!(solution, initialGuess, ctsys, control=control, tstep=τC)
            charge_den_after_minus  = w_device * z_device * (charge_density(ctsys, solution)[regionAcceptor])
            Capacitance[iCapMeas] = abs(   charge_den_after_plus - charge_den_after_minus   ) / (2*ΔV)

    end
    Plotter.figure()
    Plotter.plot(biasValues, Capacitance*1e9)
    plt.title("Dynamic Capacitance from IV curve")
    plt.xlabel("u[V]")
    plt.ylabel("C[nF]")


    testval = solution[data.index_psi, 10]
    return testval
end #  main

function test()
    testval = 1.3214196490674017
    main(test = true, unknown_storage=:dense) ≈ testval && main(test = true, unknown_storage=:sparse) ≈ testval
end

println("This message should show when this module has successfully recompiled.")


end # module



