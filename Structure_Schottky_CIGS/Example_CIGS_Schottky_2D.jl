#=
# CIGS pn junction: stationary with trap and Schottky contacts.
Simulating stationary charge transport in a pn junction with hole traps and a Schottky boundary condition.
=#

module Example_CIGS_Schottky_2D

using ChargeTransport
using ExtendableGrids
using PyPlot
using DelimitedFiles

function main(;n = 3, voltageMin=-0.5, voltageMax=0.1, Plotter = PyPlot, plotting = false, verbose = false, test = false, unknown_storage=:sparse)

    ################################################################################
    println("Set up grid and regions")
    ################################################################################

    ## region numbers
    regionAcceptorLeft  = 1                           # n doped region
    regionAcceptorRight = 2                           # p doped region
    regions             = [regionAcceptorLeft, regionAcceptorRight]
    numberOfRegions     = length(regions)

    ## boundary region numbers
    bregionAcceptorLeft     = 1
    bregionAcceptorRight    = 2
    bregionNoFlux           = 3
    bregions                = [bregionAcceptorRight, bregionAcceptorLeft, bregionNoFlux]
    numberOfBoundaryRegions = length(bregions)

    ## grid
    h_pdoping       = 1 * μm
    h_intrinsic     = 0.00e-5 * cm
    h_ndoping       = 1 * μm
    height          = 1 * μm

    x0              = 0.0 * μm
    δ               = 3*n        # the larger, the finer the mesh
    t               = 0.5*(cm)/δ # tolerance for geomspace and glue (with factor 10)

    # coord_p_u       = collect(range(x0, h_pdoping, step=h_pdoping/(0.3*δ)))
    coord_p_u       = collect(range(x0, h_pdoping, step=h_pdoping/(δ)))
    coord_n_u       = collect(range(h_pdoping, h_pdoping+h_ndoping, step=h_pdoping/(δ)))

    coord_length    = glue(coord_p_u,    coord_n_u,  tol=10*t)

    height_L        = geomspace(0.0, height/2, height/(0.5*δ), height/(0.5*δ))
    height_R        = geomspace(height/2, height, height/(0.5*δ), height/(0.5*δ))
    coord_height    = glue(height_L, height_R, tol = 10*t)

    grid            = simplexgrid(coord_length, coord_height)

    ## specify inner regions
    cellmask!(grid, [0.0, 0.0],                     [h_pdoping, height],                           regionAcceptorLeft, tol = 1.0e-18)  # p-doped region   = 1
    #cellmask!(grid, [h_pdoping, 0.0],               [h_pdoping + h_intrinsic, height],             regionIntrinsic, tol = 1.0e-18) # intrinsic region = 2
    cellmask!(grid, [h_pdoping , 0.0], [h_pdoping  + h_ndoping, height], regionAcceptorRight, tol = 1.0e-18)     # n-doped region   = 3

    ## specifiy outer regions
    ## metal interfaces
    bfacemask!(grid, [0.0, 0.0], [0.0, height], bregionAcceptorLeft) # BregionNumber = 1
    bfacemask!(grid, [h_pdoping  + h_ndoping, 0.0], [h_pdoping  + h_ndoping, height], bregionAcceptorRight) # BregionNumber = 2

    ## no flux interfaces [xmin, ymin], [xmax, ymax]
    bfacemask!(grid, [0.0, 0.0], [h_pdoping  + h_ndoping, 0.0], bregionNoFlux) # BregionNumber = 3
    bfacemask!(grid, [0.0, height], [h_pdoping  + h_ndoping, height], bregionNoFlux) # # BregionNumber = 3


    if plotting
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
    mun_CIGS          = 100.0                * (cm^2) / (V * s)
    mup_CIGS          = 25                   * (cm^2) / (V * s)
    εr_CIGS           = 13.6                 *  1.0                       
    T                 = 300.0                *  K

    An                = 4 * pi * q * mₑ * kB^2 / Planck_constant^3
    Ap                = 4 * pi * q * mₑ * kB^2 / Planck_constant^3
    vn                = An * T^2 / (q*Nc)
    vp                = Ap * T^2 / (q*Nv)
    barrier_right     = 0.7 * eV
    barrier_left      = 0.5 * eV

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
    data.boundaryType[bregionAcceptorRight]    = OhmicContact   
    
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
    params.doping[iphip, regionAcceptorLeft]             = Na        
    params.doping[iphip, regionAcceptorRight]            = Na        

    ## boundary doping
    params.bDoping[iphip, bregionAcceptorRight]          = Na        
    params.bDoping[iphip, bregionAcceptorLeft]           = Na   
    
    ## values for the schottky contacts
    params.SchottkyBarrier[bregionAcceptorLeft]             = barrier_left
    params.bVelocity[iphin,bregionAcceptorLeft]             = vn 
    params.bVelocity[iphip,bregionAcceptorLeft]             = vp 

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

        set_contact!(ctsys, bregionAcceptorRight, Δu = Δu)

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
        push!(chargeDensities,w_device * z_device *(charge_density(ctsys,solution)[regionAcceptorLeft]+charge_density(ctsys,solution)[regionAcceptorRight]))

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
        set_contact!(ctsys, bregionAcceptorRight, Δu = Δu)
    
        ## solve time step problems with timestep Δt
        solve!(solution, initialGuess, ctsys, control  = control, tstep = Inf)

        ## save IV data
        current = get_current_val(ctsys, solution)
        push!(IV, w_device * z_device * current)

        ## store uncompensated CHARGE DENSITY in CIGS
        push!(chargeDensities,w_device * z_device *(charge_density(ctsys,solution)[regionAcceptorLeft]+charge_density(ctsys,solution)[regionAcceptorRight]))

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

end #  main

function test()
    testval = 1.3214196490674017
    main(test = true, unknown_storage=:dense) ≈ testval && main(test = true, unknown_storage=:sparse) ≈ testval
end

println("This message should show when this module has successfully recompiled.")


end # module



