#=
# CIGS pn junction: stationary with trap and Schottky contacts.
([source code](SOURCE_URL))

Simulating stationary charge transport in a pn junction with hole traps and a Schottky boundary condition.
=#

module CIGS_stationary_with_traps_2D_simple_grid

using ChargeTransport
using ExtendableGrids
using PyPlot
using SimplexGridFactory
using Triangulate

function main(;n = 3, Plotter = PyPlot, plotting = false, verbose = false, test = false, unknown_storage=:sparse)

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    ## region numbers
    regionZnO            = 1                           # n doped region
    regionCIGS           = 2                           # p doped region
    regionGrainBoundary  = 3                           # p doped region with trap
    regions             = [regionZnO, regionCIGS, regionGrainBoundary]
    numberOfRegions     = length(regions)

    ## boundary region numbers
    bregionZnO              = 1
    bregionCIGSRight        = 2
    bregionNoFlux           = 3  
    bregions                = [bregionZnO, bregionCIGSRight, bregionNoFlux]
    numberOfBoundaryRegions = length(bregions)

    ## grid
    width_ZnO             = 0.5 * μm 
    width_CIGS_Left       = 1.0 * μm 
    width_grain_boundary  = 0.01 * μm
    width_CIGS_Right      = 1.0 * μm
    height                = 1.0 * μm

    x0              = 0.0 * μm
    δ               = 3*n        # the larger, the finer the mesh
    t               = 0.5*(cm)/δ # tolerance for geomspace and glue (with factor 10)

    coord_ZnO       = collect(range(x0, width_ZnO, step=width_ZnO/(δ)))
    coord_CIGS_Left = collect(range(width_ZnO, 
                                    width_ZnO + width_CIGS_Left,
                                    step=width_CIGS_Left/(δ)))
    coord_grain_bou = collect(range(width_ZnO + width_CIGS_Left, 
                                    width_ZnO + width_CIGS_Left+width_grain_boundary,
                                    step=width_grain_boundary/(δ)))
    coord_CIGS_Right= collect(range(width_ZnO + width_CIGS_Left+width_grain_boundary, 
                                    width_ZnO + width_CIGS_Left+width_grain_boundary+width_CIGS_Right, 
                                    step=width_CIGS_Right/(δ)))

    coord_length    = glue(coord_ZnO   ,    coord_CIGS_Left ,  tol=10*t)
    coord_length    = glue(coord_length,    coord_grain_bou ,  tol=10*t)
    coord_length    = glue(coord_length,    coord_CIGS_Right,  tol=10*t)

    height_L        = geomspace(0.0, height/2, height/(0.5*δ), height/(0.5*δ))
    height_R        = geomspace(height/2, height, height/(0.5*δ), height/(0.5*δ))
    coord_height    = glue(height_L, height_R, tol = 10*t)

    grid            = simplexgrid(coord_length, coord_height)

    ## specify inner regions
    cellmask!(grid, [0.0                                           , 0.0],[width_ZnO                                                                       , height],regionZnO          , tol = 1.0e-18)  # p-doped region   = 1
    cellmask!(grid, [width_ZnO                                     , 0.0],[width_ZnO+width_CIGS_Left                                                       , height],regionCIGS         , tol = 1.0e-18)  # p-doped region   = 1
    cellmask!(grid, [width_ZnO+width_CIGS_Left                     , 0.0],[width_ZnO+width_CIGS_Left+width_CIGS_Right+width_grain_boundary                 , height],regionGrainBoundary, tol = 1.0e-18)     # n-doped region   = 3
    cellmask!(grid, [width_ZnO+width_CIGS_Left+width_grain_boundary, 0.0],[width_ZnO+width_CIGS_Left+width_CIGS_Right+width_grain_boundary+width_CIGS_Right, height],regionCIGS         , tol = 1.0e-18) 

    ## specifiy outer regions
    ## metal interfaces
    bfacemask!(grid, [0.0, 0.0], [0.0, height], bregionZnO) # BregionNumber = 1
    bfacemask!(grid, [width_CIGS_Left+width_CIGS_Right+width_grain_boundary+width_CIGS_Right, 0.0],
                     [width_CIGS_Left+width_CIGS_Right+width_grain_boundary+width_CIGS_Right, height],
                     bregionCIGSRight) # BregionNumber = 2

    ## no flux interfaces [xmin, ymin], [xmax, ymax]
    bfacemask!(grid, [0.0, 0.0], 
                     [width_CIGS_Left+width_CIGS_Right+width_grain_boundary+width_CIGS_Right, 0.0], 
                     bregionNoFlux) # BregionNumber = 3
    bfacemask!(grid, [0.0, height], 
                     [width_CIGS_Left+width_CIGS_Right+width_grain_boundary+width_CIGS_Right, height], 
                     bregionNoFlux) # # BregionNumber = 3

    if plotting
        Plotter.figure()
        gridplot(grid, Plotter = Plotter, legend=:lt)
        Plotter.title("Grid")

    end
    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define physical parameters and model")
    end
    ################################################################################

    iphin             = 1 # index electron quasi Fermi potential
    iphip             = 2 # index hole quasi Fermi potential
    # iphit             = 3 # index trap quasi Fermi potential
    numberOfCarriers  = 2 # electrons, holes and this time: no traps

    ## physical data
    # Ec                = 1.424                *  eV
    # Ev                = 0.0                  *  eV
    # Et                = 0.6                  *  eV

    Ec_CIGS           = 3.4                  *  eV
    Ev_CIGS           = 2.3                  *  eV
    Ec_ZnO            = 3.4                  *  eV
    Ev_ZnO            = 0.0                  *  eV
    Et                = 2.8                  *  eV

    Nc                = 4.351959895879690e17 / (cm^3)
    Nv                = 9.139615903601645e18 / (cm^3)
    Nt                = 1e18                   / (cm^3)
    #Nt_low            = Nt#/1e3
    mun_CIGS          = 100.0                * (cm^2) / (V * s)
    mup_CIGS          = 25                   * (cm^2) / (V * s)
    mun_ZnO           = 100                  * (cm^2) / (V * s)
    mup_ZnO           = 25                   * (cm^2) / (V * s)
    mut               = 0                    * (cm^2) / (V * s)  # no flux for traps
    εr_CIGS           = 13.6                 *  1.0
    εr_ZnO            = 9                    *  1.0
    T                 = 300.0                *  K

    An                = 4 * pi * q * mₑ * kB^2 / Planck_constant^3
    Ap                = 4 * pi * q * mₑ * kB^2 / Planck_constant^3
    vn                = An * T^2 / (q*Nc)
    vp                = Ap * T^2 / (q*Nv)
    barrier           = 0.7* eV

    ## recombination parameters
    ni_CIGS           = sqrt(Nc * Nv) * exp(-(Ec_CIGS - Ev_CIGS) / (2 * kB * T)) # intrinsic concentration
    n0_CIGS           = Nc * Boltzmann( (Et-Ec_CIGS) / (kB*T) )             # Boltzmann equilibrium concentration
    p0_CIGS           = ni_CIGS^2 / n0_CIGS                                      # Boltzmann equilibrium concentration
    ni_ZnO            = sqrt(Nc * Nv) * exp(-(Ec_ZnO - Ev_ZnO) / (2 * kB * T)) # intrinsic concentration
    n0_ZnO            = Nc * Boltzmann( (Et-Ec_ZnO) / (kB*T) )             # Boltzmann equilibrium concentration
    p0_ZnO            = ni_ZnO^2 / n0_ZnO                                      # Boltzmann equilibrium concentration
    Auger             = 1.0e-29  * cm^6 / s          # 1.0e-41 m^6 / s
    SRH_LifeTime      = 1.0e-3   * ns
    Radiative         = 1.0e-10  * cm^3 / s          # 1.0e-16 m^3 / s
    G                 = 1.0e20   / (cm^3 * s)
    A_CIGS            = 1.0e5    / cm
    A_ZnO             = 0.0      / cm
    N0                = 1e17     / cm^2/s

    ## doping -- trap doping will not be set and thus automatically zero
    Nd                = 1.0e18 / (cm^3)
    Na                = 5.0e15 / (cm^3)

    ## we will impose this applied voltage on one boundary
    voltageAcceptor   = 1 * V

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Define System and fill in information about model")
    end
    ################################################################################

    ## initialize Data instance and fill in data
    data                                = Data(grid, numberOfCarriers)

    ## Possible choices: Stationary, Transient
    data.modelType                      = Stationary

    ## Possible choices: Boltzmann, FermiDiracOneHalfBednarczyk, FermiDiracOneHalfTeSCA FermiDiracMinusOne, Blakemore
    data.F                             .= [FermiDiracOneHalfTeSCA, FermiDiracOneHalfTeSCA]

    data.bulkRecombination              = set_bulk_recombination(;iphin = iphin, iphip = iphip,
                                                                 bulk_recomb_Auger = true,
                                                                 bulk_recomb_radiative = true,
                                                                 bulk_recomb_SRH = true)


    # pass trap data in stationary setting since there is no separate trap species
    zt = 1
    Nt_vector = similar(data.params.densityOfStates[zt,:])
    Nt_vector[regionGrainBoundary] = Nt
    NT = [0, 0, Nt]
    add_trap_density!(data=data, zt = zt, Nt = NT)

    ## Possible choices: GenerationNone, GenerationUniform, GenerationBeerLambert
    data.generationModel                = GenerationNone#GenerationBeerLambert

    ## Possible choices: OhmicContact, SchottkyContact (outer boundary) and InterfaceModelNone,
    ## InterfaceModelSurfaceReco (inner boundary).
    data.boundaryType[bregionCIGSRight]  = OhmicContact #SchottkyContact
    data.boundaryType[bregionZnO]        = OhmicContact
    #data.boundaryType[bregionInner]      = InterfaceModelNone

    ## Choose flux discretization scheme: ScharfetterGummel, ScharfetterGummelGraded,
    ## ExcessChemicalPotential, ExcessChemicalPotentialGraded, DiffusionEnhanced, GeneralizedSG
    data.fluxApproximation              .= ExcessChemicalPotential

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define Params and fill in physical parameters")
    end
    ################################################################################

    ## physical parameters
    params                                              = Params(grid, numberOfCarriers)
    params.temperature                                  = T
    params.UT                                           = (kB * params.temperature) / q
    params.chargeNumbers[iphin]                         = -1
    params.chargeNumbers[iphip]                         =  1
    #params.chargeNumbers[iphit]                         =  1  # +1: hole trap is used

    for ibreg in 1:numberOfBoundaryRegions   # boundary region data
        params.bDensityOfStates[iphin, ibreg]           = Nc
        params.bDensityOfStates[iphip, ibreg]           = Nv
        # params.bDensityOfStates[iphit, ibreg]           = Nt_low
        # params.bBandEdgeEnergy[iphit, ibreg]            = Et
    end

    params.bBandEdgeEnergy[iphin, bregionZnO]         = Ec_ZnO
    params.bBandEdgeEnergy[iphip, bregionZnO]         = Ev_ZnO
    params.bBandEdgeEnergy[iphin, bregionCIGSRight]      = Ec_CIGS
    params.bBandEdgeEnergy[iphip, bregionCIGSRight]      = Ev_CIGS

    for ireg in [regionZnO, regionCIGS, regionGrainBoundary]          # interior region data

        params.dielectricConstant[ireg]                 = εr_CIGS*ε0

        ## effective DOS, band-edge energy and mobilities
        params.densityOfStates[iphin, ireg]             = Nc
        params.densityOfStates[iphip, ireg]             = Nv
        params.bandEdgeEnergy[iphin, ireg]              = Ec_CIGS
        params.bandEdgeEnergy[iphip, ireg]              = Ev_CIGS
        # params.bandEdgeEnergy[iphit, ireg]              = Et
        params.mobility[iphin, ireg]                    = mun_CIGS
        params.mobility[iphip, ireg]                    = mup_CIGS
        # params.mobility[iphit, ireg]                    = mut

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

    ## overwrite parameters in ZnO donor region
    params.generationUniform[regionZnO]                  = 0.0      # only used if for "generation_uniform"
    params.generationAbsorption[regionZnO]               = A_ZnO    # only used if for "generation_beer_lambert"
    params.generationIncidentPhotonFlux[regionZnO]       = N0
    params.recombinationSRHTrapDensity[iphin, regionZnO] = n0_ZnO
    params.recombinationSRHTrapDensity[iphip, regionZnO] = p0_ZnO
    params.bandEdgeEnergy[iphin, regionZnO]              = Ec_ZnO
    params.bandEdgeEnergy[iphip, regionZnO]              = Ev_ZnO
    params.dielectricConstant[regionZnO]                 = εr_ZnO*ε0
    params.mobility[iphin, regionZnO]                    = mun_ZnO
    params.mobility[iphip, regionZnO]                    = mup_ZnO

    ## hole trap density only high in grain
    # params.densityOfStates[iphit, regionZnO]          = Nt_low
    # params.densityOfStates[iphit, regionCIGS]   = Nt_low
    # params.densityOfStates[iphit, regionGrainBoundary]   = Nt
    # params.densityOfStates[iphit, regionCIGS]  = Nt_low

    ## doping -- since we do not set any doping for the traps it is automatically zero
    params.doping[iphin, regionZnO]             = Nd
    params.doping[iphip, regionCIGS]            = Na
    params.doping[iphip, regionGrainBoundary]   = Na
    params.doping[iphip, regionCIGS]            = Na

    ## boundary doping
    params.bDoping[iphin, bregionZnO]                 = Nd
    params.bDoping[iphip, bregionCIGSRight]              = Na

    ## values for the schottky contacts
    params.SchottkyBarrier[bregionCIGSRight]             = barrier
    params.bVelocity[iphin,bregionCIGSRight]             = vn
    params.bVelocity[iphip,bregionCIGSRight]             = vp

    data.params                                         = params
    ctsys                                               = System(grid, data, unknown_storage=unknown_storage)

    if test == false
        show_params(ctsys)
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define outerior boundary conditions")
    end
    ################################################################################

    ## set zero voltage ohmic contacts for all charge carriers at all outerior boundaries.
    set_contact!(ctsys, bregionCIGSRight, Δu = 0.0)
    set_contact!(ctsys, bregionZnO,    Δu = 0.0)

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define control parameters for Newton solver")
    end
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

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Compute solution in thermodynamic equilibrium")
    end
    ################################################################################

    ## initialize solution and starting vectors
    initialGuess          = unknowns(ctsys)
    solution              = unknowns(ctsys)

    ## solve thermodynamic equilibrium and update initial guess
    solution              = equilibrium_solve!(ctsys, control = control, nonlinear_steps = 20)
    initialGuess         .= solution

    if test == false
        println("*** done\n")
    end

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
    if test == false
        println("Stationary bias loop")
    end
    ################################################################################

    ## set calculationType to OutOfEquilibrium for starting with respective simulation.
    data.calculationType = OutOfEquilibrium      # Rn = Rp = R, since the model type is stationary
    endVoltage           = voltageAcceptor       # final bias value

    IV         = zeros(0)
    maxBias    = voltageAcceptor
    biasSteps  = 52
    biasValues = collect(range(0, stop = maxBias, length = biasSteps))
    chargeDensities = zeros(0)

    w_device = 0.5    * μm  # width of device
    z_device = 1.0e-4 * cm  # depth of device

    ## adjust Newton parameters
    control.tol_absolute      = 1.0e-10
    control.tol_relative      = 1.0e-10
    control.tol_round         = 1.0e-7
    control.damp_initial      = 0.5
    control.damp_growth       = 1.2
    control.max_iterations    = 30
    control.max_round         = 3

    for i in eachindex(biasValues)

        Δu = biasValues[i] # bias

        ## Apply new voltage: set non equilibrium boundary conditions
        set_contact!(ctsys, bregionCIGSRight, Δu = Δu)

        ## increase generation rate with bias
        ctsys.data.λ2 = 10.0^(-biasSteps + i)

        if test == false
            println("bias: Δu = $(Δu)")
        end

        ## solve time step problems with timestep Δt
        solve!(solution, initialGuess, ctsys, control  = control, tstep = Inf)

        ## save IV data
        current = get_current_val(ctsys, solution)
        push!(IV, w_device * z_device * current)

        ## store charge density in donor region (ZnO)
        push!(chargeDensities,charge_density(ctsys,solution)[regionZnO])

        initialGuess .= solution

    end # bias loop

    if test == false
        println("*** done\n")
    end

    ## compute static capacitance: check this is correctly computed
    staticCapacitance = diff(chargeDensities) ./ diff(biasValues)

    ## plot solution and IV curve
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
        
        Plotter.title("Band Edge Energies and qFermi levels at \$ $(endVoltage) \$V")
        Plotter.xlabel("length [m]")
        Plotter.ylabel("height [m]")
        Plotter.zlabel("Energy [eV]")
        Plotter.tight_layout()

        #Plot densities
        Plotter.figure()
        #TO DO: make the labels visible 
        Plotter.surf(X[:], Y[:], log10.(1.0e-6 .*Nc.*exp.((-(Ec_CIGS/q .- solution[ipsi, :])-solution[iphin, :])/(params.UT)))) # electron density
        Plotter.surf(X[:], Y[:], log10.(1.0e-6 .*Nv.*exp.(( (Ev_CIGS/q .- solution[ipsi, :])+solution[iphip, :])/(params.UT)))) #hole density
        
        Plotter.title("electrons and holes densities at \$ $(endVoltage) \$V")
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
    testval = 23.96747612965515
    main(test = true, unknown_storage=:dense) ≈ testval && main(test = true, unknown_storage=:sparse) ≈ testval
end

if test == false
    println("This message should show when this module has successfully recompiled.")
end


end # module

