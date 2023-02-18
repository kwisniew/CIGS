#=
# CIGS pn junction: stationary with trap and Schottky contacts.
Simulating stationary charge transport in a pn junction with hole traps and a Schottky boundary condition.
=#

module DLTS_measutements_simulations

using ChargeTransport
using ExtendableGrids
using PyPlot
using DelimitedFiles
using FileIO, JLD2

## function to initialize the grid for a possible extension to other p-i-n devices.
function initialize_pin_grid(refinementfactor, h_pdoping, h_pdoping_right)
    coord_pdoping_left = collect(range(0.0, stop=h_pdoping, length=3 * refinementfactor))
    coord_pdoping_right = collect(range(h_pdoping, stop=(h_pdoping + h_pdoping_right), length=3 * refinementfactor))

    coord = glue(coord_pdoping_left, coord_pdoping_right)

    return coord
end

#In this script voltage is applied to the left side of the domain, what means that "+" means reverse voltage, and "-" means forward voltage
function main(; n=3, voltageStep=0.5, Plotter=PyPlot, plotting=false, verbose=false, test=false, unknown_storage=:sparse, makegif=false)

    ################################################################################
    println("Set up grid and regions")
    ################################################################################


    # Plotter.plot(1:10,1:10)

    # return
    ## region numbers
    regionAcceptor = 1  # p doped region
    regions = [regionAcceptor]
    numberOfRegions = length(regions)

    ## boundary region numbers
    bregionAcceptorLeft = 1
    bregionAcceptorRight = 2
    bregions = [bregionAcceptorRight, bregionAcceptorLeft]
    numberOfBoundaryRegions = length(bregions)

    ## grid
    refinementfactor = 2^(n - 1)
    h_pdoping = 2.5 * μm

    coord = initialize_pin_grid(refinementfactor,
        h_pdoping / 2,
        h_pdoping / 2,
    )

    grid = simplexgrid(coord)

    ## set different regions in grid, doping profiles do not intersect
    ## p doped                    
    cellmask!(grid, [0.0 * μm], [h_pdoping], regionAcceptor)


    if plotting
        gridplot(grid, Plotter=Plotter, legend=:lt)
        Plotter.title("Grid")
        Plotter.figure()

    end

    println("*** done\n")

    ################################################################################
    println("Define physical parameters and model")
    ################################################################################
    ## We will impose this applied voltage on one boundary. The positive sign corresponds to "reverse" direction
    voltageStep = voltageStep * V
    #Inverse of measurements frequency (typical 1MHz)
    τC = 1e-6*s
    #Time when the capacitance is measured 1st time
    t1 = 1e-5 * s
    #Time when the capacitance is measured 2nd time
    t2 = 1e-3 * s
    #Time for voltage to achieved final value (since in DLTS measurement setup it is almost step function: the faster, the better) 
    timeStep_for_bias = 1e-11 * s
    #Temperatures at which the measurements where performed
    Temperatures = range(50,stop=250,length=51) * K

    #Initialization:
    #Initiate vactor of DLTS measurements points
    DLTS_values  = zeros(length(Temperatures))
    #How many capacitance point we want at each temperature (typically 2: C(t1) and C(t2))
    NumCapPoints = 2
    #Capacitance at t1:
    C1=0
    #Capacitance at t2:
    C2=0
    

    println("*** done\n")

    ################################################################################
    println("Define physical parameters and model")
    ################################################################################

    iphin = 1 # index electron quasi Fermi potential
    iphip = 2 # index hole quasi Fermi potential
    numberOfCarriers = 2 # electrons, holes and no traps

    Ec_CIGS = 1.1 * eV
    Ev_CIGS = 0.0 * eV

    Nc = 4.351959895879690e17 / (cm^3)
    Nv = 9.139615903601645e18 / (cm^3)
    Nt = 5e14 / (cm^3)
    mun_CIGS = 100.0 * (cm^2) / (V * s)
    mup_CIGS = 25 * (cm^2) / (V * s)
    εr_CIGS = 13.6 * 1.0
    T = Temperatures[1]

    An = 4 * pi * q * mₑ * kB^2 / Planck_constant^3
    Ap = 4 * pi * q * mₑ * kB^2 / Planck_constant^3
    vn = An * T^2 / (q * Nc)
    vp = Ap * T^2 / (q * Nv)
    barrier_right = 0.7 * eV
    barrier_left = 0.1 * eV

    ## recombination parameters
    Auger = 1.0e-29 * cm^6 / s          # 1.0e-41 m^6 / s
    SRH_LifeTime = 1.0e-3 * ns
    Radiative = 1.0e-10 * cm^3 / s          # 1.0e-16 m^3 / s
    G = 1.0e20 / (cm^3 * s)
    ## ???
    A_CIGS = 1.0e5 / cm
    N0 = 1e17 / cm^2 / s

    ## doping -- trap doping will not be set and thus automatically zero
    Na = 1.0e15 / (cm^3)

    # # biasVal=0.5
    biasValues = range(0.0, stop=voltageStep, length=11)
    tend_bias = (length(biasValues) - 1) * timeStep_for_bias
    tend = t2
    # tend = 1e-3 * s #1e-7
    ## Define scan protocol function
    function DLTSscanProtocol(t)

        if 0.0 <= t && t <= tend_bias
            biasVal = voltageStep * (t / tend_bias)
        elseif t > tend_bias && t <= tend
            biasVal = voltageStep#+(voltageStep*t/tend)
        else
            biasVal = 0.0
        end

        return biasVal

    end

    function RDLTSscanProtocol(t)

        if 0.0 <= t && t <= tend_bias
            biasVal = voltageStep * (1 - (t / tend_bias) )
        elseif t > tend_bias && t <= tend
            biasVal = 0.0
        else
            biasVal = voltageStep
        end

        return biasVal

    end

    # contactVoltageFunction = [zeroVoltage, scanProtocol]
    DLTScontactVoltageFunction  = [ DLTSscanProtocol, zeroVoltage]
    RDLTScontactVoltageFunction = [RDLTSscanProtocol, zeroVoltage]

    #Other physical dimension of a sample
    w_device = 1.0 * cm  # width of device
    z_device = 1.0 * cm  # depth of device
    # test = 0.0:1e-13:tend
    # PyPlot.plot(test, scanProtocol.(test))
    # return

    println("*** done\n")
    ################################################################################
    println("Define System and fill in information about model")
    ################################################################################

    ## initialize Data instance and fill in data
    data = Data(grid, numberOfCarriers, contactVoltageFunction=DLTScontactVoltageFunction)
    ## Possible choices: Stationary, Transient
    data.modelType = Transient

    ## possible choices: Boltzmann, FermiDiracOneHalfBednarczyk, FermiDiracOneHalfTeSCA FermiDiracMinusOne, Blakemore
    data.F .= [FermiDiracOneHalfTeSCA, FermiDiracOneHalfTeSCA]

    data.bulkRecombination = set_bulk_recombination(; iphin=iphin, iphip=iphip,
        bulk_recomb_Auger=false,
        bulk_recomb_radiative=false,
        bulk_recomb_SRH=false)

    #enable_traps!(data)

    ## Possible choices: GenerationNone, GenerationUniform, GenerationBeerLambert
    # NOTE: GenerationBeerLambert doesn't work in 2D case!
    data.generationModel = GenerationNone#GenerationBeerLambert

    ## Possible choices: OhmicContact, SchottkyContact (outer boundary) and InterfaceModelNone,
    ## InterfaceModelSurfaceReco (inner boundary).
    data.boundaryType[bregionAcceptorLeft ] = SchottkyContact
    data.boundaryType[bregionAcceptorRight] = SchottkyContact#OhmicContact   

    ## Choose flux discretization scheme: ScharfetterGummel, ScharfetterGummelGraded,
    ## ExcessChemicalPotential, ExcessChemicalPotentialGraded, DiffusionEnhanced, GeneralizedSG
    data.fluxApproximation .= ExcessChemicalPotential

    println("*** done\n")

    ################################################################################
    println("Define Params and fill in physical parameters")
    ################################################################################

    ## physical parameters
    params = Params(grid, numberOfCarriers)
    params.temperature = T
    params.UT = (kB * params.temperature) / q
    params.chargeNumbers[iphin] = -1
    params.chargeNumbers[iphip] = 1

    for ibreg in 1:numberOfBoundaryRegions   # boundary region data
        params.bDensityOfStates[iphin, ibreg] = Nc
        params.bDensityOfStates[iphip, ibreg] = Nv
        # params.bDensityOfStates[iphit, ibreg]           = Nt_low
        # params.bBandEdgeEnergy[iphit, ibreg]            = Et
    end

    params.bBandEdgeEnergy[iphin, bregionAcceptorRight] = Ec_CIGS
    params.bBandEdgeEnergy[iphip, bregionAcceptorRight] = Ev_CIGS
    params.bBandEdgeEnergy[iphin, bregionAcceptorLeft] = Ec_CIGS
    params.bBandEdgeEnergy[iphip, bregionAcceptorLeft] = Ev_CIGS

    for ireg in 1:numberOfRegions           # interior region data

        params.dielectricConstant[ireg] = εr_CIGS * ε0

        ## effective DOS, band-edge energy and mobilities
        params.densityOfStates[iphin, ireg] = Nc
        params.densityOfStates[iphip, ireg] = Nv
        params.bandEdgeEnergy[iphin, ireg] = Ec_CIGS
        params.bandEdgeEnergy[iphip, ireg] = Ev_CIGS
        params.mobility[iphin, ireg] = mun_CIGS
        params.mobility[iphip, ireg] = mup_CIGS

        ## recombination parameters
        params.recombinationRadiative[ireg] = Radiative
        params.recombinationSRHLifetime[iphin, ireg] = SRH_LifeTime
        params.recombinationSRHLifetime[iphip, ireg] = SRH_LifeTime
        params.recombinationAuger[iphin, ireg] = Auger
        params.recombinationAuger[iphip, ireg] = Auger

        ## generation parameters
        params.generationAbsorption[ireg] = A_CIGS
        params.generationIncidentPhotonFlux[ireg] = N0
        params.generationUniform[ireg] = G

    end

    ## doping -- since we do not set any doping for the traps it is automatically zero
    params.doping[iphip, regionAcceptor] = Na

    ## boundary doping
    params.bDoping[iphip, bregionAcceptorRight] = Na
    params.bDoping[iphip, bregionAcceptorLeft ] = Na

    ## values for the schottky contacts
    params.SchottkyBarrier[bregionAcceptorLeft ] = barrier_left
    params.bVelocity[iphin, bregionAcceptorLeft] = vn
    params.bVelocity[iphip, bregionAcceptorLeft] = vp

    params.SchottkyBarrier[bregionAcceptorRight] = barrier_right
    params.bVelocity[iphin, bregionAcceptorRight] = vn
    params.bVelocity[iphip, bregionAcceptorRight] = vp


    data.params = params
    ctsys = System(grid, data, unknown_storage=unknown_storage)

    show_params(ctsys)
    println("*** done\n")

    ################################################################################
    println("Define outerior boundary conditions and enabled layers")
    ################################################################################

    ## set ohmic contact in bregionAcceptorRight and schottky contact in bregionAcceptorLeft
    set_contact!(ctsys, bregionAcceptorRight, Δu=0.0)
    set_contact!(ctsys, bregionAcceptorLeft, Δu=0.0)


    println("*** done\n")

    ################################################################################
    println("Define control parameters for Newton solver")
    ################################################################################

    control = NewtonControl()
    control.verbose = verbose
    control.damp_initial = 0.5
    control.damp_growth = 1.21    #>= 1
    control.max_iterations = 500
    control.tol_absolute = 1.0e-10
    control.tol_relative = 1.0e-10
    control.handle_exceptions = true
    control.tol_round = 1.0e-10
    control.max_round = 5

    control.Δt = timeStep_for_bias
    control.Δt_min = timeStep_for_bias/10
    control.Δt_max = tend / 10
    control.Δt_grow = 1.1

    println("*** done\n")


    ################################################################################
    println("Compute DLTS")
    ################################################################################

    for iTemp in eachindex(Temperatures)

        println(Temperatures[iTemp])
        println("\n")
        println("\n")

        println("Setup transient parameters and define system")
        data = Data(grid, numberOfCarriers, contactVoltageFunction=DLTScontactVoltageFunction)
        data.modelType = Transient
        data.F .= [FermiDiracOneHalfTeSCA, FermiDiracOneHalfTeSCA]
        data.bulkRecombination = set_bulk_recombination(; iphin=iphin, iphip=iphip,
            bulk_recomb_Auger=false,
            bulk_recomb_radiative=false,
            bulk_recomb_SRH=false)
        data.generationModel = GenerationNone
        data.boundaryType[bregionAcceptorLeft ] = SchottkyContact
        data.boundaryType[bregionAcceptorRight] = SchottkyContact#OhmicContact   
        data.fluxApproximation .= ExcessChemicalPotential

        params.temperature = Temperatures[iTemp]
        params.UT = (kB * params.temperature) / q

        data.params = params
        ctsys = System(grid, data, unknown_storage=unknown_storage)
        println("*** done\n")

        ################################################################################
        println("Compute solution in thermodynamic equilibrium")
        ################################################################################

        ## initialize solution and starting vectors
        initialGuess = unknowns(ctsys)
        equilibriumSolution = unknowns(ctsys)
        solution = unknowns(ctsys)

        #data.calculationType = inEquilibrium 

        ## solve thermodynamic equilibrium and update initial guess
        solution = equilibrium_solve!(ctsys, control=control, nonlinear_steps=20)
        equilibriumSolution .= solution
        initialGuess .= solution


        println("*** done\n")


        if plotting
            ## ##### set legend for plotting routines #####
            label_solution, label_density, label_energy = set_plotting_labels(data)

            ## for electrons 
            label_energy[1, iphin] = "\$E_c-q\\psi\$"
            label_energy[2, iphin] = "\$ - q \\varphi_n\$"
            label_density[iphin] = "n"
            label_solution[iphin] = "\$ \\varphi_n\$"

            ## for holes 
            label_energy[1, iphip] = "\$E_v-q\\psi\$"
            label_energy[2, iphip] = "\$ - q \\varphi_p\$"
            label_density[iphip] = "p"
            label_solution[iphip] = "\$ \\varphi_p\$"

            ## ##### set legend for plotting routines #####
            plot_energies(Plotter, ctsys, solution, "Energies in Equilibrium", label_energy)
            Plotter.figure()
            plot_densities(Plotter, ctsys, solution, "Densities in Equilibrium", label_density)
            Plotter.figure()
        end
        println("*** done\n")

        # ################################################################################
        # println("Compute Capacitance in thermodynamic equilibrium")
        # ################################################################################

        # τC = 1e-6 *s
        # ΔV = 0.001 * V
        # data.calculationType = OutOfEquilibrium
        # set_contact!(ctsys, bregionAcceptorLeft, Δu=ΔV)
        # solve!(solution, initialGuess, ctsys, control=control, tstep=τC)
        # charge_den_after_plus   = w_device * z_device * (charge_density(ctsys, solution)[regionAcceptor])
        # set_contact!(ctsys, bregionAcceptorLeft, Δu=-ΔV)
        # solve!(solution, initialGuess, ctsys, control=control, tstep=τC)
        # charge_den_after_minus  = w_device * z_device * (charge_density(ctsys, solution)[regionAcceptor])
        # println("Pojemność w stanie rownowagi w 0.0V")
        # println( abs(   charge_den_after_plus - charge_den_after_minus   ) * 1e9 / (2*ΔV) )

        # println("*** done\n")

        ################################################################################
        println("Automatic time stepping")
        ################################################################################
        data.calculationType = OutOfEquilibrium
        sol = solve(ctsys, inival=initialGuess, times=(0.0, tend), control=control)
        println(sol.t[:])
        # return
        if plotting
            tsol = sol(tend)
            plot_densities(Plotter, ctsys, tsol, "Densities at end time", label_density)
            Plotter.legend()
            Plotter.figure()
            plot_energies(Plotter, ctsys, tsol, "Energies at end time", label_energy)
            Plotter.figure()
            plot_solution(Plotter, ctsys, tsol, "Solution at end time", label_solution)
        end

        println("*** done\n")

        ################################################################################
        println("Postprocessing: Capacitance Calculation")
        ################################################################################

        ################################################################################
        println("Capacitance Calculation: Define System")
        ################################################################################

        data = Data(grid, numberOfCarriers)
        data.modelType = Transient
        data.F .= [FermiDiracOneHalfTeSCA, FermiDiracOneHalfTeSCA]
        data.bulkRecombination = set_bulk_recombination(; iphin=iphin, iphip=iphip,
                                                        bulk_recomb_Auger=false,
                                                        bulk_recomb_radiative=false,
                                                        bulk_recomb_SRH=false)
        data.generationModel = GenerationNone
        data.boundaryType[bregionAcceptorLeft ] = SchottkyContact
        data.boundaryType[bregionAcceptorRight] = SchottkyContact#OhmicContact   
        data.fluxApproximation .= ExcessChemicalPotential

        println("*** done CC:Define System\n")

        data.params                                         = params
        ctsys                                               = System(grid, data, unknown_storage=unknown_storage)

        ################################################################################
        println("Capacitance Calculation: Defining measurements  setup & calculating capacitance")
        ################################################################################

        tC_start = t1
        tC_end   = t2
        tC = range(tC_start, stop=tC_end, length=NumCapPoints)
        t_min = minimum(sol.t[:])
        t_max = maximum(sol.t[:])
        ΔV = 0.001 * V
        Capacitance = zeros(length(tC))

        #Why do we need those 3 lines?
        set_contact!(ctsys, bregionAcceptorRight, Δu = 0.0)
        set_contact!(ctsys, bregionAcceptorLeft,  Δu = 0.0)
        solution = equilibrium_solve!(ctsys, control=control, nonlinear_steps=20)

        data.calculationType = OutOfEquilibrium
        set_contact!(ctsys, bregionAcceptorLeft, Δu=ΔV)
        solve!(solution, initialGuess, ctsys, control=control, tstep=τC)
        charge_den_after_plus   = w_device * z_device * (charge_density(ctsys, solution)[regionAcceptor])
        set_contact!(ctsys, bregionAcceptorLeft, Δu=-ΔV)
        solve!(solution, initialGuess, ctsys, control=control, tstep=τC)
        charge_den_after_minus = w_device * z_device * (charge_density(ctsys, solution)[regionAcceptor])
        println("Pojemność w stanie rownowagi w 0")
        println( abs(   charge_den_after_plus - charge_den_after_minus   ) * 1e9 / (2*ΔV) )

        for iCapMeas in eachindex(tC)
            if tC[iCapMeas] >= t_min && tC[iCapMeas] <= t_max
                initialGuess .= sol(tC[iCapMeas])
                # charge_den_before = w_device * z_device * (charge_density(ctsys, initialGuess)[regionAcceptor])
                set_contact!(ctsys, bregionAcceptorLeft, Δu=voltageStep + ΔV)
                solve!(solution, initialGuess, ctsys, control=control, tstep=τC)
                charge_den_after_plus   = w_device * z_device * (charge_density(ctsys, solution)[regionAcceptor])
                set_contact!(ctsys, bregionAcceptorLeft, Δu=voltageStep - ΔV)
                solve!(solution, initialGuess, ctsys, control=control, tstep=τC)
                charge_den_after_minus = w_device * z_device * (charge_density(ctsys, solution)[regionAcceptor])
                Capacitance[iCapMeas] = abs(   charge_den_after_plus - charge_den_after_minus   ) / (2*ΔV)
                # (dQplus-dQminus)/dV = ((charge_den_after_plus - charge_den_before) - (charge_den_after_minus - charge_den_before))/dV
            end
        end
        #TODO: Loop over all Temperatures
        DLTS_values[iTemp] = (Capacitance[1]-Capacitance[2])*1e9
        # Plotter.figure()
        # plt.xscale("log")
        # Plotter.plot(tC, Capacitance*1e9)
        # plt.title("Capacitance transient DLTS")
        # plt.xlabel("t[s]")
        # plt.ylabel("C[nF]")

        println("*** done Capaciatance Calculations\n")
    end
    println(DLTS_values)
    Plotter.figure()
    Plotter.plot(Temperatures, DLTS_values)
    plt.title("Capacitance transient DLTS")
    plt.xlabel("T[K]")
    plt.ylabel("ΔC[nF]")

    # testval = solution[data.index_psi, 10]
    return 1#testval

end #  main

function test()
    testval = 1.3214196490674017
    main(test=true, unknown_storage=:dense) ≈ testval && main(test=true, unknown_storage=:sparse) ≈ testval
end

println("This message should show when this module has successfully recompiled.")


end # module




    # function scanProtocol(t)

    #     if 0.0 <= t && t <= tend_bias
    #         biasVal = voltageStep + ΔV* (t / tend_bias)
    #     elseif t > tend_bias && t <= tend
    #         biasVal = voltageStep + ΔV#+(voltageStep*t/tend)
    #     else
    #         biasVal = voltageStep
    #     end

    #     return biasVal

    # end
    # CapVolFun(t) = 0.0

    # contactVoltageFunction = [CapVolFun, zeroVoltage]
    # data.contactVoltageFunction = contactVoltageFunction 