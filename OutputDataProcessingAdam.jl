#Postprocessing of data (plotting tests)

#declaration of modules
using CSV 
using DataFrames
using Plots
using GLM
using Polynomials
using Unitful

#declaration of arrays and reading data form files
biasValues        = zeros(0) * u"V" #voltage
staticCapacitance = zeros(0) * u"F" #capacitance
filename = "staticCapacitance.csv" 
filepath = joinpath(@__DIR__, filename)
filename2 = "biasValues.csv" 
filepath2 = joinpath(@__DIR__, filename2)
staticCapacitance        = CSV.read(filepath  , DataFrame; header=false) 
biasValues               = CSV.read(filepath2 , DataFrame; header=false) 

#creating data for plots
y    = abs.(staticCapacitance.Column1).^-2
bias = biasValues.Column1#.+biasValues.Column1[2]/2
bias = bias[1:100]
bias = abs.(bias)

data = DataFrame(X=bias[1:95],Y=y[1:95])
ols = lm(@formula(Y ~ X), data)
coef(ols)

function display_capacitance() #function to display plots

    plot1=plot(biasValues.Column1[1:95],abs.(staticCapacitance.Column1[1:95])*10^9, reuse=false, xlabel = "U[V]", ylabel="C[nF]", seriestype = :scatter)
    plot2=plot(bias[1:95],y[1:95]*10^(-18), reuse=false, xlabel = "U[V]", ylabel="1/C[nF]^2", seriestype = :scatter)
    

    model(x) = coef(ols)[1] + coef(ols)[2] * x
    plot3=plot(bias[1:95],y[1:95]*10^(-18), reuse=false, xlabel = "U[V]", ylabel="1/C[nF]^2", seriestype = :scatter, title="Linear fit:")
    plot!(data.X, model.(data.X)*10^(-18), lw=3)

    plot(plot1, plot2, plot3, layout = (3,1), legend=false) 
    
end
 function cal_doping_and_build_in_pot() #function for calculations 
    eps0       = 8.85*10^-14 * u"F/cm"
    epsRelativ = 13.6 
    S          = 0.5*0.5
    q          = 1.602176565*10^-19 * u"C"
    N_derived = 2/(coef(ols)[2]*q*epsRelativ*eps0*S^2)
    build_in_potential = coef(ols)[1]/coef(ols)[2]
    return N_derived, build_in_potential
 end

 function print_liner_fit_coefs() #test function for printing values
    printer = coef(ols)
    return printer
 end

