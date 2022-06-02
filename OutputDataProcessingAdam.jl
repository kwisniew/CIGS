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
bias = biasValues.Column1.+biasValues.Column1[2]/2
bias = bias[1:100]
bias = abs.(bias)

data = DataFrame(X=bias,Y=y)
ols = lm(@formula(Y ~ X), data)
coef(ols)

function display() #function to display plots

    plot1=plot(bias,y, reuse=false, xlabel = "Voltage", ylabel="Capacitance")

    model(x) = coef(ols)[1] + coef(ols)[2] * x
    #1.4775483873006742e7
    build_in_potential = coef(ols)[1]/coef(ols)[2]
    scatter(data.X, data.Y)
    plot2= plot(data.X, model.(data.X), legend=false, reuse = false, xlabel = "Stranger", ylabel="Things")

    plot(plot1, plot2, layout = (2,1), legend=false) 
    
end
 function calculate() #function for calculations 
    eps0       = 8.85*10^-14 * u"F/cm"
    epsRelativ = 13.6 
    S          = 0.5*0.5
    q          = 1.602176565*10^-19 * u"C"
    N_derived = 2/(coef(ols)[2]*q*epsRelativ*eps0*S^2)
    return N_derived
 end

 function print() #test function for printing values
    printer = coef(ols)
    return printer
 end

