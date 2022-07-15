#File to extract doping and build in potential from capacitance - voltage relationship
#it is also used for postprocessing of data (plotting)

#declaration of modules
using CSV 
using DataFrames
using Plots
#using PyPlot
using GLM
using Polynomials
using Unitful
using NumericIO


#declaration of arrays and reading data form files
biasValues        = zeros(0) #voltage
staticCapacitance = zeros(0) #capacitance
filename = "staticCapacitance.csv" 
filepath = joinpath(@__DIR__, filename)
filename2 = "biasValues.csv" 
filepath2 = joinpath(@__DIR__, filename2)
staticCapacitance        = CSV.read(filepath  , DataFrame; header=false) 
biasValues               = CSV.read(filepath2 , DataFrame; header=false) 

#One can get from capacitance-voltage (CV) relationship the values of "build in potential" and "doping density"
#For details see eq. 3.10 in https://in.ncu.edu.tw/ncume_ee/SchottkyDiode.htm  
#We will try to calculate "build in potential" and "doping density" from linear fit after plotting (1/C^2) vs. V
y    = abs.(staticCapacitance.Column1).^-2
#The capacitance is calculated by numerical differentiation of the uncompensated charge density with respect to voltage C(V)=dQ(V)/dV. 
#One can argue that C calculated this way should be evaluated not in V, but V+dV/2: C(V+dV/2)= dQ(V)/dV
bias = biasValues.Column1[1:(length(biasValues.Column1)-1)] .+ diff(biasValues.Column1/2)
#make sure that X and Y have the same length 
bias = bias[1:length(y)]
#bias = abs.(bias)

#The last few values of C deviates from other, so we artificially cut of the last five of it
cutoff = 0
#TO DO: Check what cause the deviation - this is potential bug. 
data = DataFrame(X=bias[1:(length(bias)-cutoff)],Y=y[1:(length(y)-cutoff)])

#linear fit
ols = lm(@formula(Y ~ X), data)
intercept = coef(ols)[1]* u"1/F^2"
slope = coef(ols)[2]* u"(1/F^2)/V"

#function to display plots
function display_capacitance() 

   #Capacitance vs Voltage Plot
   plot1=plot(biasValues.Column1[1:(length(biasValues.Column1)-cutoff)],abs.(staticCapacitance.Column1[1:(length(staticCapacitance.Column1)-cutoff)])*10^9,
              reuse=false, xlabel = "U[V]", ylabel="C[nF]", seriestype = :scatter)
   #1/C^2 vs Voltage Plot
   plot2=plot(bias[1:(length(bias)-cutoff)],y[1:(length(y)-cutoff)]*10^(-18), 
              reuse=false, xlabel = "U[V]", ylabel="1/C[nF]^2", seriestype = :scatter)
   
   #linear function based on linear fit
   model(x) = intercept + slope * (x* u"V") 
   #1/C^2 vs Voltage Plot with a line on the top of it 
   plot3=plot(bias[1:(length(bias)-cutoff)],y[1:(length(y)-cutoff)]*10^(-18), 
              reuse=false, xlabel = "U[V]", ylabel="1/C[nF]^2", seriestype = :scatter, title="Linear fit:")
   plot!(data.X, model.(data.X)*10^(-18)* u"F^2", lw=3)

   plot(plot1, plot2, plot3, layout = (3,1), legend=false) 
    
end
#function for extracting doping and build in potential
function cal_doping_and_build_in_pot()  
   eps0       = 8.85*10^-14 * u"F/cm"
   epsRelativ = 13.6 
   #ASSUMPTION: width and depth = 1cm. Check if this parameters are the same as in the script you use to get "biasValues" and "staticCapacitance"!
   #TO DO: read  width and depth from the file.
   w_device = 1.0 *u"cm"  # width of device
   z_device = 1.0 *u"μm"#cm"  # depth of device
   S        = w_device*z_device #1cm2
   q          = 1.602176565*10^-19 * u"C"
   N_derived = 2/(slope*q*epsRelativ*eps0*S^2)
   build_in_potential = intercept/slope
   return formatted(ustrip(u"cm^-3",(-N_derived)), :SCI, ndigits=3)*" cm^-3", round(typeof(1.0u"V"),-build_in_potential;digits=2)
end
#upreferred(N_derived)
#test function for printing values
function print_linear_fit_coefs() 
   printer = coef(ols)
   return printer
end

