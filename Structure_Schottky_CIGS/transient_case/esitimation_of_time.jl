using ChargeTransport
using Unitful

mup_CIGS = 25* u"(cm^2) / (V * s)"
Na       = 1.0e15 *u"1/(cm^3)"

#estimation of depletion region (you can read it from plot od energy band in equilibrium)
W = 0.75 * u"μm"

w_device = 1.0 *u"cm"  # width of device
z_device = 1.0 *u"cm"  # depth of device
S        = w_device*z_device #1cm2

Q = 1.602176565*10^-19 * u"C"
#Resistivity
ρ = 1/(Q*Na*mup_CIGS)
#Resistance (from Pouillet's law)
R= ρ * W / S

#I'm approximating the value of C from the CV plot: C was taken for V=0
C = 1.1e-8 *u"F"

τ = R*C
upreferred(τ)