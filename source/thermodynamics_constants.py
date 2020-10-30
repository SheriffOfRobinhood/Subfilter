# Last updated 29/10/2020 by Peter Clark
from numpy import pi

MONC = True

if MONC : g       =              9.81
else :    g       =              9.80665
grav = g
avogadro          =              6.02214076e23
k_boltzmann       =              1.380649e-23
planck            =              6.62607015e-34
speed_of_light    =              2.99792458e8
perf_gas_const    = k_boltzmann*avogadro
sigma_stephan     =              5.67000e-08
solar_const       =           1365.00
radius_earth      =              6.37123e06
cp_water          =           4218.00
cfus_water        =         333500.
if MONC : cvap_water        =    2.501E6
else :    cvap_water        =    2.48890e+06
mol_wt_water      =             18.0150
gas_const_water   = perf_gas_const / mol_wt_water * 1000
cp_ice            =           2106.00
rho_ice           =            917.000
freeze_pt         =            273.150
triple_pt         =            273.160
boil_pt           =            373.150
gas_vol           =              0.0224100
p_ref             =         101325.0
p_ref_um          =         100000.0
mol_wt_air        =             28.9640
gas_const_air     = perf_gas_const / mol_wt_air * 1000
cp_air            =           1005.00
cv_air            =            718.000
cp_by_cv          =           1.40000
rho_air_stp       =              1.29300
viscosity_air_stp =              1.73000e-05
thermal_cond_air  =              0.0240000
cp_water_vap      =           1875.0
secs_day          =       86400.0

epsilon           = mol_wt_water/mol_wt_air
#  originally dk=0.28557214
kappa             = gas_const_air/cp_air
rk                = 1./kappa
kappa_v = 0.28
ts                = freeze_pt
# p0                = 1000.0
p0                = p_ref_um/100.0
rho0              = rho_air_stp
kttoms            = radius_earth*2*pi/360.0/60.0/3600.0
# originally 1.2923
rw                = gas_const_air/epsilon
rhow              = 1000.0
c_virtual         = 1.0/epsilon-1.0
L_over_cp         = cvap_water / cp_air
