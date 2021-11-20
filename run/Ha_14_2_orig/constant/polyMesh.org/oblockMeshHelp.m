# Càlculs per trobar el mallat


Ha  = 14.2      # Hartmann al Manifold

display(" ")
display("CAS benchmarking per al Manifold")
display("==========================")

### CARACTERÍSTIQUES FÍSIQUES ====================================================
T	= 573.15;			            # [K] són 300ºC ve definit per disseny TEMPERATURA
k	= 1.95+0.0196*T;		        # [W/(m·K)] CONDUCTIVITAT TÈRMICA
Cp	= 195-0.009116*T;		        # [J/(Kg·K)] COEFICIENT CALORÍFIC
rho	= 10520*(1-0.000113*T)		    # [kg/m3] 9838,662206 kg/m3 a 300ºC DENSITAT
mu	= 0.000187*exp(11640/(8.314*T));# [Pa·s]    VISCOSITAT
sigma	= 100/(0.0001023+4.26e-08*T);    #   [UNITATS???]  CONDUCTIVITAT ELÈCTRICA	
nu	= mu/rho;
Pr	= mu*Cp/k;
### DADES DIBUIX: =================================================================
a   =   1;
b   =   1;
w   =   0.2;
xin =   a;
xout=   3*a;
### Condicions  =========================================================
display(" ")
display("Velocity")
display("========")
Re  =   200
L  =    w;       # characteristic length for obstacles is Cylinder width.
U   =   Re*nu/L

display(" ")
display("Magnetic Field")
display("===============")
Ha
B   =   Ha/(L*(sigma/(nu*rho))^0.5)
N   =   Ha^2/Re
B_ha_6_33   = 6.33/(L*(sigma/(nu*rho))^0.5)
N_ha_6_33   = 6.33^2/Re

ReHa        = Re/Ha
ReSqHa      = Re/sqrt(Ha)
if (ReHa>200)
    display(" ")
    display("Regim Real és Turbulent")
elseif (ReHa<200 && ReSqHa>65)
    display(" ")
    display("Regim Real és Q2D")
elseif (ReSqHa<65)
    display(" ")
    display("Regim Real és Laminar")
endif



display("================================")
display(" ")
display("Per un ...")
xlen        = (xin+xout);
tau         = 2*rho/(sigma*B_ha_6_33^2)   # ha de ser major que els timesteps ref. Tesis Elisabet p67.
# prové del Von Newman stability analysis'
display(" ")
display(".. el CFL ha de ser com a màxim:")
CFL_max     = tau*U/xlen
# desitjat, en funció de si es 1-D (CFL=1), 2D (CFL=0.5) o 3D (CFL=0.2).
#display(" ")
Nvoltes		= 10
endTime		= Nvoltes*xlen/U













