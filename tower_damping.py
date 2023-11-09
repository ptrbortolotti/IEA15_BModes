import numpy as np

# Frequency 1st FA, 2nd FA, 1st SS, 2nd SS
freqs_rna_off = np.array([0.77188E+00, 0.32533E+01, 0.77188E+00, 0.32533E+01])
freqs_rna_on = np.array([0.24655E+00, 0.11238E+01, 0.24368E+00, 0.73358E+00])

# Target log dec, natural log of the ratio of the amplitudes
#  of any two successive peaks. 
delta = np.array([0.06, 0.06, 0.06, 0.06]) 
# Convert to critical damping
# Dimensionless measure describing how oscillations 
# in a system decay after a disturbance
zeta = 1. / np.sqrt(1.+(2.*np.pi / delta)**2.)
# zeta = 0.32 * np.ones_like(delta)
omega_rna_off = freqs_rna_off*2*np.pi # (rad/s)
omega_rna_on = freqs_rna_on*2*np.pi # (rad/s)

mu_rna_off = 2.*zeta/omega_rna_off
mu_rna_on = 2.*zeta/omega_rna_on

# target_twr_damp = np.array([0.32, 0.32, 0.32, 0.32])

input_twr_damp = zeta * freqs_rna_off / freqs_rna_on
print(input_twr_damp)
delta = 2*np.pi/np.sqrt((1/zeta)**2.-1)

# Tower only HAWCStab2
delta = 0.02 # log dec
zeta = 1. / np.sqrt(1.+(2.*np.pi / delta)**2.)
mu_rna_off = 2.*zeta/omega_rna_off
