import numpy as np
import matplotlib.pyplot as plt

dirname = "./"

plt.rcParams.update({'font.size': 22})
plt.rcParams.update({'legend.fontsize': 22})

fig, ax = plt.subplots(figsize=(2,6))


import os

dat1_fullname     = os.path.join( dirname, 'data/Gw_1_1_i7.dat')
print( dat1_fullname )
dat2_fullname     = os.path.join( dirname, 'data/Gw_2_2_i7.dat')
print( dat2_fullname )

x1_freq, g1_re, g1_im = np.genfromtxt( dat1_fullname ).transpose()
x2_freq, g2_re, g2_im = np.genfromtxt( dat2_fullname ).transpose()

plt.fill( -g1_im, x1_freq, "b-", label="DMFT"+r"($\sigma=\downarrow$)", alpha=0.6 )
plt.fill( -g2_im, x2_freq, "r-", label="DMFT"+r"($\sigma=\uparrow$)"  , alpha=0.6 )

plt.plot( -g1_im, x1_freq, "b-", alpha=0.6 )
plt.plot( -g2_im, x2_freq, "r-", alpha=0.6 )

plt.legend(fontsize=8)

# ax.set_axisbelow(True)
print( "axis : ", ax.axis() )


# plt.grid(axis="y",color="grey",linestyle=":")
plt.subplots_adjust( left=0.25 )
plt.savefig("fig_lco_dos_U6.pdf")

plt.show()
