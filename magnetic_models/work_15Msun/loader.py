import mesa_reader as mr
import numpy as np
import matplotlib.pyplot as plt

data = mr.MesaData(file_name='./LOGS/history.data')
index = mr.MesaProfileIndex(file_name='./LOGS/profiles.index')

last_model = index.model_numbers[-1]
last_model_name = index.profile_with_model_number(last_model)

prof = mr.MesaData(file_name='./LOGS/profile' + str(last_model_name) + '.data')

Rsun = 7e10

opacity = prof.data('opacity')
R = Rsun * 10**prof.data('logR')
D = 10**prof.data('logRho')
T = 10**prof.data('logT')
N2 = prof.data('brunt_N2')
N2_thermal = prof.data('brunt_N2_composition_term')

R_with_zero = np.array(list(R) + [0])
dR = np.diff(R_with_zero)

sel = ((N2 < 0) & (R < 0.2 * R[0]))
N2_conv_avg = -np.sum(N2[sel] * dR[sel]) / np.sum(dR[sel])
N_conv_avg = N2_conv_avg**0.5

l_grid = 1 + np.arange(20)
f_grid = N_conv_avg * 10**np.linspace(-2,2,num=20,endpoint=True)

def make_plot(l,f):

	k_perp = (l * (l + 1))**0.5 / R
	k_r = k_perp * (N2**0.5 / f)

	vA_crit_radial = f / k_r
	B_crit_wave_radial = vA_crit_radial * (4 * np.pi * D)**0.5
	vA_crit_toroidal = f / k_perp
	B_crit_wave_toroidal = vA_crit_toroidal * (4 * np.pi * D)**0.5


	plt.plot(R / R[0], B_crit_wave_radial, label='$B_{\mathrm{crit,wave,radial}} / \mathrm{G}$', linestyle='--')
	plt.plot(R / R[0], B_crit_wave_toroidal, label='$B_{\mathrm{crit,wave,toroidal}} / \mathrm{G}$')
	plt.plot(R / R[0], vA_crit_radial, label='$v_{A,\mathrm{crit,radial}} / \mathrm{G}$', linestyle='--')
	plt.plot(R / R[0], N2**0.5, label='$N / \mathrm{s^{-1}}$')
	plt.plot(R / R[0], R*k_r, label='$r k_r$')
	plt.axhline(f, label='$\omega$',linestyle='-.')
	plt.title('l = ' + str(l) + ' ' + 'f = ' + str(round(f / N_conv_avg,2)) + ' f_conv')
	plt.legend()
	plt.xlabel('r / R')
	plt.yscale('log')
	plt.savefig(str(l) + '_' + str(round(f / N_conv_avg,2)) + '.pdf')
	plt.close()


def make_composite_plot(ls,fs,ftitles):

	plt.figure(figsize=(5,4))

	for l,f,ftitle in zip(*(ls,fs,ftitles)):

		k_perp = (l * (l + 1))**0.5 / R
		k_r = k_perp * (N2**0.5 / f)

		vA_crit_radial = f / k_r
		B_crit_wave_radial = vA_crit_radial * (4 * np.pi * D)**0.5

		plt.plot(1-R / R[0], B_crit_wave_radial, label='$\ell = ' + str(l) + ', \\nu =  ' + ftitle + '$')
	plt.legend()
	plt.xlim([1e-5,1e-2])
	plt.ylim([1e-2,1e5])
	plt.xscale('log')
	plt.xlabel('$1 - r / R_\star$')
	plt.yscale('log')
	plt.ylabel('$B_{r,\mathrm{crit}} / \mathrm{G}$')
	plt.tight_layout()
	plt.savefig('composite.pdf')
	plt.close()

make_composite_plot([1,2,3,1,2,3], 3*[N_conv_avg] + 3*[2*np.pi*5/(24*3600)], 3*['\\nu_\mathrm{conv}'] + 3*['\\nu_\mathrm{char}'])