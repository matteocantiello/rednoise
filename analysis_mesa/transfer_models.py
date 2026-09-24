"""
Simplified models for how the FeCZ could set the surface velocity field (v_macro), the
photometric amplitude (alpha0) and the characteristic frequency (nu_char).

Every candidate is a closed-form function of MESA history quantities with exponents fixed
by the physical argument; only a normalisation (and, in the tests, a floor) is fitted.
All inputs cgs unless the column name says otherwise. alpha_MLT = 1.6.

Notation (zone quantities of the FeCZ):
  v_c   zone-averaged convective velocity     v_x  maximum        v_t  averaged over the top alpha H_P
  rho_c zone-averaged density                 rho_t averaged over the top alpha H_P
  M_t   Mach number over the top alpha H_P    c_c  zone-averaged sound speed
  H_c   zone-averaged pressure scale height   dR   zone thickness   M_cz zone mass
  F_c/F max of L_conv/L in the zone           n_HP pressure scale heights between FeCZ top and surface
Surface: rho_s (outer cell), c_s (photosphere), H_s = c_s^2/(Gamma1 g), F* = sigma Teff^4.
"""
import numpy as np

ALPHA = 1.6
SIG = 5.6704e-5
G = 6.674e-8
RSUN, MSUN, DAY, LSUN = 6.957e10, 1.989e33, 86400.0, 3.828e33
GAMMA1 = 1.5   # only rescales H_s and P; enters as a constant factor


def inputs(df):
    """Physical inputs (cgs) from a mesa_ms/mesa_post dataframe."""
    q = {}
    q['R'] = 10**df.logR.values * RSUN
    q['g'] = 10**df.logg.values
    q['Teff'] = 10**df.logTeff.values
    q['L'] = 10**df.logL.values * LSUN
    q['M'] = df.M.values * MSUN
    q['Fstar'] = SIG * q['Teff']**4
    q['v_x'] = df.FeCZ_vmax_kms.values * 1e5
    q['v_c'] = df.FeCZ_vaver_kms.values * 1e5
    q['v_t'] = df.FeCZ_vahp_kms.values * 1e5
    q['rho_c'] = df.FeCZ_rho.values
    q['rho_t'] = df.FeCZ_rho_ahp.values
    q['M_t'] = df.FeCZ_mach_ahp.values
    q['M_x'] = df.FeCZ_mach_max.values
    q['c_c'] = df.FeCZ_cs_kms.values * 1e5
    q['H_c'] = df.FeCZ_hp_Rsun.values * RSUN
    q['dR'] = df.FeCZ_thick_R.values * q['R']
    q['M_cz'] = df.FeCZ_mass_Msun.values * MSUN
    q['FcF'] = df.FeCZ_FcF.values
    q['t_c'] = df.FeCZ_tc_s.values
    q['depth'] = df.FeCZ_depth_top_R.values * q['R']
    q['nHP'] = df.FeCZ_nHP_top.values
    q['rho_s'] = df.rho_surf.values
    q['c_s'] = df.cs_surf_kms.values * 1e5
    q['H_s'] = q['c_s']**2 / (GAMMA1 * q['g'])
    q['r_top'] = q['R'] - q['depth']
    q['P_t'] = q['rho_t'] * q['c_c']**2 / GAMMA1          # pressure at the FeCZ top (approx.)
    q['Sigma'] = q['P_t'] / q['g']                          # column mass above the FeCZ
    q['v_C09'] = df.FeCZ_vsurf_C09_kms.values * 1e5
    # core convection (control: core-excited waves)
    q['v_core'] = df.core_vmax_kms.values * 1e5
    q['rho_core'] = df.core_rho.values
    q['M_core'] = df.core_mach_max.values
    q['H_core'] = df.core_hp_Rsun.values * RSUN
    q['r_core'] = df.core_r_Rsun.values * RSUN
    return q


def wave_flux(q, kind):
    """Wave energy flux leaving the FeCZ top, diluted to the surface by (r_top/R)^2."""
    Fkin = q['rho_t'] * q['v_t']**3
    eff = {'GK': q['M_t'],                 # Goldreich & Kumar 1990: F_w ~ M F_kin (IGW, sharp interface)
           'LQ': q['M_t']**(5 / 8),        # Lecoanet & Quataert 2013, smooth interface
           'Li': q['M_t']**5}[kind]        # Lighthill acoustic emission
    return eff * Fkin * (q['r_top'] / q['R'])**2


def core_wave_flux(q):
    return q['M_core'] * q['rho_core'] * q['v_core']**3 * (q['r_core'] / q['R'])**2


def v_from_flux(F, q, kind):
    """Surface velocity carrying a wave energy flux F = 1/2 rho_s v^2 v_group."""
    if kind == 'cs':     # group velocity ~ surface sound speed
        return np.sqrt(2 * F / (q['rho_s'] * q['c_s']))
    if kind == 'self':   # group velocity ~ the velocity itself (turbulent cascade)
        return (2 * F / q['rho_s'])**(1 / 3)
    raise ValueError(kind)


def velocity_models(q):
    """Candidate predictors of v_macro [km/s], with the physical argument."""
    V = {}
    V['v_c,max'] = (q['v_x'], 'maximum MLT velocity in the FeCZ')
    V['v_c,aver'] = (q['v_c'], 'zone-averaged MLT velocity')
    V['v_c,top'] = (q['v_t'], 'MLT velocity averaged over the top alpha H_P')
    V['C09 microturb'] = (q['v_C09'], 'Cantiello+09: 1/2 rho_s v^2 = M_t 1/2 rho_t v_t^2')
    V['KE density'] = (q['v_t'] * np.sqrt(q['rho_t'] / q['rho_s']), 'kinetic-energy density conserved, no Mach factor')
    for fk, fl in (('GK', 'M F_kin'), ('LQ', 'M^5/8 F_kin'), ('Li', 'M^5 F_kin (acoustic)')):
        F = wave_flux(q, fk)
        V[f'wave {fk}, v_g=c_s'] = (v_from_flux(F, q, 'cs'), f'F_w = {fl} carried at c_s')
        V[f'wave {fk}, v_g=v'] = (v_from_flux(F, q, 'self'), f'F_w = {fl} carried at v')
    Fc = q['FcF'] * q['Fstar']
    for mk, M in (('M_x', q['M_x']), ('M_t', q['M_t'])):
        V[f'wave {mk} F_c, v_g=c_s'] = (v_from_flux(M * Fc * (q['r_top'] / q['R'])**2, q, 'cs'),
                                     f'F_w = {mk} x convective flux, carried at c_s')
        V[f'wave {mk} F_c, v_g=v'] = (v_from_flux(M * Fc * (q['r_top'] / q['R'])**2, q, 'self'),
                                   f'F_w = {mk} x convective flux, carried at v')
    V['v_c,max x sqrt(M_x)'] = (q['v_x'] * np.sqrt(q['M_x']), 'energy-density transfer with Mach efficiency, no density change')
    V['flux velocity'] = ((q['FcF'] * q['Fstar'] / q['rho_s'])**(1 / 3), 'v needed at rho_s to carry the FeCZ convective flux')
    V['Mach x c_s'] = (q['M_x'] * q['c_s'], 'FeCZ Mach number times surface sound speed')
    V['c_s (photosphere)'] = (q['c_s'], 'reference: photospheric sound speed alone')
    # saturation: whatever arrives with velocity v_in is capped near the photospheric sound speed,
    # v = (v_in^-2 + c_s^-2)^-1/2 (smooth minimum); the fitted normalisation k sets the plateau k c_s
    for vk in ('v_c,max', 'v_c,top', 'C09 microturb', 'KE density', 'wave GK, v_g=c_s', 'wave M_x F_c, v_g=v',
               'wave M_t F_c, v_g=v', 'wave M_x F_c, v_g=c_s', 'flux velocity'):
        vin = V[vk][0]
        V[f'sat[{vk}]'] = ((vin**-2 + q['c_s']**-2)**-0.5, f'{vk} saturated at the photospheric sound speed')
    V['core IGW, v_g=c_s'] = (v_from_flux(core_wave_flux(q), q, 'cs'), 'control: core-excited waves, M_core F_kin,core')
    return {k: (v / 1e5, d) for k, (v, d) in V.items()}


def amplitude_models(q, V):
    """Candidate predictors of alpha0 (fractional flux perturbation, dimensionless)."""
    A = {}
    A['F_c/F'] = (q['FcF'], 'convective flux fraction')
    A['F_kin/F*'] = (q['rho_t'] * q['v_t']**3 / q['Fstar'] * (q['r_top'] / q['R'])**2, 'kinetic energy flux fraction')
    for fk in ('GK', 'LQ', 'Li'):
        A[f'F_w({fk})/F*'] = (wave_flux(q, fk) / q['Fstar'], f'wave energy flux fraction, {fk}')
    cell = ALPHA * q['H_c'] / q['R']        # 1/sqrt(N_cells) for incoherent cells of size alpha H_c
    A['F_c/F x cells'] = (q['FcF'] * cell, 'F_c/F diluted by N_cells^-1/2, cells of size alpha H_c')
    A['F_w(GK)/F* x cells'] = (wave_flux(q, 'GK') / q['Fstar'] * cell, 'wave flux fraction diluted by N_cells^-1/2')
    for vk in ('C09 microturb', 'wave GK, v_g=c_s', 'wave GK, v_g=v', 'KE density'):
        Ms = V[vk][0] * 1e5 / q['c_s']
        A[f'M_s [{vk}]'] = (Ms, f'surface Mach number (linear density perturbation), v from {vk}')
        A[f'M_s^2 [{vk}]'] = (Ms**2, f'surface Mach^2 (pressure/turbulent-pressure perturbation), v from {vk}')
        A[f'M_s x cells [{vk}]'] = (Ms * cell, 'surface Mach diluted by N_cells^-1/2')
    for mk, M in (('M_x', q['M_x']), ('M_t', q['M_t'])):
        A[f'{mk} F_c/F'] = (M * q['FcF'], f'wave flux fraction F_w = {mk} F_c')
        A[f'{mk} F_c/F x cells'] = (M * q['FcF'] * cell, f'{mk} F_c/F diluted by N_cells^-1/2')
    A['F_c/F x cells(H_s)'] = (q['FcF'] * q['H_s'] / q['R'], 'F_c/F diluted by N_cells^-1/2, cells of size H_s at the surface')
    A['F_c/F x dR/H_c'] = (q['FcF'] * q['dR'] / q['H_c'], 'F_c/F x zone thickness in scale heights')
    A['F_c/F x cells x dR/H_c'] = (q['FcF'] * cell * q['dR'] / q['H_c'], 'cells-diluted F_c/F x zone thickness in H_P')
    A['F_c/F x depth/R'] = (q['FcF'] * q['depth'] / q['R'], 'F_c/F x depth of the FeCZ below the surface')
    A['M_cz/M_above'] = (q['M_cz'] / (4 * np.pi * q['R']**2 * q['Sigma']), 'convective mass over overlying mass')
    A['F_c/F x e^-nHP'] = (q['FcF'] * np.exp(-q['nHP']), 'F_c/F attenuated by overlying scale heights')
    A['F_c/F / nHP'] = (q['FcF'] / q['nHP'], 'F_c/F diluted linearly by depth in H_P')
    A['core IGW flux/F*'] = (core_wave_flux(q) / q['Fstar'], 'control: core-excited wave flux fraction')
    return A


def frequency_models(q, V):
    """Candidate predictors of nu_char [1/d]."""
    N = {}
    N['1/(2pi t_c)'] = (1 / (2 * np.pi * q['t_c']), 'FeCZ turnover, t_c = alpha H_c / v_c (Cantiello+21)')
    N['v_max/(2pi aH_c)'] = (q['v_x'] / (2 * np.pi * ALPHA * q['H_c']), 'turnover with the maximum velocity')
    N['v_top/(2pi aH_c)'] = (q['v_t'] / (2 * np.pi * ALPHA * q['H_c']), 'turnover of the top layer')
    N['v_c/(2pi dR)'] = (q['v_c'] / (2 * np.pi * q['dR']), 'crossing time of the zone thickness')
    N['v_top/(2pi depth)'] = (q['v_t'] / (2 * np.pi * q['depth']), 'crossing time of the overlying layer')
    for vk in ('C09 microturb', 'wave GK, v_g=c_s'):
        N[f'v_s/(2pi H_s) [{vk}]'] = (V[vk][0] * 1e5 / (2 * np.pi * q['H_s']), 'surface eddy turnover at the surface scale height')
    Eabove = 4 * np.pi * q['R']**2 * q['Sigma'] * 3 * q['P_t'] / q['rho_t']
    N['thermal, overlying'] = (q['L'] / (2 * np.pi * Eabove), 'inverse thermal time of the layer above the FeCZ')
    N['dynamical'] = (np.sqrt(G * q['M'] / q['R']**3) / (2 * np.pi), 'control: dynamical frequency')
    numax = 3090e-6 * (q['g'] / 27400) * (q['Teff'] / 5777)**-0.5
    N['nu_max scaling'] = (numax, 'control: solar-like nu_max ~ g Teff^-1/2')
    N['core turnover'] = (q['v_core'] / (2 * np.pi * ALPHA * q['H_core']), 'control: core convective turnover')
    return {k: (v * DAY, d) for k, (v, d) in N.items()}


def all_models(df):
    q = inputs(df)
    V = velocity_models(q)
    A = amplitude_models(q, V)
    N = frequency_models(q, V)
    return q, V, A, N
