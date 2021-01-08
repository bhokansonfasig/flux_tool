"""Functions for calculating flux limits and neutrino counts"""

import numpy as np

from .neutrino import NeutrinoInteraction
from . import units



def veff_to_aeff(energies, effective_volume):
    # Get average interaction lengths
    # (harmonic mean, since average should be in cross section)
    int_len = np.zeros(len(energies))
    for i, e in enumerate(energies):
        nu = NeutrinoInteraction(neutrino_type='nu_e',
                                 interaction_type=None,
                                 energy=e)
        nubar = NeutrinoInteraction(neutrino_type='nu_e_bar',
                                    interaction_type=None,
                                    energy=e)
        avg_int = 2/((1/nu.total_interaction_length)+
                     (1/nubar.total_interaction_length))
        int_len[i] = avg_int # water equivalent interaction length

    ice_density = 0.92 # relative to water density
    return effective_volume * ice_density / int_len


def flux_sensitivity(energies, effective_area, stations=1, livetime=1,
                     limit_factor=2.44):
    # Get number of energy bins per decade
    log_energy = np.log10(energies)
    d_log_energy = np.diff(log_energy)
    for d_log in d_log_energy:
        if not np.isclose(d_log, d_log_energy[0]):
            raise ValueError("Energies should be evenly spaced in log-10-space")
    bins_per_decade = 1/d_log_energy[0]

    factors = (limit_factor / livetime * stations *
               bins_per_decade / np.log(10) / energies)

    return factors / effective_area


def neutrino_count(energies, sensitivities, model, model_band=False):
    """Count the number of neutrinos observed for a given model at each energy"""
    log_energy = np.log10(energies)
    step = np.diff(log_energy)[0]

    if model_band:
        flux_min = lambda e: np.interp(e, model.energies, model.band_min)
        flux_max = lambda e: np.interp(e, model.energies, model.band_max)
        mean_fluxes = np.zeros((len(energies), 2))
    else:
        flux = lambda e: np.interp(e, model.energies, model.fluxes)
        mean_fluxes = np.zeros(len(energies))

    for i, log_e in enumerate(log_energy):
        e_range = np.logspace(log_e-step/2, log_e+step/2, 101)
        log_e_range = np.linspace(log_e-step/2, log_e+step/2, 101)
        if model_band:
            mean_fluxes[i, 0] = np.trapz(flux_min(e_range), x=log_e_range) / step
            mean_fluxes[i, 1] = np.trapz(flux_max(e_range), x=log_e_range) / step
        else:
            mean_fluxes[i] = np.trapz(flux(e_range), x=log_e_range) / step

    # Transposes are used to make sure the operation works for the band values
    # as well as simple flux values
    return (mean_fluxes.T / sensitivities).T
