"""Functions for calculating flux limits and neutrino counts"""

import numpy as np

from .neutrino import NeutrinoInteraction
from . import units



def veff_to_aeff(energies, effective_volume):
    """
    Convert effective volumes into effective areas.

    Parameters
    ----------
    energies : array_like
        Energies (GeV) at which effective volumes were calculated.
    effective_volume : array_like
        Effective volumes (m^3 sr) for the detector/station in question.

    Returns
    -------
    effective_area : ndarray
        Effective areas (m^2 sr) for the detector/station in question.

    """
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
    return np.asarray(effective_volume) * ice_density / int_len


def flux_sensitivity(energies, effective_area, livetime=units.yr, limit_factor=2.44):
    """
    Convert effective areas into sensitivities.

    Parameters
    ----------
    energies : array_like
        Energies (GeV) at which effective volumes were calculated.
    effective_area : array_like
        Effective areas (m^2 sr) for the detector/station in question.
    livetime : float, optional
        Total observation time (s) of the detector/station in question.
    limit_factor : float, optional
        Additional factor to multiply resulting sensitivity. See notes.

    Returns
    -------
    sensitivity : ndarray
        Sensitivities (GeV^-1 s^-1 m^-2 sr^-1) for the detector/station.

    Notes
    -----
    By default the limit factor is 2.44. This value is the factor which is used
    for calculating a Feldman-Cousins upper limit with zero background events.
    Other potentially useful values are 1 for a single-event-sensitivity
    calculation, and 2.3 for a Neyman upper limit with zero background events.

    """
    # Get number of energy bins per decade
    log_energy = np.log10(energies)
    d_log_energy = np.diff(log_energy)
    for d_log in d_log_energy:
        if not np.isclose(d_log, d_log_energy[0]):
            raise ValueError("Energies should be evenly spaced in log-10-space")
    bins_per_decade = 1/d_log_energy[0]

    factors = (limit_factor / livetime *
               bins_per_decade / np.log(10) / np.asarray(energies))

    return factors / np.asarray(effective_area)


def neutrino_count(model, energies, effective_area, livetime=units.yr, model_band=False):
    """
    Count the number of neutrinos observed for a given model at each energy.

    Parameters
    ----------
    model : Model
        A flux model object representing the desired source model.
    energies : array_like
        Energies (GeV) at which effective volumes were calculated.
    effective_area : array_like
        Effective areas (m^2 sr) for the detector/station in question.
    livetime : float, optional
        Total observation time (s) of the detector/station in question.
    model_band : bool, optional
        Whether to use the lower and upper bands of the model rather than
        the central flux value when calculating neutrino counts.

    Returns
    -------
    neutrino_counts : ndarray
        The number of neutrinos expected from the given model at each energy.
        If `model_band` was `False`, the shape of the output will match the
        shape of `energies`. If `model_band` was `True`, the shape of the
        output will be (N, 2), where N is the length of `energies`. The first
        column contains neutrino counts for the lower band of the model and the
        second column contains neutrino counts for the upper band.

    """
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

    sensitivities = flux_sensitivity(energies, effective_area,
                                     livetime=livetime, limit_factor=1)

    # Transposes are used to make sure the operation works for the band values
    # as well as simple flux values
    return (mean_fluxes.T / sensitivities).T
