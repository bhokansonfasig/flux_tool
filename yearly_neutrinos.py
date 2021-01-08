from collections import namedtuple
import numpy as np
import matplotlib.pyplot as plt

import flux_tool
from limit_figure import LimitFigure


def plot_neutrinos_by_year_ara(models, energies, veff_yearly, styles=None,
                               years=[2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022],
                               stations=["ARA 1", "ARA 1-3", "ARA 1-3", "ARA 1-3", "ARA 1-3", "ARA 1-3", "ARA 1-5"],
                               xlim=None, save_name=None, color_damping=0):
    plt.figure(figsize=(7, 6))
#     plt.axhline(5, ls=':', c='k')
#     plt.grid(ls='--', c='dimgray', axis='y')

    aeff_yearly = flux_tool.veff_to_aeff(energies, veff_yearly)
    sensitivity_yearly = flux_tool.flux_sensitivity(energies, aeff_yearly, limit_factor=1)

    detections = {}
    for name, model in models.items():
        detections[name] = [0]
        for i, year in enumerate(years[1:]):
            nus = flux_tool.neutrino_count(energies, sensitivity_yearly[i], model[0], model_band=bool(model[1]))
            if model[1]:
                nus = nus[:, model[1]-1]
            detections[name].append(np.sum(nus))

    uniq = []
    for name in stations:
        if name not in uniq:
            uniq.append(name)
    i = 0
    c = 0
    while i<min(len(stations), len(years)):
        xmin = years[i]
        st = stations[i]
        j = i+1
        while j<len(stations) and stations[j]==st:
            j += 1
        if j==len(stations):
            xmax = years[-1]
        else:
            xmax = years[j]
        cval = 1-c/(len(uniq)+color_damping)
        plt.axvspan(xmin, xmax, color=(cval, cval, cval, 0.15))
        plt.annotate(str(stations[i]),
                     xy=(i/(len(years)-1)+0.1/(len(years)-1), 0.4), xycoords='axes fraction',
                     color='k', rotation=90, verticalalignment='bottom', fontsize=12)
        i = j
        c += 1

    for model, counts in detections.items():
        if styles is None:
            style = '.-'
            color = None
            ls = None
        else:
            style = styles[model]
            color = style.rstrip('.-:&')
            style = style[len(color):]
            ls = None
            if style.endswith('&'):
                style = style[:-1]
                ls = (0, (5, 10))
        line, = plt.plot(years, np.cumsum(counts), style, color=color, label=model)
        if ls is not None:
            line.set_linestyle(ls)
    plt.xlabel("Year", fontsize=16)
    plt.ylabel("Expected Number of Neutrinos", fontsize=16)
    if xlim is None:
        plt.xlim(years[0], years[-1])
    else:
        plt.xlim(*xlim)
    # plt.ylim(0, plt.ylim()[1])
    plt.ylim(0, 4)
    plt.axvline(2019, c='k', ls='--')
    ticks = plt.xticks()[0]
    plt.xticks(ticks[ticks%1==0], [str(int(tick)) for tick in ticks[ticks%1==0]])
    plt.legend(loc=2, fontsize=12)
    ax1 = plt.gca()
    ax2 = ax1.twinx()
    ax2.set_ylim(ax1.get_ylim())
#     ax2.set_yticklabels([])
    for item in (ax1.get_xticklabels() + ax1.get_yticklabels() + ax2.get_yticklabels()):
        item.set_fontsize(14)
    plt.tight_layout()
    if save_name is not None:
        plt.savefig(save_name, dpi=300)
    plt.show()




if __name__=="__main__":
    ara_energies = 10**np.array([16.5, 17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0]) * flux_tool.units.eV
    ara_100m = np.array([8.818E-2, 4.006E-1, 1.381E+0, 3.143E+0, 6.155E+0, 1.045E+1, 1.565E+1, 2.248E+1]) * flux_tool.units.km**3*flux_tool.units.sr
    ara_200m = np.array([1.105E-1, 5.195E-1, 1.826E+0, 5.259E+0, 1.106E+1, 2.077E+1, 3.431E+1, 5.073E+1]) * flux_tool.units.km**3*flux_tool.units.sr
    ara_200mpa = np.array([2.618E-1, 1.003E+0, 3.096E+0, 7.406E+0, 1.466E+1, 2.578E+1, 4.039E+1, 5.073E+1]) * flux_tool.units.km**3*flux_tool.units.sr

    ara_by_year = np.array([284*ara_100m/366, (217+233)*ara_200m/365, (124*ara_100m+(313+304)*ara_200m)/365,
                           (345+250)*ara_200m/365, (127*ara_100m+(323+297)*ara_200m)/366, 137*ara_200m/365,
                           (316*ara_100m+(303+316+313)*ara_200m+186*ara_200mpa)/365,
                           ara_100m+3*ara_200m+ara_200mpa, ara_100m+3*ara_200m+ara_200mpa, ara_100m+3*ara_200m+ara_200mpa]) * flux_tool.units.yr

    # Treat the i3_nu_fit function as a model object
    FauxModel = namedtuple('FauxModel', ['energies', 'fluxes'])
    i3_nu_fit_model = FauxModel(ara_energies, LimitFigure._i3_nu_fit(ara_energies, slope=-2.19, offset=1.01))

    # Models dictionary structure - "model name": (model_object, column_int)
    # Where column_int is 0 for a pure flux, 1 for a minimum bound,
    # and 2 for a maximum bound
    models = {
        "SFR $E_{max}=10^{21.5}$, Kotera et al.": (flux_tool.Model('kotera_high_e'), 0),
        "SFR $E_{max}=10^{20.5}$, Kotera et al.": (flux_tool.Model('kotera_mid'), 0),
        "IceCube $E^{-2.19}$ Extrapolated": (i3_nu_fit_model, 0),
        "100% protons, Ahlers & Halzen": (flux_tool.Model('ahlers_100'), 0),
        "10% protons, Ahlers & Halzen": (flux_tool.Model('ahlers_10'), 0),
        "1% protons, Ahlers & Halzen": (flux_tool.Model('ahlers_1'), 0),
    }

    styles = {
        "SFR $E_{max}=10^{21.5}$, Kotera et al.": "deeppink.--",
        "SFR $E_{max}=10^{20.5}$, Kotera et al.": "darkmagenta.--",
        "IceCube $E^{-2.19}$ Extrapolated": "dodgerblue.-",
        "100% protons, Ahlers & Halzen": "mediumblue.-.",
        "10% protons, Ahlers & Halzen": "royalblue.-.",
        "1% protons, Ahlers & Halzen": "cornflowerblue.-.",
    }

    plot_neutrinos_by_year_ara(models, ara_energies, ara_by_year, styles=styles,
                               years=[2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022],
                               stations=["ARA 1", "ARA 1-3", "ARA 1-3", "ARA 1-3", "ARA 1-3", "ARA 1-3", "ARA 1-5"],
                               save_name='yearly_neutrinos_ara_fluxes.pdf')



    models = {
        "AGN, Murase": (flux_tool.Model('murase_agn'), 0),
        "GRB afterglow-late prompt, Murase": (flux_tool.Model('murase_grb_late_prompt'), 0),
        "GRB afterglow-wind, Murase": (flux_tool.Model('murase_grb_wind'), 0),
        "GRB afterglow-ISM, Murase": (flux_tool.Model('murase_grb_ism'), 0),
        "Clusters, Fang & Murase": (flux_tool.Model('fang_cluster'), 0),
    }

    styles = {
        "AGN, Murase": "red.--",
        "GRB afterglow-late prompt, Murase": "saddlebrown.-.",
        "GRB afterglow-wind, Murase": "goldenrod.-.",
        "GRB afterglow-ISM, Murase": "gold.-.",
        "Clusters, Fang & Murase": "mediumvioletred.&",
    }

    plot_neutrinos_by_year_ara(models, ara_energies, ara_by_year, styles=styles,
                               years=[2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022],
                               stations=["ARA 1", "ARA 1-3", "ARA 1-3", "ARA 1-3", "ARA 1-3", "ARA 1-3", "ARA 1-5"],
                               save_name='yearly_neutrinos_ara_sources.pdf')