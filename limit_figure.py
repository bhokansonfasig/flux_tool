"""Class for creating limit/sensitivity figure"""

import numpy as np
import matplotlib.pyplot as plt

import flux_tool


DEFAULT_ENERGY_UNITS = flux_tool.units.GeV
DEFAULT_FLUX_UNITS = (flux_tool.units.GeV**-1 * flux_tool.units.cm**-2 *
                      flux_tool.units.s**-1 * flux_tool.units.sr**-1)


# This class based on Anna Nelles's plotting script:
# https://github.com/nu-radio/NuRadioMC/blob/138f8419e2db935bd07cb41d88ff2ea1b9ee99e1/NuRadioMC/examples/Sensitivities/E2_fluxes2.py
class LimitFigure:
    def __init__(self, figsize=(7, 6), xlims=(1e5, 1e11), ylims=(1e-11, 1e-5),
                 energy_units=DEFAULT_ENERGY_UNITS, flux_units=DEFAULT_FLUX_UNITS,
                 e_bins_per_decade=1, e_power=2, font_size=12, tick_size=12):
        self.fig, self.ax = plt.subplots(1, 1, figsize=figsize)

        self.ax.set_xscale('log')
        self.ax.set_yscale('log')

        self.ax.set_xlabel(r'Neutrino Energy [GeV]')
        if e_power==2:
            self.ax.set_ylabel(r'$E^2\Phi$ [GeV cm$^{-2}$ s$^{-1}$ sr$^{-1}$]')
        elif e_power==1:
#             self.ax.set_ylabel(r'$E\Phi$ [cm$^{-2}$ s$^{-1}$ sr$^{-1}$]')
            self.ax.set_ylabel(r'$E\ dN/(dE\ dA\ d\Omega\ dt)$ [cm$^{-2}$ s$^{-1}$ sr$^{-1}$]')
        else:
            raise ValueError(f"Invalid power ({e_power})")

        self.ax.set_xlim(*xlims)
        self.ax.set_ylim(*ylims)

        self.ax.title.set_fontsize(font_size)

        plt.tight_layout()

        self.ax.xaxis.label.set_fontsize(font_size)
        self.ax.yaxis.label.set_fontsize(font_size)
        for label in self.ax.get_xticklabels() + self.ax.get_yticklabels():
            label.set_fontsize(tick_size)

        self.e_power = e_power
        self.e_unit = energy_units
        self.f_unit = flux_units * energy_units**e_power
        self.e_bins = e_bins_per_decade
        self.font_size = font_size
        
        self.neutrino_models = []
        self.custom_limits = []

    @staticmethod
    def _i3_nu_fit(energy, slope=-2.13, offset=0.9):
        flux = 3 * offset * (energy / (100 * flux_tool.units.TeV))**slope * 1e-18
        return flux * (flux_tool.units.GeV**-1 * flux_tool.units.cm**-2 *
                       flux_tool.units.s**-1 * flux_tool.units.sr**-1)

    @classmethod
    def _get_i3_mu_range(cls):
        energy = np.arange(1e5, 5e6, 1e5) * flux_tool.units.GeV
        upper = np.maximum(cls._i3_nu_fit(energy, offset=0.9, slope=-2.),
                           cls._i3_nu_fit(energy, offset=1.2, slope=-2.13))
        lower = np.minimum(cls._i3_nu_fit(energy, offset=0.9, slope=-2.26),
                           cls._i3_nu_fit(energy, offset=0.63, slope=-2.13))
        return energy, upper, lower

    @classmethod
    def _get_i3_hese_range(cls):
        energy = np.arange(1e5, 5e6, 1e5) * flux_tool.units.GeV
        upper = np.maximum(cls._i3_nu_fit(energy, offset=2.46, slope=-2.63),
                           cls._i3_nu_fit(energy, offset=2.76, slope=-2.92))
        lower = np.minimum(cls._i3_nu_fit(energy, offset=2.46, slope=-3.25),
                           cls._i3_nu_fit(energy, offset=2.16, slope=-2.92))
        return energy, upper, lower


    def get_data(self, filename):
        model = flux_tool.Model(filename)
        energies = model.energies / self.e_unit
        fluxes = model.fluxes * model.energies**self.e_power / self.f_unit
        band_min = model.band_min * model.energies**self.e_power / self.f_unit
        band_max = model.band_max * model.energies**self.e_power / self.f_unit
        return energies, fluxes, band_min, band_max


    def add_model(self, name):
        if name=='heinze':
            energy, flux, band_min, band_max = self.get_data('heinze_cr')
            heinze, = self.ax.plot(energy, flux,
                                   color='black', linestyle='-.',
                                   label=r'Best fit, Heinze et al.')
#                                    label=r'Best fit UHECR ($\pm$ 3$\sigma$), Heinze et al.')
#             self.ax.fill_between(energy, band_min, band_max,
#                                  color='0.8')
            self.neutrino_models.append(heinze)

        elif name=='van_vliet':
            energy, flux, band_min, band_max = self.get_data('van_vliet_10')
            prot10, = self.ax.plot(energy, flux,
                                   color='orchid', linestyle='-.',
                                   label=r'10% protons, van Vliet et al.')
#             prot = self.ax.fill_between(energy, band_min, band_max,
#                                         color='orchid', alpha=0.25,
#                                         label=r'not excluded from UHECRs')
#             self.neutrino_models.append(prot)
            self.neutrino_models.append(prot10)

        elif name=='ahlers':
            energy, flux, _, _ = self.get_data('ahlers_100')
            prot100, = self.ax.plot(energy, flux,
                                    color='mediumblue', linestyle='-.',
                                    label=r'100% protons, Ahlers & Halzen')
            energy, flux, _, _ = self.get_data('ahlers_10')
            prot10, = self.ax.plot(energy, flux,
                                   color='royalblue', linestyle='-.',
                                   label=r'10% protons, Ahlers & Halzen')
            energy, flux, _, _ = self.get_data('ahlers_1')
            prot1, = self.ax.plot(energy, flux,
                                  color='cornflowerblue', linestyle='-.',
                                  label=r'1% protons, Ahlers & Halzen') # (1208.4181)
            self.neutrino_models.extend([prot100, prot10, prot1])

        elif name=='kotera':
            energy, _, band_min, band_max = self.get_data('kotera_olinto')
            compositions = self.ax.fill_between(energy, band_min, band_max,
                                                color='gray', alpha=0.25,
                                                label=r'UHECR, Olinto et al.')

            energy, flux, _, _ = self.get_data('kotera_high_e')
            kotera_high, = self.ax.plot(energy, flux,
                                        color='deeppink', linestyle='--',
                                        label=r'SFR $E_{max}=10^{21.5}$, Kotera et al.') # (1009.1382)

            energy, flux, _, _ = self.get_data('kotera_mid')
            kotera, = self.ax.plot(energy, flux,
                                   color='darkmagenta', linestyle='--',
                                   label=r'SFR $E_{max}=10^{20.5}$, Kotera et al.') # (1009.1382)
            self.neutrino_models.extend([compositions, kotera_high, kotera])

        elif name=='fang_merger':
            energy, flux, _, _ = self.get_data('fang_ns_merger')
            ns_merger, = self.ax.plot(energy, flux,
                                      color='palevioletred', linestyle=(0, (3, 5, 1, 5)),
                                      label='NS-NS merger, Fang & Metzger') # (1707.04263)
            self.neutrino_models.append(ns_merger)

        elif name=='fang_pulsar':
            energy, _, band_min, band_max = self.get_data('fang_pulsar')
            p_pulsar = self.ax.fill_between(energy, band_min, band_max,
                                            color='pink', alpha=0.5,
                                            label="Pulsar, Fang et al.") # (1311.2044)
            self.neutrino_models.append(p_pulsar)
            
        elif name=='fang_cluster':
            energy, flux, _, _ = self.get_data('fang_cluster')
            p_cluster, = self.ax.plot(energy, flux,
                                      color="mediumvioletred", zorder=1, linestyle=(0, (5, 10)),
                                      label="Clusters, Fang & Murase") # (1704.00015)
            self.neutrino_models.append(p_cluster)

        elif name=='biehl':
            energy, flux, band_min, band_max = self.get_data('biehl_tde')
            self.ax.fill_between(energy, band_min, band_max,
                                 color='thistle', alpha=0.5)
            p_tde, = self.ax.plot(energy, flux,
                                  color='darkmagenta', linestyle=':',zorder=1,
                                  label="TDE, Biehl et al.") # (1711.03555)
            self.neutrino_models.append(p_tde)

        elif name=='boncioli':
            energy, flux, band_min, band_max = self.get_data('boncioli_llgrb')
            self.ax.fill_between(energy, band_min, band_max,
                                 color='0.8')
            p_ll_grb, = self.ax.plot(energy, flux,
                                     linestyle='-.', c='k', zorder=1,
                                     label="LLGRB, Boncioli et al.") # (1808.07481)
            self.neutrino_models.append(p_ll_grb)

        elif name=='murase_agn':
            energy, flux, _, _ = self.get_data('murase_agn')
            agn, = self.ax.plot(energy, flux,
                                color="red", linestyle='--',
                                label="AGN, Murase") # (1511.01590)
            self.neutrino_models.append(agn)

        elif name=='murase_grb':
            energy, flux, _, _ = self.get_data('murase_grb_late_prompt')
            late, = self.ax.plot(energy, flux,
                                 color="saddlebrown", linestyle='-.',
                                 label="GRB afterglow-late prompt, Murase") # (0707.1140)
            energy, flux, _, _ = self.get_data('murase_grb_wind')
            wind, = self.ax.plot(energy, flux,
                                 color="goldenrod", linestyle='-.',
                                 label="GRB afterglow-wind, Murase") # (0707.1140)
            energy, flux, _, _ = self.get_data('murase_grb_ism')
            ism, = self.ax.plot(energy, flux,
                                color="gold", linestyle='-.',
                                label="GRB afterglow-ISM, Murase") # (0707.1140)
            self.neutrino_models.extend([late, wind, ism])

        else:
            raise ValueError(f"Unrecognized data set '{str(name)}'")


    def add_experiment(self, name):
        if 'ice_cube' in name:
            if self.e_power==2:
                self.ax.annotate('IceCube',
                                 xy=(1e6, 5e-8), xycoords='data',
                                 horizontalalignment='center', color='dodgerblue', rotation=0)
            if self.e_power==1:
                self.ax.annotate('IceCube',
                                 xy=(1.1e7, 2.5e-15), xycoords='data',
                                 horizontalalignment='center', color='dodgerblue', rotation=0)

        if name=='grand_10k':
            energy, flux, _, _ = self.get_data('experiments/grand_10k.txt')
            self.ax.plot(energy, flux,
                         color='saddlebrown', linestyle="--")
            if self.e_power==2:
                self.ax.annotate('GRAND 10k',
                                 xy=(5e9, 2e-8*self.e_bins), xycoords='data',
                                 horizontalalignment='left', color='saddlebrown', rotation=40)
            if self.e_power==1:
                self.ax.annotate('GRAND 10k',
                                 xy=(2e9, 5e-18*self.e_bins), xycoords='data',
                                 horizontalalignment='left', color='saddlebrown', rotation=-10)

        elif name=='grand_200k':
            energy, flux, _, _ = self.get_data('experiments/grand_200k.txt')
            self.ax.plot(energy, flux,
                         color='saddlebrown', linestyle="--")
            if self.e_power==2:
                self.ax.annotate('GRAND 200k',
                                 xy=(1e10, 3e-9*self.e_bins), xycoords='data',
                                 horizontalalignment='left', color='saddlebrown', rotation=40)

        elif name=='radar':
            energy, _, band_min, band_max = self.get_data('experiments/radar.txt')
            self.ax.fill_between(energy, band_min, band_max,
                                 facecolor='None', edgecolor='0.8', hatch='x')
            if self.e_power==2:
                self.ax.annotate('Radar',
                                 xy=(1e9, 3e-8*self.e_bins), xycoords='data',
                                 horizontalalignment='left', color='0.7', rotation=45)

        elif name=='ice_cube_ehe':
            energy, flux, _, _ = self.get_data('experiments/ice_cube_ehe.txt')
            self.ax.plot(energy, flux,
                         color='dodgerblue')

        elif name=='ice_cube_hese_data':
            energy, flux, err_min, err_max = self.get_data('experiments/ice_cube_hese.txt')
            uplimit = err_max-flux
            uplimit[np.where(err_max-flux == 0)] = 1
            uplimit[np.where(err_max-flux != 0)] = 0

            self.ax.errorbar(energy, flux*3,
                             yerr=np.asarray([flux-err_min, err_max-flux])*3, uplims=uplimit,
                             color='dodgerblue', marker='o', ecolor='dodgerblue', linestyle='None')

        elif name=='ice_cube_hese_fit':
            ice_cube_hese_range = self._get_i3_hese_range()
            energy = ice_cube_hese_range[0] / self.e_unit
            band_min = ice_cube_hese_range[1] * ice_cube_hese_range[0]**self.e_power / self.f_unit
            band_max = ice_cube_hese_range[2] * ice_cube_hese_range[0]**self.e_power / self.f_unit
            self.ax.fill_between(energy, band_min, band_max,
                                 hatch='\\', edgecolor='dodgerblue', facecolor='azure')
            flux = self._i3_nu_fit(ice_cube_hese_range[0], offset=2.46, slope=-2.92) * ice_cube_hese_range[0]**self.e_power / self.f_unit
            self.ax.plot(energy, flux,
                         color='dodgerblue')

        elif name=='ice_cube_mu_fit':
            ice_cube_mu_range = self._get_i3_mu_range()
            energy = ice_cube_mu_range[0] / self.e_unit
            band_min = ice_cube_mu_range[1] * ice_cube_mu_range[0]**self.e_power / self.f_unit
            band_max = ice_cube_mu_range[2] * ice_cube_mu_range[0]**self.e_power / self.f_unit
            self.ax.fill_between(energy, band_min, band_max,
                                 edgecolor='dodgerblue', facecolor='azure')
            flux = self._i3_nu_fit(ice_cube_mu_range[0], offset=1.01, slope=-2.19) * ice_cube_mu_range[0]**self.e_power / self.f_unit
            self.ax.plot(energy, flux,
                         color='dodgerblue')

        elif name=='anita':
            energy, flux, _, _ = self.get_data('experiments/anita.txt')
            self.ax.plot(energy, flux,
                         color='darkorange')
            if self.e_power==2:
                self.ax.annotate('ANITA I - III',
                                 xy=(4e9, 1.2e-6*self.e_bins), xycoords='data',
                                 horizontalalignment='left', color='darkorange')
            if self.e_power==1:
                self.ax.annotate('ANITA I - III',
                                 xy=(3e9, 1e-15*self.e_bins), xycoords='data',
                                 horizontalalignment='left', color='darkorange')

        elif name=='auger':
            energy, flux, _, _ = self.get_data('experiments/auger.txt')
            self.ax.plot(energy, flux,
                         color='forestgreen')
            if self.e_power==2:
                self.ax.annotate('Auger',
                                 xy=(1.2e8, 5e-8*self.e_bins), xycoords='data',
                                 horizontalalignment='left', color='forestgreen', rotation=-40)
            if self.e_power==1:
                self.ax.annotate('Auger',
                                 xy=(3e10, 8e-18*self.e_bins), xycoords='data',
                                 horizontalalignment='left', color='forestgreen', rotation=-8)

        elif name=='ara':
            energy, flux, _, _ = self.get_data('experiments/ara_ming.txt')
            self.ax.plot(energy, flux,
                         color='#2288AA')
            energy, flux, _, _ = self.get_data('experiments/ara_projected.txt')
            self.ax.plot(energy, flux,
                         color='#2288AA', linestyle='--')
            if self.e_power==2:
                self.ax.annotate('ARA (2x1yr)',
                                 xy=(3e8, 4e-7*self.e_bins), xycoords='data',
                                 horizontalalignment='left', color='#2288AA', rotation=-10)
                self.ax.annotate('ARA (2x4yr)',
                                 xy=(3e8, 1e-7*self.e_bins), xycoords='data',
                                 horizontalalignment='left', color='#2288AA', rotation=-10)
            if self.e_power==1:
                self.ax.annotate('ARA (2x1yr)',
                                 xy=(2.5e8, 1.5e-15*self.e_bins), xycoords='data',
                                 horizontalalignment='left', color='#2288AA', rotation=-50)
                self.ax.annotate('ARA (2x4yr)',
                                 xy=(2e8, 2e-16*self.e_bins), xycoords='data',
                                 horizontalalignment='left', color='#2288AA', rotation=-50)

        else:
            raise ValueError(f"Unrecognized data set '{str(name)}'")


    def build_base_plot(self, group='clean', experiments=None, models=None):
        if group=='all':
            if experiments is None:
                experiments = ['grand_10k', 'grand_200k', 'radar', 'anita', 'auger', 'ara',
                               'ice_cube_ehe', 'ice_cube_hese_data', 'ice_cube_hese_fit', 'ice_cube_mu_fit']
            if models is None:
                models = ['heinze', 'ahlers', 'kotera', 'van_vliet']
        elif group=='clean':
            if experiments is None:
                experiments = ['grand_10k', 'anita', 'auger', 'ara',
                               'ice_cube_ehe', 'ice_cube_hese_data', 'ice_cube_mu_fit']
            if models is None:
                models = ['kotera', 'ahlers']
        elif group=='ice_cube':
            if experiments is None:
                experiments = ['ice_cube_ehe', 'ice_cube_hese_data', 'ice_cube_hese_fit', 'ice_cube_mu_fit']
            if models is None:
                models = []
        elif group=='ming':
            if experiments is None:
                experiments = ['ice_cube_ehe', 'ice_cube_hese_data', 'ice_cube_mu_fit', 'ara', 'anita', 'auger']
            if models is None:
                models = ['ahlers', 'kotera']
        elif group=='rno':
            if experiments is None:
                experiments = ['grand_10k', 'anita', 'auger', 'ara',
                               'ice_cube_ehe', 'ice_cube_hese_data', 'ice_cube_mu_fit']
            if models is None:
                models = ['heinze', 'ahlers']
        elif group=='ara':
            if experiments is None:
                experiments = ['anita', 'auger', 'ara',
                               'ice_cube_ehe', 'ice_cube_hese_data', 'ice_cube_mu_fit']
            if models is None:
                models = ['kotera', 'ahlers']
        elif group=='ara_src':
            if experiments is None:
                experiments = ['anita', 'auger', 'ara',
                               'ice_cube_ehe', 'ice_cube_hese_data', 'ice_cube_mu_fit']
            if models is None:
                models = ['fang_pulsar', 'murase_agn', 'murase_grb', 'fang_cluster']
        elif group=='test':
            if experiments is None:
                experiments = ['grand_10k', 'anita', 'auger', 'ara',
                               'ice_cube_ehe', 'ice_cube_hese_data', 'ice_cube_mu_fit']
            if models is None:
                models = ['biehl', 'boncioli', 'fang_merger', 'fang_pulsar', 'fang_cluster'] 
        else:
            if experiments is None:
                experiments = []
            if models is None:
                models = []
        for name in models:
            self.add_model(name)
        for name in experiments:
            self.add_experiment(name)


    def add_limit(self, name, energies, veffs, stations=1, years=1, color=None, linestyle=None, label=None):
        print(f"Energies: {energies}")
        print(f"Veffs: {veffs}")
        aeffs = flux_tool.veff_to_aeff(energies, veffs)
        print(f"Aeffs: {aeffs}")
        print(f"Total Aeffs: {aeffs*stations*years}")
        limits = flux_tool.flux_sensitivity(energies, aeffs*stations,
                                            livetime=years*flux_tool.units.yr)
        print(f"Limits: {limits}")
        limits *= energies**self.e_power # adjust for energy scaling of plot
        print(f"Scaled Limits: {limits}")

        # Scale by energy bin size
        log_energy = np.log10(energies)
        d_log_energy = np.diff(log_energy)
        bins_per_decade = 1/d_log_energy[0]
        limits *= self.e_bins / bins_per_decade

        if label is None:
            label = f"{name}: {stations} stations, {years} years"

        # Plot limit
        _plt, = self.ax.plot(energies / self.e_unit,
                             limits / self.f_unit,
                             color=color, linestyle=linestyle,
                             label=label,
                             linewidth=3,
                             zorder=100+len(self.custom_limits))
        self.custom_limits.append(_plt)


    def title(self, title, size=None):
        self.ax.set_title(title)
        if size is None:
            size = self.font_size
        self.ax.title.set_fontsize(size)

    def show(self, legend_size=12, save_name=None, *args, **kwargs):
        if self.e_power==2:
            self.ax.add_artist(plt.legend(handles=self.neutrino_models, loc=4, fontsize=legend_size))
            self.ax.add_artist(plt.legend(handles=self.custom_limits, loc=2, fontsize=legend_size))
        elif self.e_power==1:
            self.ax.add_artist(plt.legend(handles=self.neutrino_models, loc=3, fontsize=legend_size))
            self.ax.add_artist(plt.legend(handles=self.custom_limits, loc=1, fontsize=legend_size))
        plt.tight_layout()
        if save_name is not None:
            plt.savefig(save_name, *args, **kwargs)
        plt.show()




if __name__=='__main__':
    my_data = {
        "energy": 10**np.array([16.5, 17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0]) * flux_tool.units.eV,
        "veff": np.array([1.105E-1, 5.195E-1, 1.826E+0, 5.259E+0,
                          1.106E+1, 2.077E+1, 3.431E+1, 5.073E+1]) * flux_tool.units.km**3*flux_tool.units.sr,
    }

    figure = LimitFigure(e_power=1, xlims=(1e6, 1e11), ylims=(1e-19, 2e-14), font_size=16, tick_size=14)
    figure.build_base_plot('ara')
    figure.add_limit('My Data', my_data['energy'], my_data['veff'],
                     stations=100, years=5, color='black', linestyle=None)
    figure.title("Trigger Level Sensitivity")
    figure.show(legend_size=10, save_name='my_sensitivity_plot_e1.pdf')

    figure = LimitFigure(e_power=2, xlims=(1e5, 1e11), ylims=(1e-11, 2e-6), font_size=16, tick_size=14)
    figure.build_base_plot('rno')
    figure.add_limit('My Data', my_data['energy'], my_data['veff'],
                     stations=100, years=5, color='black', linestyle=None)
    figure.title("Trigger Level Sensitivity")
    figure.show(legend_size=10, save_name='my_sensitivity_plot_e2.pdf')
