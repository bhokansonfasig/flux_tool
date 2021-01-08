"""Class for neutrino interaction physics calculations"""

from enum import Enum

import numpy as np
import scipy.constants

from . import units


def get_from_enum(value, enum):
    """
    Find the enum value given some representation of it.

    Transforms the given `value` into the corresponding value from the `enum`
    by checking the type of `value` given.

    Parameters
    ----------
    value
        Representation of the desired `enum` value. If already a member of
        `enum`, no change. If ``str``, assumed to be a name in the `enum`.
        Otherwise, assumed to be a value type of the `enum`.
    enum : Enum
        Python ``Enum`` to compare names values with.

    Returns
    -------
    Enum value
        Value in the `enum` represented by the given `value`.

    Examples
    --------
    >>> from enum import Enum
    >>> class Color(Enum):
    ...     red = 1
    ...     green = 2
    ...     blue = 3
    >>> get_from_enum(Color.red, Color)
    <Color.red: 1>
    >>> get_from_enum("green", Color)
    <Color.green: 2>
    >>> get_from_enum(3, Color)
    <Color.blue: 3>

    """
    if isinstance(value, enum):
        return value
    elif isinstance(value, str):
        return enum[value]
    else:
        return enum(value)


class NeutrinoInteraction:
    """
    Class for storing and calculating neutrino interaction parameters.

    Parameters
    ----------
    neutrino_type
        Identification value of the neutrino type. Values should be from the
        ``NeutrinoInteraction.NeutrinoType`` enum, but integer or string values
        may work if carefully chosen.
    interaction_type
        Identification value of the neutrino's interaction type. Values should
        be from the ``NeutrinoInteraction.InteractionType`` enum, but integer
        or string values may work if carefully chosen.
    energy : float
        Energy (GeV) of the neutrino.

    Attributes
    ----------
    neutrino : NeutrinoInteraction.NeutrinoType
        Identification value of the neutrino type.
    interaction : NeutrinoInteraction.InteractionType
        Identification value of the neutrino's interaction type.
    energy : float
        Energy (GeV) of the neutrino.

    """
    class NeutrinoType(Enum):
        """
        Enum containing possible neutrino types.

        Values based on the PDG particle numbering scheme.
        http://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf

        Attributes
        ----------
        nu_e, electron_neutrino
        nu_e_bar, electron_antineutrino
        nu_mu, muon_neutrino
        nu_mu_bar, muon_antineutrino
        nu_tau, tau_neutrino
        nu_tau_bar, tau_antineutrino
        unknown, undefined

        """
        undefined = 0
        unknown = 0
        electron_neutrino = 12
        nu_e = 12
        electron_antineutrino = -12
        nu_e_bar = -12
        muon_neutrino = 14
        nu_mu = 14
        muon_antineutrino = -14
        nu_mu_bar = -14
        tau_neutrino = 16
        nu_tau = 16
        tau_antineutrino = -16
        nu_tau_bar = -16

    class InteractionType(Enum):
        """
        Enum containing possible interaction types.

        Attributes
        ----------
        cc, charged_current
        nc, neutral_current
        unknown, undefined

        """
        undefined = 0
        unknown = 0
        charged_current = 1
        cc = 1
        neutral_current = 2
        nc = 2

    def __init__(self, neutrino_type, interaction_type, energy):
        self.neutrino = neutrino_type
        self.interaction = interaction_type
        self.energy = energy

    @property
    def neutrino(self):
        """
        Identification value of the neutrino type.

        Should always be a value from the ``NeutrinoType`` enum. Setting with
        integer or string values may work if carefully chosen.

        """
        return self._neutrino_id

    @neutrino.setter
    def neutrino(self, neutrino_id):
        if neutrino_id is None:
            self._neutrino_id = self.NeutrinoType.undefined
        else:
            self._neutrino_id = get_from_enum(neutrino_id, self.NeutrinoType)

    @property
    def interaction(self):
        """
        Identification value of the neutrino's interaction type.

        Should always be a value from the ``InteractionType`` enum. Setting
        with integer or string values may work if carefully chosen.

        """
        return self._interaction_id

    @interaction.setter
    def interaction(self, interaction_id):
        if interaction_id is None:
            self._interaction_id = self.InteractionType.undefined
        else:
            self._interaction_id = get_from_enum(interaction_id, self.InteractionType)


    @property
    def total_cross_section(self):
        """
        The total cross section of the neutrino.

        Calculation is determined by whether the neutrino is a neutrino
        or antineutrino and is dependent on the energy of the neutrino.
        Combines the charged-current and neutral-current cross sections.
        Based on Equation 7 and Table III of the CTW 2011 paper.

        """
        # Total cross section should be sum of nc and cc cross sections

        # Neutrino
        if self.neutrino.value>0:
            c_0_cc = -1.826
            c_0_nc = -1.826
            c_1_cc = -17.31
            c_1_nc = -17.31
            c_2_cc = -6.406
            c_2_nc = -6.448
            c_3_cc = 1.431
            c_3_nc = 1.431
            c_4_cc = -17.91
            c_4_nc = -18.61
        # Antineutrino
        elif self.neutrino.value<0:
            c_0_cc = -1.033
            c_0_nc = -1.033
            c_1_cc = -15.95
            c_1_nc = -15.95
            c_2_cc = -7.247
            c_2_nc = -7.296
            c_3_cc = 1.569
            c_3_nc = 1.569
            c_4_cc = -17.72
            c_4_nc = -18.30
        else:
            raise ValueError("Unable to calculate cross section without a"+
                             " particle type")
        # Calculate cross section based on CTW 2011
        eps = np.log10(self.energy / units.GeV)
        log_term_cc = np.log(eps - c_0_cc)
        power_cc = (c_1_cc + c_2_cc*log_term_cc + c_3_cc*log_term_cc**2
                    + c_4_cc/log_term_cc)
        log_term_nc = np.log(eps - c_0_nc)
        power_nc = (c_1_nc + c_2_nc*log_term_nc + c_3_nc*log_term_nc**2
                    + c_4_nc/log_term_nc)
        return (10**power_cc + 10**power_nc) * units.cm**2

    @property
    def cross_section(self):
        """
        The cross section of the neutrino.

        Calculation is determined by whether the neutrino is a neutrino
        or antineutrino and what type of interaction it produces, and is
        dependent on the energy of the neutrino. Based on Equation 7 and
        Table III of the CTW 2011 paper.

        """
        # Neutrino
        if self.neutrino.value>0:
            if self.interaction==self.InteractionType.charged_current:
                c_0 = -1.826
                c_1 = -17.31
                c_2 = -6.406
                c_3 = 1.431
                c_4 = -17.91
            elif self.interaction==self.InteractionType.neutral_current:
                c_0 = -1.826
                c_1 = -17.31
                c_2 = -6.448
                c_3 = 1.431
                c_4 = -18.61
            else:
                raise ValueError("Unable to calculate cross section without an"
                                 +" interaction type")
        # Antineutrino
        elif self.neutrino.value<0:
            if self.interaction==self.InteractionType.charged_current:
                c_0 = -1.033
                c_1 = -15.95
                c_2 = -7.247
                c_3 = 1.569
                c_4 = -17.72
            elif self.interaction==self.InteractionType.neutral_current:
                c_0 = -1.033
                c_1 = -15.95
                c_2 = -7.296
                c_3 = 1.569
                c_4 = -18.30
            else:
                raise ValueError("Unable to calculate cross section without an"
                                 +" interaction type")
        else:
            raise ValueError("Unable to calculate cross section without a"+
                             " neutrino type")
        # Calculate cross section based on CTW 2011
        eps = np.log10(self.energy / units.GeV)
        log_term = np.log(eps - c_0)
        power = c_1 + c_2*log_term + c_3*log_term**2 + c_4/log_term
        return (10**power) * units.cm**2


    @property
    def total_interaction_length(self):
        """
        The interaction length of the neutrino.

        The interaction length is calculated in water equivalent material.
        Calculation is determined by whether the neutrino is a neutrino or
        antineutrino and is dependent on the energy of the neutrino. Combines
        the charged-current and neutral-current interaction lengths.

        """
        # Water equivalent density is 1 g/cm^3
        # The approximate number of nucleons per gram is Avogadro's number,
        # so the nucleon density is approximately 1 nucleon/cm^3 in water.
        # Ultimately then the interaction length can be caluclated as
        # 1 / NA / cross_section
        return 1 / (scipy.constants.N_A/units.cm**3) / self.total_cross_section

    @property
    def interaction_length(self):
        """
        The interaction length of the neutrino interaction.

        The interaction length is calculated in water equivalent material.
        Calculation is determined by whether the neutrino is a neutrino or
        antineutrino and what type of interaction it produces, and is dependent
        on the energy of the neutrino.

        """
        return 1 / (scipy.constants.N_A/units.cm**3) / self.cross_section
