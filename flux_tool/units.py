"""Standardized units for calculation"""


############################### UNIT DEFINITIONS ###############################

base_values = {
    "meter": 1,
    "second": 1,
    "electronvolt": 1e-9,
    "gram": 1,
}

base_shorthands = {
    "m": "meter",
    "s": "second",
    "eV": "electronvolt",
    "g": "gram",
}

modifier_factors = {
    "yocto": 1e-24,
    "zepto": 1e-21,
    "atto": 1e-18,
    "femto": 1e-15,
    "pico": 1e-12,
    "nano": 1e-9,
    "micro": 1e-6,
    "milli": 1e-3,
    "centi": 1e-2,
    "deci": 1e-1,
    "deca": 1e1,
    "hecto": 1e2,
    "kilo": 1e3,
    "mega": 1e6,
    "giga": 1e9,
    "tera": 1e12,
    "peta": 1e15,
    "exa": 1e18,
    "zetta": 1e21,
    "yotta": 1e24,
}

modifier_shorthands = {
    "y": "yocto",
    "z": "zepto",
    "a": "atto",
    "f": "femto",
    "p": "pico",
    "n": "nano",
    "u": "micro",
    "m": "milli",
    "c": "centi",
    "d": "deci",
    "da": "deca",
    "h": "hecto",
    "k": "kilo",
    "M": "mega",
    "G": "giga",
    "T": "tera",
    "P": "peta",
    "E": "exa",
    "Z": "zetta",
    "Y": "yotta",
}

special_values = {
    "minute": 60 * base_values['second'],
    "hour": 60 * 60 * base_values['second'],
    "day": 24 * 60 * 60 * base_values['second'],
    "year": 365 * 24 * 60 * 60 * base_values['second'],
    "steradian": 1,
}

special_shorthands = {
    "min": "minute",
    "hr": "hour",
    "yr": "year",
    "sr": "steradian"
}

############################# END UNIT DEFINITIONS #############################


shorthand_base_values = {}
for short, long in base_shorthands.items():
    shorthand_base_values[short] = base_values[long]

shorthand_modifier_factors = {}
for short, long in modifier_shorthands.items():
    shorthand_modifier_factors[short] = modifier_factors[long]

shorthand_special_values = {}
for short, long in special_shorthands.items():
    shorthand_special_values[short] = special_values[long]


def get_modified_value(name, bases, prefixes, allow_plurals=False):
    # Check for matching base name
    for base, value in bases.items():
        if name.endswith(base) or (allow_plurals and name.endswith(base+"s")):
            if allow_plurals and name[-1]=='s':
                name = name[:-1]
            prefix = name[:-len(base)]
            if prefix=="":
                return value
            # Check for prefix modifiers
            for match, mod in prefixes.items():
                if prefix==match:
                    return mod*value

def get_unmodified_value(name, bases, allow_plurals=False):
    # Check for matching name
    for base, value in bases.items():
        if name==base or (allow_plurals and name==base+"s"):
            return value


def __getattr__(name):
    value = None

    # Check base values for matching base name
    value = get_modified_value(name,
                               bases=base_values,
                               prefixes=modifier_factors,
                               allow_plurals=True)
    if value is not None:
        return value

    # Check shorthand base values for matching base name
    value = get_modified_value(name,
                               bases=shorthand_base_values,
                               prefixes=shorthand_modifier_factors,
                               allow_plurals=False)
    if value is not None:
        return value

    # Check special values for matching name
    value = get_unmodified_value(name,
                                 bases=special_values,
                                 allow_plurals=True)
    if value is not None:
        return value

    # Check special value shorthands for matching name
    value = get_unmodified_value(name,
                                 bases=shorthand_special_values,
                                 allow_plurals=False)
    if value is not None:
        return value

    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")
