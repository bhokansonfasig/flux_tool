# flux_tool

A tool for calculating sensitivities and neutrino counts, given effective volumes or areas of a detector.


## Python Package

The `flux_tool` directory should act as a python package, so you can `import flux_tool` to get access to the necessary classes and functions. Each of these also has a docstring if you peek into the code, which should help a bit with explaining what parameters to pass. Below I describe the basic functionality:

`veff_to_aeff` is a function for converting effective volumes into effective areas, which will be necessary as the other functions take effective areas as input.

`flux_sensitivity` is a function for calculating the sensitivity of a detector or station over some livetime (i.e. the upper bound that would be set by the detector).

`neutrino_count` is a function for calculating the expected number of neutrinos which will be observed by a detector or station over some livetime, given a source model of the neutrinos.

`Model` is a class for creating model objects, which store the necessary information regarding the expected flux of neutrinos from a given source.


## Examples

There are two example scripts included in this repository which show a detailed example of how these functions can be used:

The first is `limit_figure.py` which can be used to show the limit which will be set by a given detector (or a handful of detectors) against other detectors in the field and some neutrino models. Plots like this are frequently used in proposals and papers, so being able to produce such plots is handy. The class in this script can handle two types of energy weighting of the fluxes as well: E^2 and E^1. Different people like one or the other of these, so being able to produce both is also nice.

The second example is `yearly_neutrinos.py`, which is a bit more specialized in its current state. This is used to plot the expected number of neutrinos observed by a detector (in this case ARA1-5) over time, for various models. This plot is mostly only interesting for detectors which change over time (e.g. deploying more stations), but for a mostly stable detector it is probably sufficient to provide a simple "neutrinos per year" number, which can be done with the `neutrino_count` function much more simply.


## Units

There is some units handling baked into the package, which can be useful or cumbersome depending on your perspective. If you want to ignore it altogether, the base units are 1 meter, 1 second, 1 GeV, and 1 gram. If you have your inputs in these units, you can get away with it and your outputs will be in terms of these units as well.

If you do use the units system, then when you are inputting a number (say an energy of 10 PeV), you need to multiply your number (10) by the unit (units.PeV). Then when you get an output, you need to divide that number by the units you want to express it in.
