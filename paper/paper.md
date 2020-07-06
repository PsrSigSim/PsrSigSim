---
title: 'The Pulsar Signal Simulator: A Python package for simulating radio signals from pulsars'
tags:
  - Python
  - pulsars
  - astronomy
  - gravitational waves
  - pulsar timing arrays
authors:
  - name: Jeffrey S. Hazboun
    orcid: 0000-0003-2742-3321
    affiliation: 1
  - name: Brent Shapiro-Albert
    orcid:
    affiliation: 2
  - name: Paul T. Baker
    orcid:
    affiliation: 3
  - name:
    orcid:
    affiliation: 4

affiliations:
 - name: University of Washington Bothell
   index: 1
 - name: West Virginia University
   index: 2
 - name: Widener College
   index: 3
date: 16 September 2019
bibliography: paper.bib
---

# Summary

Gravitational waves are quickly changing the way that we view the wider
universe, enabling observations of compact objects in highly relativistic
scenarios. Gravitational-wave detectors measure the minuscule, time-dependent
perturbations to the spacetime metric. These detectors have long been
characterized by a sensitivity curve, a plot in the frequency domain, which
summarizes their ability to *detect* a given signal. Pulsar timing arrays
(PTAs) are collections of highly precise millisecond pulsars regularly
monitored for shifts in the spin period of pulsars indicative of gravitational
waves in the nanohertz regime. See @hobbs and @burke-spolaor for a review of
pulsar timing arrays and the astrophysics of nanohertz gravitational waves. The
sensitivity curves for PTAs are often overly simplified in the literature,
lacking detailed information about the fit to a pulsar's timing parameters and
assuming identical pulsar noise characteristics.

The ``PsrSigSim`` is a Python package for calculating and building accurate PTA
sensitivity curves, largely based on the formalism presented in [@hazboun:2019].
This software is designed for use by astronomers interested in sources of
nanohertz gravitational waves and PTA data analysts alike.
It uses standard Python packages, such as ``Numpy`` [@numpy] and ``Astropy``
[@astropy] to build sensitivity curves from generic PTAs of individually
constructed pulsars.

<!--- The source code for ``Hasasia`` has been archived to Zenodo with the linked DOI: [@zenodo] --->

# Acknowledgements

JSH , and  acknowledge subawards from the University of Wisconsin-Milwaukee for the NSF NANOGrav Physics Frontier Center (NSF PFC-1430284). Finally, we thank Joseph D. Romano,  for useful discussions and sharing preliminary code.

# References
