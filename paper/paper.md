---
title: 'The Pulsar Signal Simulator: A Python package for simulating radio signal data from pulsars'
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
  - name: Amelia Henkel
    orcid:
    affiliation: 4
  - name: Cassidy Wagner
    orcid:
    affiliation: 4
  - name: Paul Brook
    orcid:
    affiliation: 4
  - name: Jacob Hesse
    orcid:
    affiliation: 4
  - name: Maura Mclaughlin
    orcid:
    affiliation: 2
  - name: Michael T. Lam
    orcid:
    affiliation: 4

affiliations:
 - name: University of Washington Bothell
   index: 1
 - name: West Virginia University
   index: 2
 - name: Widener College
   index: 3
 - name: Widener College
   index: 4
 - name: Widener College
   index: 5
 - name: Rochester Institute of Technology
   index: 6
date: 16 July 2020
bibliography: paper.bib
---

# Summary

Pulsar observations have been a hallmark of radio astronomy since Jocelyn Bell
discovered the repeating signal from PSR B19 in 1967. Many radio telescopes have
been designed with these observations in mind. The phenomenal precision of pulsar
rotational periods has allowed for many notable observations stemming from strong
field general relativity.

Gravitational waves are quickly changing the way that we view the wider
universe, enabling observations of compact objects in highly relativistic
scenarios. Gravitational-wave detectors measure the minuscule, time-dependent
perturbations to the spacetime metric. These detectors have long been
characterized by a sensitivity curve, a plot in the frequency domain, which
summarizes their ability to *detect* a given signal. Pulsar timing arrays
(PTAs) are collections of highly precise millisecond pulsars regularly
monitored for shifts in the spin period of pulsars indicative of gravitational
waves in the nanohertz regime. See @hobbs and @burke-spolaor for a review of
pulsar timing arrays and the astrophysics of nanohertz gravitational waves.

The ``PsrSigSim`` is a Python package for simulating radio telescope data from pulsar timing observations, largely based on the formalism presented in [@shapiro-key:2020].
It uses standard Python packages, such as ``Numpy`` [@numpy] and ``Astropy``
[@astropy] to simulate radio pulses from a neutron, the propagation effects upon that signal
from the interstellar medium and the various Doppler shifts that come from the motion of the pulsar and the Earth.  

# Acknowledgements

JSH , BJS, PTB,  and  acknowledge subawards from the University of Wisconsin-Milwaukee for the NSF NANOGrav Physics Frontier Center (NSF PFC-1430284). Finally, we thank Joseph D. Romano, Tim Pennucci and for useful discussions and sharing preliminary code.

# References
