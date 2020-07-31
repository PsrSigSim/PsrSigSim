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
    orcid: 0000-0002-7283-1124
    affiliation: 2
  - name: Paul T. Baker
    orcid: 0000-0003-2745-753X
    affiliation: 3
  - name: Amelia M. Henkel
    orcid:
    affiliation: 4
  - name: Cassidy M. Wagner
    orcid: 0000-0002-1186-2082
    affiliation: 5
  - name: Paul Brook
    orcid: 0000-0003-3053-6538
    affiliation: 2
  - name: Jacob Hesse
    orcid:
    affiliation: 1
  - name: Maura Mclaughlin
    orcid: 0000-0001-7697-7422
    affiliation: 2
  - name: Michael T. Lam
    orcid: 0000-0003-0721-651X
    affiliation: 6

affiliations:
 - name: Physical Sciences Division, University of Washington Bothell, 18115 Campus Way NE, Bothell, WA 98011, USA
   index: 1
 - name: Department of Physics and Astronomy, West Virginia University, P.O. Box 6315, Morgantown, WV 26506, USA
   index: 2
 - name: Department of Physics and Astronomy, Widener University, One University Place, Chester, PA 19013, USA
   index: 3
 - name: Widener College
   index: 4
 - name: Department of Astronomy, University of Illinois, 1002 West Green Street, Urbana IL 61802, USA
   index: 5
 - name: School of Physics and Astronomy, Rochester Institute of Technology, Rochester, NY 14623, USA
   index: 6
 - name: Laboratory for Multiwavelength Astronomy, Rochester Institute of Technology, Rochester, NY 14623, USA
   index: 7
date: 16 July 2020
bibliography: paper.bib
---

# Summary

Pulsar observations have been a hallmark of radio astronomy since Jocelyn Bell-Burnell
discovered the repeating signal from PSR B1919+21 in 1967. Many radio telescopes have
been designed with observations of these sources in mind. The phenomenal precision of pulsar
rotational periods has allowed for many notable observations stemming from strong
field general relativity, including a Nobel prize for the discovery of the first binary neutron star, which led to the first evidence for gravitational waves [@taylor:1982].

Pulsars continue to be used as probes of strong field gravity as the constituents of pulsar timing arrays (PTAs). PTAs are collections of highly precise millisecond pulsars regularly
monitored for shifts in their pulse arrival times indicative of gravitational
waves in the nanohertz regime. See @hobbs and @burke-spolaor for a review of
pulsar timing arrays and the astrophysics of nanohertz gravitational waves.

In order to use these neutron stars as a galactic-scale gravitational wave observatory it is important to understand as many sources of noise in the timing of pulses. This timing depends on many parameters in a timing process and deterministic model that includes aspects of the pulse profile, intrinsic flux, orbital parameters of the pulsar, proper motion across the sky, effects of the interstellar medium (ISM), position of the observatory with respect to the solar system barycenter, as well as various aspects of the telescope and backend processing. In order to understand the covariances and noise terms in this model an end-to-end simulation is needed.

The ``PsrSigSim`` is a Python package for simulating radio telescope data from pulsar timing observations, largely based on the formalism presented in [@shapiro-key:2020].
It uses standard Python packages, such as ``Numpy`` [@numpy] and ``Astropy``
[@astropy] to simulate radio pulses from a neutron star, the propagation effects upon that signal
from the interstellar medium and the various Doppler shifts that come from the motion of the pulsar and the Earth. ``Signals``-objects are passed to ``Pulsar``-objects and `Telescope`-objects to build a realistic pulsar signal from scratch. `ISM`-objects use various signal processing techniques, such as the Fourier shift theorem and convolution of pulse profiles with decaying exponentials in order to simulate the frequency dependent effects of the interstellar medium on pulse arrival times.

The modularity of the simulation allows for enumerable options for investigating the processing of pulsar signals. The signals from different pulsars can be treated by the same telescope in order to investigate how various pulsars would look using a particular telescope configuration. A user can either use a provided telescope, or build one with their own specifications. This modularity also makes the `PsrSigSim` an excellent tool for educating students about the physics of pulsar signals and how they are observed and processed.

The documentation for the `PsrSigSim` includes a number of beginner level tutorials for getting started using the package.

# Acknowledgements

JSH , BJS, PTB,  and  acknowledge subawards from the NSF NANOGrav Physics Frontier Center (NSF PFC-1430284). Finally, we thank Joseph D. Romano, Tim Pennucci and for useful discussions and sharing preliminary code.

# References
