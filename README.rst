=========
PsrSigSim
=========


.. image:: https://img.shields.io/pypi/v/psrsigsim.svg
        :target: https://pypi.python.org/pypi/psrsigsim

.. image:: https://img.shields.io/travis/PsrSigSim/PsrSigSim.svg
        :target: https://travis-ci.org/PsrSigSim/PsrSigSim

.. image:: https://readthedocs.org/projects/psrsigsim/badge/?version=latest
        :target: https://psrsigsim.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

The NANOGrav Pulsar Signal Simulator
------------------------------------

* Free software: MIT license
* Documentation: https://psrsigsim.readthedocs.io.

The Pulsar Signal Simulator (`PsrSigSim`) is a Python package developed by the
North American Nanohertz Observatory for Gravitational Waves (NANOGrav). This
software is free to use and is designed as a tool for simulating realistic
pulsar signals, as well as educating students and researchers about those
signals. Various models from the copious pulsar literature are used for
simulating emission, propagation effects and data processing artifacts.

Goals
-----

* **Investigate sources of time-of-arrival errors:** Simulate various sources of time of arrival errors, including interstellar medium effects, intrinsic pulsar noise, various pulsar emission mechanisms and gravitational waves. Simulate instrumental noise sources, including radio frequency interference, radiometer noise and backend sampling effects.
* **Education and Outreach:** Allow users to build pulsar signals piece-by-piece to demonstrate to students how pulsar signals change as they propagate and how they are changed by the signal processing done at the telescope. Make materials for talks and outreach including signal diagrams and sound files of sonified pulsars.
* **Search algorithms and search training:** Test new search algorithms on signals with known parameters. Use simulated signals for search training. Assess the sensitivity of current search algorithms with simulated signals.

Code of Conduct
---------------
The `PsrSigSim` community expects contributors to follow a `Code of Conduct`_ outlined with our documentation.

Credits
-------
Development Team: Jeffrey S. Hazboun, Brent Shapiro-Albert, Paul T. Baker

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _`Code of Conduct`: https://psrsigsim.readthedocs.io/en/latest/code_of_conduct.html
.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
