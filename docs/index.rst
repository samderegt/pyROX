.. pyROX: Rapid Opacity X-sections documentation master file, created by
   sphinx-quickstart on Mon Apr  7 13:53:59 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pyROX: Rapid Opacity X-sections
===============================

Welcome to the pyROX documentation! pyROX is a Python package for calculating opacity cross-sections for sub-stellar atmospheres. It can handle line opacities using various databases (ExoMol, HITRAN, Kurucz) and also supports the calculation of collision-induced absorption coefficients (HITRAN, Borysow).


Installation
------------

To install pyROX, clone the `repository from GitHub <https://github.com/samderegt/pyROX>`_:

.. code-block:: bash

   git clone https://github.com/samderegt/pyROX

Next, navigate to the cloned directory and install the package using pip:

.. code-block:: bash

   cd pyROX
   pip install .


License and Acknowledgements
----------------------------
pyROX is available under the `MIT License <https://github.com/samderegt/pyROX/LICENSE>`_. 

The package is developed and maintained by Sam de Regt, with contributions from
Siddharth Gandhi, Louis Siebenaler, and Darío González Picos.

If you used pyROX for your research, please cite the following paper:

.. code-block:: bibtex

   @ARTICLE{2025arXiv251020870D,
         author = {{de Regt}, Sam and {Gandhi}, Siddharth and {Siebenaler}, Louis and {Gonz{\'a}lez Picos}, Dar{\'\i}o},
         title = "{pyROX: Rapid Opacity X-sections}",
         journal = {arXiv e-prints},
      keywords = {Instrumentation and Methods for Astrophysics, Earth and Planetary Astrophysics},
            year = 2025,
         month = oct,
            eid = {arXiv:2510.20870},
         pages = {arXiv:2510.20870},
            doi = {10.48550/arXiv.2510.20870},
   archivePrefix = {arXiv},
         eprint = {2510.20870},
   primaryClass = {astro-ph.IM},
         adsurl = {https://ui.adsabs.harvard.edu/abs/2025arXiv251020870D},
         adsnote = {Provided by the SAO/NASA Astrophysics Data System}
   }
   

.. toctree::
   :maxdepth: 2
   :caption: Guide

   tutorial
   databases
   parallelisation
   contributing

.. toctree::
   :maxdepth: 2
   :caption: Code Documentation

   GitHub <https://github.com/samderegt/pyROX>