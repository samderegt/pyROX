.. pyROX: Rapid Opacity X-sections documentation master file, created by
   sphinx-quickstart on Mon Apr  7 13:53:59 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pyROX: Rapid Opacity X-sections
===============================

Welcome to the pyROX documentation. pyROX is a Python package for calculating
opacity cross-sections for sub-stellar atmospheres. It can handle line-by-line
opacities using various databases (ExoMol, HITRAN, Kurucz) and also supports
the calculation of collision-induced absorption.


Installation
------------

To install pyROX, clone the `repository from GitHub <https://github.com/samderegt/pyROX>`_:

.. code-block:: bash

   git clone https://github.com/samderegt/pyROX

Next, navigate to the cloned directory and install the package using pip:

.. code-block:: bash

   cd pyROX
   pip install -e .

The ``-e`` will install pyROX in editable mode, allowing you to make changes to the
source code and see them reflected in your Python environment immediately.


License and Acknowledgements
----------------------------
pyROX is available under the `MIT License <https://github.com/samderegt/pyROX/LICENSE>`_. 

The package is developed and maintained by Sam de Regt, with contributions from
Siddharth Gandhi, Louis Siebenaler, and Darío González Picos.

If you used pyROX for your research, please cite the following paper:

.. code-block:: bibtex

   PLACEHOLDER
   

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