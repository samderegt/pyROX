===============
Parallelisation
===============

Parallelising by wavenumber (ExoMol)
------------------------------------

Some line lists consist of millions or even billions of transitions. The ExoMol database separates these by wavenumber into multiple ``.trans``-files. pyROX can be executed on these files in parallel and subsequently combine the results into a final output-file.

In this example, we'll show how to run pyROX in parallel on ExoMol's H\ :sub:`2`\ S `line list <https://www.exomol.com/data/molecules/H2S/1H2-32S/AYT2/>`_. The configuration parameters below are also found in the ``examples``-directory on `GitHub <https://github.com/samderegt/pyROX/tree/main/examples/exomol_h2s/exomol_h2s.py>`_.

.. literalinclude:: ../examples/exomol_h2s/exomol_h2s.py
   :language: python

First, we'll need to download the data (~0.9GB). Fortunately, the ``.def.json``-file contains all the necessary information to download the ``.trans``-files so we won't need to provide each URL manually. Navigate to the ``examples/exomol_h2s``-directory and run the following command:

.. code-block:: bash

   cd examples/exomol_h2s
   pyROX exomol_h2s.py -d

Note that the ``transitions`` key in the ``files`` dictionary above is a list of with the filenames of the 34 ``.trans``-files. pyROX would loop over these files and run the calculations sequentially when called like ``pyROX exomol_h2s.py -c``. However, we can provide the ``--files_range`` command-line argument to specify which files to run. The `shell-script <https://github.com/samderegt/pyROX/tree/main/examples/exomol_h2s/run_parallel.sh>`_ ``run_parallel.sh`` provides a convenient way to run 8 pyROX-calls simultaneously.

.. literalinclude:: ../examples/exomol_h2s/run_parallel.sh
   :language: bash

The command-outputs will be logged in the ``logs``-directory. The script should take ~40 minutes to complete for all files (depending on the number of parallel tasks) and can be executed via:

.. code-block:: bash

   sh ./run_parallel.sh &

.. note:: By default, the temporary output-files will be named ``xsec_<transition_filename>.hdf5`` which avoids overwriting the files. 

The temporary output-files can then be combined into a final output-file and plotted with:

.. code-block:: bash

   pyROX exomol_h2s.py -s -p

.. note:: If the temporary or final output-files already exist, pyROX prompts you to overwrite them which can result in an error when running from a shell-script. The command-line argument ``--overwrite`` (or ``-o``) can be used to overwrite the files without prompting.


Parallelising by PT-grid
------------------------

In a similar way, pyROX can run multiple pressure-temperature points in parallel, which can be useful for the HITRAN/HITEMP line lists. The ``examples/hitemp_co``-directory provides a `configuration file <https://github.com/samderegt/pyROX/tree/main/examples/hitemp_co/hitemp_co.py>`_ and a `shell-script <https://github.com/samderegt/pyROX/tree/main/examples/hitemp_co/run_parallel.sh>`_ to run pyROX on 4 temperature points in parallel. 

.. literalinclude:: ../examples/hitemp_co/run_parallel.sh
   :language: bash

.. note:: In the above shell-script we provide new basenames for the output-files (``--tmp_output_basename`` or ``-out``) to avoid overwriting the temporary output-files.

.. note:: HITRAN/HITEMP requires a login to download data, which will raise an error when downloading the ``05_HITEMP2019.par.bz2``-file. Please `download this manually <https://hitran.org/files/HITEMP/bzip2format/05_HITEMP2019.par.bz2>`_ and move it into the ``input_data_dir`` directory. 

The partition function and broadening files can be downloaded with the ``--download`` (or ``-d``) command-line argument. After downloading the data, the shell-script can be executed and will take ~5 minutes to complete for all temperature points. Finally, the temporary output-files can be combined into a final output-file and plotted.

.. code-block:: bash
    
   cd examples/hitemp_co
   pyROX hitemp_co.py -d
   sh ./run_parallel.sh
   pyROX hitemp_co.py -s -p

Similar to the ``--T_grid`` command-line argument, the ``--P_grid`` command-line argument can be used to specify the pressure points. In that way, the user can also run pyROX, in parallel, on multiple pressure points. 