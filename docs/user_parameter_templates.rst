User Parameter Templates
========================

This page shows the simplest way to create your own parameter folder and load it into ``epcsaft``.

Create a template
-----------------

Use ``create_parameter_template(...)`` to generate a dataset folder with the expected layout:

- ``user_options.json``
- ``pure/``
- ``mixed/binary_interaction/``
- ``mixed/rel_perm/``

If you do not pass arguments, the function will prompt for:

- a location where the folder should be created
- a folder name
- the species list to include in the template

Example:

.. code-block:: python

   from epcsaft import create_parameter_template

   template_root = create_parameter_template(
       location=r"C:\Users\Tanner\Documents\my_epcsaft_data",
       folder_name="water_salt_case",
       species=["H2O", "Na+", "Cl-"],
   )

What gets created
-----------------

The generated folder contains blank CSV files in the same shape the loader expects:

- ``pure/<set>.csv`` with one row per species
- ``mixed/binary_interaction/k_ij.csv``
- ``mixed/binary_interaction/l_ij.csv``
- ``mixed/binary_interaction/k_hb_ij.csv``
- ``mixed/rel_perm/parameters.csv``
- ``user_options.json``

You can open those files in a spreadsheet editor or text editor and fill in your parameter values.

For the meaning of ``user_options.json`` and the supported settings, see ``user_options``.

Load the template
-----------------

Once the files are filled in, point ``ePCSAFTMixture.from_dataset(...)`` at the folder:

.. code-block:: python

   import numpy as np
   from epcsaft import ePCSAFTMixture

   mixture = ePCSAFTMixture.from_dataset(
       template_root,
       ["H2O", "Na+", "Cl-"],
       np.asarray([0.98, 0.01, 0.01]),
       298.15,
   )

If you are working from a checkout of this repository, you can also inspect the example folders in ``data/epcsaft_parameters/`` by name. If you pass a path instead, the loader uses that folder directly.

The example folders in ``data/epcsaft_parameters/`` are the same kind of folders you can create for yourself with this helper, but they are source-checkout examples rather than installed-package assets.

Tips
----

- Keep the species order consistent between the template files and the mixture call.
- Fill ``user_options.json`` only with the options you want to override.
- If a component is a solvent with a known name, the template helper will choose a matching pure-file name such as ``water.csv``.
