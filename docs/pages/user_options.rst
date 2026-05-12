User Options
============

The ``user_options.json`` file lets you change the model settings used for a dataset.
It is optional, and it can stay very small. Put only the settings you want to change.

How it works
------------

When a dataset is loaded, the package starts with its default settings and then
overlays anything you put in ``user_options.json``. That means:

- you can leave out any option you do not want to touch
- you only need to write the values you want to override
- nested settings stay nested, so you can update just one part of the model

A minimal file looks like this:

.. code-block:: json

   {
     "elec_model": {
       "include_born_model": false
     }
   }

Supported top-level options
---------------------------

- ``solvated_ion_diameter_mixing_rule``: enable or disable the ion-diameter mixing rule for solvated ions
- ``ion_dispersion_mixing_rule``: enable or disable the ion-dispersion mixing rule
- ``elec_model``: nested electrical-model settings

Legacy top-level electrolyte shorthands are no longer accepted. Put all dielectric,
Debye-Hückel, and Born configuration under ``elec_model`` using the nested blocks below.

Electrostatic model options
---------------------------

Use ``elec_model`` to adjust the dielectric, Debye-Hückel, and Born settings.
Only the documented nested keys are supported.

``rel_perm``
~~~~~~~~~~~~~

- ``rule``: dielectric-mixing rule
  - ``1`` or ``linear``
  - ``0`` or ``constant``
  - ``2`` or ``linear-massfraction``
  - ``3`` or ``combined``
  - ``4`` or ``empirical``
  - ``7`` or ``linear-saltfraction``
  - ``8`` or ``aqueous-organic``
  - ``9`` or ``salt-free-massfraction`` for a solvent-only mass-fraction average that ignores ions in the dielectric average
- ``differential_mode``: choose ``auto`` (default), ``analytical``, or ``autodiff``

``hc_model``, ``disp_model``, ``assoc_model``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``dadx_differential_mode``: choose ``auto`` (default), ``analytical``, or ``autodiff``

``DH_model``
~~~~~~~~~~~~

- ``d_ion_mode``: choose ``0``, ``1``, or ``2``
- ``bjeruum_treatment``: enable or disable Bjerrum treatment
- ``mu_DH_model.differential_mode``: choose ``auto`` (default), ``analytical``, or ``autodiff``
- ``mu_DH_model.comp_dep_rel_perm``: include composition-dependent permittivity in the derivative model
- ``mu_DH_model.include_sum_term``: include the sum term in the derivative model

``born_model``
~~~~~~~~~~~~~~

- ``d_Born_mode``: choose ``0``, ``1``, ``2``, ``3``, or a named alias such as ``fitted_param``
- ``solvation_shell_model``: enable or disable the solvation-shell model
- ``dielectric_saturation``: enable or disable dielectric saturation
- ``bulk_mode``: choose ``mix`` or ``solvent``
- ``mu_born_model.differential_mode``: choose ``auto`` (default), ``analytical``, or ``autodiff``
- ``mu_born_model.comp_dep_rel_perm``: include composition-dependent permittivity in the derivative model
- ``mu_born_model.include_sum_term``: include the sum term in the derivative model
- ``mu_born_model.comp_dep_delta_d``: include the delta-d composition term

Common examples
---------------

Turn off the Born contribution:

.. code-block:: json

   {
     "elec_model": {
       "include_born_model": false
     }
   }

Use constant dielectric mixing:

.. code-block:: json

   {
     "elec_model": {
       "rel_perm": {
         "rule": "constant"
       }
     }
   }

Use a more advanced Born setup:

.. code-block:: json

   {
     "elec_model": {
       "rel_perm": {
         "rule": "empirical",
          "differential_mode": "auto"
       },
       "born_model": {
         "d_Born_mode": 3,
         "solvation_shell_model": true,
         "dielectric_saturation": true,
         "mu_born_model": {
            "differential_mode": "auto",
           "comp_dep_delta_d": true
         }
       }
     }
   }

Tips
----

- Start from a small file and add only the settings you actually need.
- If a value is left out, the package keeps its default.
- The default derivative policy is ``auto``: use validated analytical derivatives where they already exist, then autodiff for implemented derivative paths. Unsupported paths raise clear ``backend_unavailable`` errors.
- ``autodiff`` is stricter than ``auto``. It requests the autodiff backend directly and raises for unsupported paths instead of selecting an analytical formula.
- If you are unsure, keep ``user_options.json`` empty and add settings one at a time.
