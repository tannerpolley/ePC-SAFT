Parameter Schema
================

The public package still accepts legacy parameter dictionaries through
``ePCSAFTMixture.from_params(...)``. New code can also use canonical records:

* ``ComponentIdentifier``
* ``AssociationSite``
* ``PureRecord``
* ``BinaryRecord``
* ``PermittivityRecord``
* ``ParameterSet``

``ParameterSet`` can be constructed from records, legacy dictionaries, or
packaged datasets, then converted back to the native legacy payload.
Canonical ``PureRecord.molar_mass`` values are always kg/mol, matching the
legacy native ``MW`` field. If a source table reports g/mol, use
``PureRecord.from_g_per_mol(...)`` so the conversion is explicit.

.. code-block:: python

   import epcsaft

   params = epcsaft.ParameterSet.from_records(
       [
           epcsaft.PureRecord("H2O", molar_mass=0.01801528, m=1.2047, sigma=2.7927, epsilon_k=353.95),
           epcsaft.PureRecord.from_g_per_mol("Ethanol", molar_mass_g_per_mol=46.07, m=2.0, sigma=4.0, epsilon_k=250.0),
       ],
       [epcsaft.BinaryRecord(("H2O", "Ethanol"), k_ij=0.02)],
   )

   mixture = epcsaft.ePCSAFTMixture.from_params(params)

Record Validation
-----------------

``ParameterSet.validate()`` checks required pure records, duplicate component
labels, positive pure-component values, kg/mol molar-mass scale, charged-species
Born diameters, and binary records that reference unknown species.

Legacy Compatibility
--------------------

``ParameterSet.to_legacy_dict()`` emits the array dictionary consumed by the
native runtime. Existing dictionary payloads keep working:

.. code-block:: python

   mixture = epcsaft.ePCSAFTMixture.from_params(legacy_dict, species=["A", "B"])

Use canonical records when provenance, component identity, association sites,
or binary-interaction records need to be explicit.

Template Generation
-------------------

``create_parameter_template(...)`` defaults to the legacy CSV layout so existing
``ePCSAFTMixture.from_dataset(...)`` workflows keep working:

.. code-block:: python

   root = epcsaft.create_parameter_template(
       "data/reference/epcsaft_parameters",
       "2026_User",
       ["H2O", "Na+", "Cl-"],
   )

New datasets can request a canonical JSON scaffold instead:

.. code-block:: python

   root = epcsaft.create_parameter_template(
       "data/reference/epcsaft_parameters",
       "2026_User_Canonical",
       ["H2O", "Na+", "Cl-"],
       schema="canonical",
   )

The canonical scaffold writes ``parameter_set.json`` plus ``user_options.json``.
Its pure records use the same names as ``PureRecord`` and include
``molar_mass_units: "kg/mol"`` so source-table g/mol values are converted before
they become runtime parameters.
