Parameter Schema
================

The public package still accepts native array parameter dictionaries through
``ePCSAFTMixture.from_params(...)``. New code can also use canonical records:

* ``ComponentIdentifier``
* ``AssociationSite``
* ``PureRecord``
* ``BinaryRecord``
* ``PermittivityRecord``
* ``ParameterSet``
* ``ParameterSource``

``ParameterSet`` can be constructed from records, native array dictionaries, or
packaged datasets, then converted back to the native array payload.
Canonical ``PureRecord.molar_mass`` values are always kg/mol, matching the
native ``MW`` field. If a source table reports g/mol, use
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

Dictionary Compatibility
------------------------

``ParameterSet.to_runtime_dict()`` is the canonical emitter for the array
dictionary consumed by the native runtime. ``ParameterSet.to_legacy_dict()`` and
``ParameterSet.from_dict(...)`` remain compatibility APIs for existing callers,
but new code should use the runtime/canonical names:

.. code-block:: python

   runtime_dict = params.to_runtime_dict()
   mixture = epcsaft.ePCSAFTMixture.from_params(runtime_dict, species=["A", "B"])

Runtime options are copied into the emitted payload only after the canonical
parameter arrays and matrices are built. They cannot override parameter payload
keys such as ``m``, ``k_ij``, ``k_hb``, or structural keys such as
``components``. Regression and reactive-regression helpers that accept an
in-memory ``ParameterSet`` use ``ParameterSet.to_runtime_dict()`` as the single
runtime-payload boundary.

``ParameterSource`` is the resolver behind dataset paths, direct
``ParameterSet`` objects, and direct runtime dictionaries. It keeps source labels
stable for diagnostics and applies per-call runtime options at the same payload
boundary:

.. code-block:: python

   source = epcsaft.ParameterSource(params, species=["A", "B"])
   runtime_dict = source.to_runtime_dict(user_options={"source_tag": "train"})

Existing dictionary payloads keep working:

.. code-block:: python

   mixture = epcsaft.ePCSAFTMixture.from_params(legacy_dict, species=["A", "B"])

Use canonical records when provenance, component identity, association sites,
or binary-interaction records need to be explicit.

Template Generation
-------------------

``create_parameter_template(...)`` defaults to the native CSV layout so existing
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
they become runtime parameters. Once the JSON records are filled in, the folder
is a normal canonical dataset path:

.. code-block:: python

   params = epcsaft.ParameterSet.from_dataset(
       root,
       ["H2O", "Na+", "Cl-"],
       [0.98, 0.01, 0.01],
       298.15,
   )
   mixture = epcsaft.ePCSAFTMixture.from_dataset(
       root,
       ["H2O", "Na+", "Cl-"],
       [0.98, 0.01, 0.01],
       298.15,
   )
