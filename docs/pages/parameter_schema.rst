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

.. code-block:: python

   import epcsaft

   params = epcsaft.ParameterSet.from_records(
       [
           epcsaft.PureRecord("A", molar_mass=18.0, m=1.0, sigma=3.0, epsilon_k=200.0),
           epcsaft.PureRecord("B", molar_mass=46.0, m=2.0, sigma=4.0, epsilon_k=250.0),
       ],
       [epcsaft.BinaryRecord(("A", "B"), k_ij=0.02)],
   )

   mixture = epcsaft.ePCSAFTMixture.from_params(params)

Record Validation
-----------------

``ParameterSet.validate()`` checks required pure records, duplicate component
labels, positive pure-component values, charged-species Born diameters, and
binary records that reference unknown species.

Legacy Compatibility
--------------------

``ParameterSet.to_legacy_dict()`` emits the array dictionary consumed by the
native runtime. Existing dictionary payloads keep working:

.. code-block:: python

   mixture = epcsaft.ePCSAFTMixture.from_params(legacy_dict, species=["A", "B"])

Use canonical records when provenance, component identity, association sites,
or binary-interaction records need to be explicit.
