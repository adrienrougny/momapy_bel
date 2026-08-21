"""Smoke tests for the `momapy_bel.core` model shape.

These pin the properties the BEL reader relies on: value equality of frozen
model elements, hashability of annotations and definitions, the builder round
trip, and the fact that a nested frozenset cannot live on a model element.
"""

import dataclasses

import frozendict
import pytest

import momapy.builder
import momapy.core

import momapy_bel.core


def _protein(identifier):
    return momapy_bel.core.ProteinAbundance(namespace="HGNC", identifier=identifier)


def _statement():
    return momapy_bel.core.Increases(
        source=_protein("APP"), target=_protein("BACE1")
    )


def test_generic_annotation_hashes_without_definition():
    annotation = momapy_bel.core.BELGenericAnnotation(
        name="Citation", args=("PubMed", "12345")
    )
    assert annotation.definition is None
    assert hash(annotation) == hash(
        momapy_bel.core.BELGenericAnnotation(
            name="Citation", args=("PubMed", "12345")
        )
    )


def test_generic_annotation_hashes_with_definition():
    definition = momapy_bel.core.BELGenericAnnotationDefinition(
        name="Species",
        as_type=momapy_bel.core.BELAsDefinitionType.URL,
        as_=("http://example.org/species.belanno",),
    )
    annotation = momapy_bel.core.BELGenericAnnotation(
        name="Species", definition=definition, args=("9606",)
    )
    assert hash(annotation)
    assert annotation.name == "Species"


def test_definition_with_tuple_as_hashes():
    definition = momapy_bel.core.BELGenericAnnotationDefinition(
        name="Confidence",
        as_type=momapy_bel.core.BELAsDefinitionType.LIST,
        as_=("High", "Medium", "Low"),
    )
    assert hash(definition)
    namespace_definition = momapy_bel.core.BELNamespaceDefinition(
        name="Custom",
        as_type=momapy_bel.core.BELAsDefinitionType.LIST,
        as_=("A", "B"),
    )
    assert hash(namespace_definition)


def test_annotation_table_is_keyed_by_value():
    contexts = frozenset(
        [
            frozenset(
                [momapy_bel.core.BELGenericAnnotation(name="Species", args=("9606",))]
            ),
            frozenset(
                [momapy_bel.core.BELGenericAnnotation(name="Species", args=("10090",))]
            ),
        ]
    )
    table = frozendict.frozendict({_statement(): contexts})
    # An equal but distinct statement finds the same entry: `id_` does not
    # take part in equality or hashing.
    other = _statement()
    assert other is not next(iter(table))
    assert table[other] == contexts


def test_builder_round_trip():
    model = momapy_bel.core.BELModel(
        statements=frozenset([_statement()]),
        namespace_definitions=frozenset(
            [
                momapy_bel.core.BELNamespaceDefinition(
                    name="HGNC",
                    as_type=momapy_bel.core.BELAsDefinitionType.URL,
                    as_=("hgnc.belns",),
                )
            ]
        ),
    )
    assert (
        momapy.builder.object_from_builder(momapy.builder.builder_from_object(model))
        == model
    )


def test_namespace_definition_order_does_not_matter():
    first = momapy_bel.core.BELNamespaceDefinition(
        name="HGNC",
        as_type=momapy_bel.core.BELAsDefinitionType.URL,
        as_=("hgnc.belns",),
    )
    second = momapy_bel.core.BELNamespaceDefinition(
        name="CHEBI",
        as_type=momapy_bel.core.BELAsDefinitionType.URL,
        as_=("chebi.belns",),
    )
    one = momapy_bel.core.BELModel(namespace_definitions=frozenset([first, second]))
    two = momapy_bel.core.BELModel(namespace_definitions=frozenset([second, first]))
    assert one == two


def test_is_submodel_is_directional():
    small = momapy_bel.core.BELModel(statements=frozenset([_statement()]))
    large = momapy_bel.core.BELModel(
        statements=frozenset([_statement(), _protein("MAPT")]),
        namespace_definitions=frozenset(
            [
                momapy_bel.core.BELNamespaceDefinition(
                    name="HGNC",
                    as_type=momapy_bel.core.BELAsDefinitionType.URL,
                    as_=("hgnc.belns",),
                )
            ]
        ),
    )
    assert small.is_submodel(large)
    assert not large.is_submodel(small)


def test_complex_members_are_order_insensitive():
    one = momapy_bel.core.ComplexAbundance(
        members=frozenset([_protein("APP"), _protein("BACE1")])
    )
    two = momapy_bel.core.ComplexAbundance(
        members=frozenset([_protein("BACE1"), _protein("APP")])
    )
    assert one == two


def test_complex_members_deduplicate():
    homodimer = momapy_bel.core.ComplexAbundance(
        members=frozenset([_protein("NFKB1"), _protein("NFKB1")])
    )
    assert len(homodimer.members) == 1


def test_nested_frozenset_field_breaks_the_builder():
    """A `frozenset[frozenset[...]]` field cannot live on a model element.

    `momapy.builder` rebuilds collection fields eagerly, which evaluates
    `set([set(...)])` and raises. This is why the annotation contexts live in
    the `ReaderResult` rather than on `BELModel` or on a statement.
    """

    @dataclasses.dataclass(frozen=True, kw_only=True)
    class _Nested(momapy_bel.core.BELModelElement):
        contexts: frozenset = dataclasses.field(default_factory=frozenset)

    element = _Nested(contexts=frozenset([frozenset(["a"])]))
    with pytest.raises(TypeError):
        momapy.builder.builder_from_object(element)


def _every_core_class():
    for name in dir(momapy_bel.core):
        obj = getattr(momapy_bel.core, name)
        if (
            isinstance(obj, type)
            and dataclasses.is_dataclass(obj)
            and obj.__module__ == "momapy_bel.core"
            and not momapy.builder.issubclass_or_builder(obj, momapy.builder.Builder)
        ):
            yield name, obj


def test_every_core_class_instantiates_and_builds():
    """`get_or_make_builder_cls` is lazy, so every class must be exercised."""
    protein = _protein("APP")
    fusion = momapy_bel.core.Fusion(
        namespace5="HGNC",
        identifier5="TMPRSS2",
        range5="r.1_79",
        namespace3="HGNC",
        identifier3="ERG",
        range3="r.312_5034",
    )
    values = {
        "namespace": "HGNC",
        "identifier": "APP",
        "abundance": protein,
        "source": protein,
        "target": protein,
        "start_stop": "5_20",
        "descriptor": "p.Tyr633X",
        "fusion": fusion,
        "namespace5": "HGNC",
        "identifier5": "TMPRSS2",
        "range5": "r.1_79",
        "namespace3": "HGNC",
        "identifier3": "ERG",
        "range3": "r.312_5034",
        "from_namespace": "GO",
        "from_identifier": "intracellular",
        "to_namespace": "GO",
        "to_identifier": "extracellular space",
        "name": "HGNC",
        "as_type": momapy_bel.core.BELAsDefinitionType.URL,
        "as_": ("hgnc.belns",),
        "args": ("PubMed", "12345"),
    }
    seen = 0
    for name, cls in _every_core_class():
        kwargs = {}
        for field in dataclasses.fields(cls):
            if field.name == "id_":
                continue
            if (
                field.default is dataclasses.MISSING
                and field.default_factory is dataclasses.MISSING
            ):
                assert field.name in values, f"no test value for {name}.{field.name}"
                kwargs[field.name] = values[field.name]
        instance = cls(**kwargs)
        assert hash(instance) == hash(cls(**kwargs))
        if isinstance(instance, (momapy.core.ModelElement, momapy.core.Model)):
            assert (
                momapy.builder.object_from_builder(
                    momapy.builder.builder_from_object(instance)
                )
                == instance
            )
        seen += 1
    assert seen > 40
