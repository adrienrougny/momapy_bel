"""Reader tests over the hand-written corpus in `tests/data`.

One file per hazard, kept small so a failure points at a single construct.
"""

import dataclasses
import pathlib

import pytest

import momapy.io
import momapy.io.core
import momapy_bel.core
import momapy_bel.io.bel


DATA_DIR = pathlib.Path(__file__).parent.parent / "data"


def _read(name, **options):
    return momapy_bel.io.bel.BELReader.read(DATA_DIR / name, **options)


def _model(name, **options):
    return _read(name, return_type="model", **options).obj


def _protein(identifier, **kwargs):
    return momapy_bel.core.ProteinAbundance(
        namespace="HGNC", identifier=identifier, **kwargs
    )


def _annotation_names(result, element):
    return {
        annotation.name
        for context in result.element_to_annotations[element]
        for annotation in context
    }


@pytest.mark.parametrize("return_type", ["map", "model"])
def test_return_type(return_type):
    result = _read("glued.bel", return_type=return_type)
    if return_type == "map":
        assert isinstance(result.obj, momapy_bel.core.BELMap)
        assert result.obj.layout is None
        model = result.obj.model
    else:
        assert isinstance(result.obj, momapy_bel.core.BELModel)
        model = result.obj
    assert len(model.statements) == 3


def test_return_type_layout_is_not_supported():
    with pytest.raises(NotImplementedError):
        _read("glued.bel", return_type="layout")


def test_return_type_is_validated():
    with pytest.raises(ValueError):
        _read("glued.bel", return_type="nonsense")


def test_check_file():
    assert momapy_bel.io.bel.BELReader.check_file(DATA_DIR / "definitions.bel")
    assert momapy_bel.io.bel.BELReader.check_file(DATA_DIR / "document.bel")
    assert not momapy_bel.io.bel.BELReader.check_file(DATA_DIR / "relations.bel")
    assert not momapy_bel.io.bel.BELReader.check_file(DATA_DIR / "missing.bel")


def test_glued_tokens_and_relations():
    model = _model("glued.bel")
    assert (
        momapy_bel.core.IsA(
            source=momapy_bel.core.Abundance(
                namespace="CHEBI", identifier="methoctramine"
            ),
            target=momapy_bel.core.Abundance(
                namespace="CHEBI", identifier="antagonist"
            ),
        )
        in model.statements
    )
    assert (
        momapy_bel.core.Increases(source=_protein("APP"), target=_protein("BACE1"))
        in model.statements
    )
    assert model.namespace_definitions == frozenset(
        [
            momapy_bel.core.BELNamespaceDefinition(
                name="HGNC",
                as_type=momapy_bel.core.BELAsDefinitionType.URL,
                as_=("https://example.org/hgnc.belns",),
            )
        ]
    )


def test_set_forms():
    result = _read("set_forms.bel", return_type="model")
    statement = next(iter(result.obj.statements))
    assert result.element_to_annotations[statement] == frozenset(
        [
            frozenset(
                [
                    momapy_bel.core.BELGenericAnnotation(
                        name="Confidence", args=("Low",)
                    ),
                    momapy_bel.core.BELGenericAnnotation(
                        name="Species", args=("9606",)
                    ),
                ]
            )
        ]
    )


def test_definitions():
    result = _read("definitions.bel", return_type="model")
    assert result.obj.namespace_definitions == frozenset(
        [
            momapy_bel.core.BELNamespaceDefinition(
                name="DBSNP", as_type=momapy_bel.core.BELAsDefinitionType.PATTERN, as_=(r"^rs\d+$",)
            ),
            momapy_bel.core.BELNamespaceDefinition(
                name="HGNC",
                as_type=momapy_bel.core.BELAsDefinitionType.URL,
                as_=("https://example.org/hgnc.belns",),
            ),
        ]
    )
    # Declared but never SET: the reader result carries it anyway.
    assert result.annotation_definitions == frozenset(
        [
            momapy_bel.core.BELGenericAnnotationDefinition(
                name="Confidence",
                as_type=momapy_bel.core.BELAsDefinitionType.LIST,
                as_=("High", "Medium", "Low"),
            ),
            momapy_bel.core.BELGenericAnnotationDefinition(
                name="Species",
                as_type=momapy_bel.core.BELAsDefinitionType.URL,
                as_=("https://example.org/species.belanno",),
            ),
        ]
    )


def test_fragments():
    model = _model("fragments.bel")
    assert (
        _protein(
            "APP",
            fragment=momapy_bel.core.Fragment(start_stop="?", descriptor="41kD"),
        )
        in model.statements
    )
    assert (
        _protein(
            "MAPT",
            fragment=momapy_bel.core.Fragment(
                start_stop="672_713", descriptor=None
            ),
        )
        in model.statements
    )
    assert len(model.statements) == 3


def test_two_variants_are_both_kept():
    model = _model("variants.bel")
    protein = next(
        statement
        for statement in model.statements
        if getattr(statement, "identifier", None) == "App"
    )
    assert protein.variants == (
        momapy_bel.core.Variant(descriptor="K,670,N"),
        momapy_bel.core.Variant(descriptor="M,671,L"),
    )


def test_protein_modification_arguments_are_positional():
    model = _model("variants.bel")
    assert (
        _protein(
            "SNCA",
            modifications=(
                momapy_bel.core.ProteinModification(
                    namespace="", identifier="Ph", amino_acid="P", residue="213"
                ),
            ),
        )
        in model.statements
    )
    assert (
        _protein(
            "TAU",
            modifications=(
                momapy_bel.core.ProteinModification(
                    namespace="MOD", identifier="PhosphoSerine"
                ),
            ),
        )
        in model.statements
    )


def test_quoted_identifiers():
    model = _model("quoted_identifiers.bel")
    identifiers = {statement.identifier for statement in model.statements}
    assert identifiers == {
        "sphingomyelin 34:2",
        "1,2-dimethoxybenzene",
        "dotarizine",
        "PPARG",
    }


def test_term_only_lines_and_member_deduplication():
    model = _model("terms_only.bel")
    assert len(model.statements) == 3
    composite = next(
        statement
        for statement in model.statements
        if isinstance(statement, momapy_bel.core.CompositeAbundance)
    )
    # WASF2 is listed twice and collapses.
    assert len(composite.members) == 2


def test_containers():
    model = _model("containers.bel")
    has_members = next(
        statement
        for statement in model.statements
        if isinstance(statement, momapy_bel.core.HasMembers)
    )
    assert isinstance(has_members.target, momapy_bel.core.List)
    assert len(has_members.target.elements) == 2
    # `rxn` and `reaction` are the same claim, so they collapse.
    reactions = [
        statement
        for statement in model.statements
        if isinstance(statement, momapy_bel.core.Reaction)
    ]
    assert len(reactions) == 1
    translocation = next(
        statement
        for statement in model.statements
        if isinstance(statement, momapy_bel.core.Translocation)
    )
    assert translocation.from_identifier == "intracellular"
    assert translocation.to_identifier == "extracellular space"


def test_unset_forms():
    result = _read("unset_forms.bel", return_type="model")
    first = momapy_bel.core.Increases(
        source=_protein("APP"), target=_protein("BACE1")
    )
    second = momapy_bel.core.Increases(
        source=_protein("APP"), target=_protein("MAPT")
    )
    third = _protein("MAPT")
    assert _annotation_names(result, first) == {
        "Confidence",
        "Species",
        "Citation",
    }
    assert _annotation_names(result, second) == {"Citation"}
    assert result.element_to_annotations[third] == frozenset([frozenset()])


def test_support_and_evidence_share_a_slot():
    result = _read("support_evidence.bel", return_type="model")
    first = momapy_bel.core.Increases(
        source=_protein("APP"), target=_protein("BACE1")
    )
    second = momapy_bel.core.Increases(
        source=_protein("APP"), target=_protein("MAPT")
    )
    assert result.element_to_annotations[first] == frozenset(
        [
            frozenset(
                [
                    momapy_bel.core.BELGenericAnnotation(
                        name="Support", args=("some evidence",)
                    )
                ]
            )
        ]
    )
    # `UNSET Evidence` clears what `SET Support` set, and `UNSET NeverSet` is
    # not an error.
    assert result.element_to_annotations[second] == frozenset([frozenset()])


def test_document_annotations_comments_and_continuations():
    result = _read("document.bel", return_type="model")
    assert result.element_to_annotations[result.obj] == frozenset(
        [
            frozenset(
                [
                    momapy_bel.core.BELDocumentAnnotation(
                        name="Test Document",
                        version="1.0",
                        contact_info="a@example.org",
                    )
                ]
            )
        ]
    )
    # The comment mentioning `p(HGNC:X)` is not a statement, and the joined
    # continuation is one statement, not two.
    assert len(result.obj.statements) == 1
    statement = next(iter(result.obj.statements))
    support = next(
        annotation
        for context in result.element_to_annotations[statement]
        for annotation in context
    )
    assert support.args == ("set of related networks     spanning two lines",)


def test_relation_spellings():
    model = _model("relations.bel")
    assert {type(statement) for statement in model.statements} == {
        momapy_bel.core.EquivalentTo,
        momapy_bel.core.NoCorrelation,
        momapy_bel.core.RateLimitingStepFor,
        momapy_bel.core.TranslatedTo,
        momapy_bel.core.TranscribedTo,
    }
    # `eq` and `equivalentTo` are the same claim and collapse.
    assert len(model.statements) == 5


def test_unknown_function_is_reported():
    with pytest.raises(momapy_bel.io.bel.BELReadingError) as excinfo:
        _model("bad_function.bel")
    assert "gmod" in str(excinfo.value)
    assert excinfo.value.failures[0][0] == 1


def test_with_annotations_false():
    result = _read("unset_forms.bel", return_type="model", with_annotations=False)
    assert result.element_to_annotations is None
    assert result.annotation_definitions is None
    assert len(result.obj.statements) == 3


def test_with_model_false():
    result = _read("glued.bel", with_model=False)
    assert result.obj.model is None


def test_registered_reader_round_trips_through_momapy_io():
    result = momapy.io.read(
        DATA_DIR / "glued.bel", reader="bel", return_type="model"
    )
    assert isinstance(result.obj, momapy_bel.core.BELModel)


def test_reader_result_is_a_reader_result_subclass():
    result = _read("definitions.bel", return_type="model")
    assert isinstance(result, momapy_bel.io.bel.BELReaderResult)
    assert isinstance(result, momapy.io.core.ReaderResult)
    # The field is declared, not attached: it shows up in the dataclass field
    # list, so anything reflecting over the result sees it.
    assert "annotation_definitions" in {
        field.name for field in dataclasses.fields(result)
    }
