"""Round-trip tests for the BEL reader and writer.

`read -> write -> read` compares models and annotation tables, not bytes:
the writer normalises `reaction` to `rxn` and `--` to `association`, and sorts
set-valued content. `write -> read -> write` is the byte-level direction, on a
hand-built model.
"""

import pathlib

import pytest

import momapy_bel.core
import momapy_bel.io.bel


DATA_DIR = pathlib.Path(__file__).parent.parent / "data"


def _read(path):
    return momapy_bel.io.bel.BELReader.read(path, return_type="model")


def _write(result, path):
    return momapy_bel.io.bel.BELWriter.write(
        result.obj,
        path,
        element_to_annotations=result.element_to_annotations,
        annotation_definitions=result.annotation_definitions,
        # A faithful round trip must re-emit every statement line, including
        # the term-only ones the writer drops by default.
        with_abundances_as_statements=True,
        with_biological_processes_as_statements=True,
    )


def _assert_read_write_read(source_path, tmp_path):
    first = _read(source_path)
    written = tmp_path / "written.bel"
    _write(first, written)
    second = _read(written)
    assert first.obj == second.obj
    assert first.element_to_annotations == second.element_to_annotations
    assert first.annotation_definitions == second.annotation_definitions
    return first, second


@pytest.mark.parametrize(
    "name",
    [
        "glued.bel",
        "set_forms.bel",
        "definitions.bel",
        "fragments.bel",
        "variants.bel",
        "quoted_identifiers.bel",
        "terms_only.bel",
        "containers.bel",
        "unset_forms.bel",
        "support_evidence.bel",
        "document.bel",
        "relations.bel",
    ],
)
def test_read_write_read(name, tmp_path):
    _assert_read_write_read(DATA_DIR / name, tmp_path)


def _hand_built_model():
    app = momapy_bel.core.ProteinAbundance(namespace="HGNC", identifier="APP")
    bace1 = momapy_bel.core.ProteinAbundance(
        namespace="HGNC",
        identifier="BACE1",
        variants=(momapy_bel.core.Variant(descriptor="P,S,396"),),
        fragment=momapy_bel.core.Fragment(start_stop="?", descriptor="41kD"),
    )
    disease = momapy_bel.core.Pathology(
        namespace="MESH", identifier="Alzheimer Disease"
    )
    reaction = momapy_bel.core.Reaction(
        reactants=frozenset(
            [momapy_bel.core.Abundance(namespace="CHEBI", identifier="ATP")]
        ),
        products=frozenset(
            [momapy_bel.core.Abundance(namespace="CHEBI", identifier="ADP")]
        ),
    )
    translocation = momapy_bel.core.Translocation(
        abundance=app,
        from_namespace="GO",
        from_identifier="intracellular",
        to_namespace="GO",
        to_identifier="extracellular space",
    )
    statements = frozenset(
        [
            momapy_bel.core.Increases(source=app, target=bace1),
            momapy_bel.core.Association(source=bace1, target=disease),
            momapy_bel.core.EquivalentTo(source=app, target=bace1),
            reaction,
            translocation,
        ]
    )
    model = momapy_bel.core.BELModel(
        statements=statements,
        namespace_definitions=frozenset(
            [
                momapy_bel.core.BELNamespaceDefinition(
                    name="HGNC",
                    as_type=momapy_bel.core.BELAsDefinitionType.URL,
                    as_=("https://example.org/hgnc.belns",),
                ),
                momapy_bel.core.BELNamespaceDefinition(
                    name="DBSNP",
                    as_type=momapy_bel.core.BELAsDefinitionType.PATTERN,
                    as_=(r"^rs\d+$",),
                ),
            ]
        ),
    )
    citation = momapy_bel.core.BELGenericAnnotation(
        name="Citation", args=("PubMed", "12345")
    )
    other_citation = momapy_bel.core.BELGenericAnnotation(
        name="Citation", args=("PubMed", "67890")
    )
    confidence = momapy_bel.core.BELGenericAnnotation(
        name="Confidence", args=("High",)
    )
    element_to_annotations = {
        # The same claim asserted under two distinct contexts.
        momapy_bel.core.Increases(source=app, target=bace1): frozenset(
            [
                frozenset([citation, confidence]),
                frozenset([other_citation]),
            ]
        ),
        momapy_bel.core.Association(source=bace1, target=disease): frozenset(
            [frozenset([citation])]
        ),
        model: frozenset(
            [
                frozenset(
                    [
                        momapy_bel.core.BELDocumentAnnotation(
                            name="Hand built", version="1.0"
                        )
                    ]
                )
            ]
        ),
    }
    annotation_definitions = frozenset(
        [
            momapy_bel.core.BELGenericAnnotationDefinition(
                name="Confidence",
                as_type=momapy_bel.core.BELAsDefinitionType.LIST,
                as_=("High", "Medium", "Low"),
            )
        ]
    )
    return model, element_to_annotations, annotation_definitions


def test_write_read_write_is_byte_stable(tmp_path):
    model, element_to_annotations, annotation_definitions = _hand_built_model()
    first_path = tmp_path / "first.bel"
    momapy_bel.io.bel.BELWriter.write(
        model,
        first_path,
        element_to_annotations=element_to_annotations,
        annotation_definitions=annotation_definitions,
        with_abundances_as_statements=True,
        with_biological_processes_as_statements=True,
    )
    result = _read(first_path)
    assert result.obj == model
    second_path = tmp_path / "second.bel"
    _write(result, second_path)
    assert first_path.read_text() == second_path.read_text()


def test_url_and_pattern_survive_a_round_trip(tmp_path):
    source = tmp_path / "source.bel"
    source.write_text(
        'DEFINE NAMESPACE HGNC AS URL "hgnc.belns"\n'
        'DEFINE NAMESPACE GO AS PATTERN "GO:\\d+"\n'
        "p(HGNC:APP)\n"
    )
    written = tmp_path / "written.bel"
    _write(_read(source), written)
    text = written.read_text()
    assert 'AS URL "hgnc.belns"' in text
    assert 'AS PATTERN "GO:\\d+"' in text


def test_write_returns_a_writer_result(tmp_path):
    model, _, _ = _hand_built_model()
    path = tmp_path / "out.bel"
    result = momapy_bel.io.bel.BELWriter.write(model, path)
    assert result.obj is model
    assert result.file_path == path

