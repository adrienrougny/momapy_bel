import dataclasses

import momapy.core
import momapy.core.map
import momapy.builder


@dataclasses.dataclass(frozen=True, kw_only=True)
class BELModelElement(momapy.core.ModelElement):
    pass


@dataclasses.dataclass(frozen=True, kw_only=True)
class List(BELModelElement):
    elements: tuple[BELModelElement] = dataclasses.field(default_factory=tuple)


@dataclasses.dataclass(frozen=True, kw_only=True)
class Location(BELModelElement):
    namespace: str
    identifier: str


@dataclasses.dataclass(frozen=True, kw_only=True)
class Abundance(BELModelElement):
    namespace: str
    identifier: str
    location: Location | None = None


@dataclasses.dataclass(frozen=True, kw_only=True)
class MolecularActivity(BELModelElement):
    namespace: str
    identifier: str


@dataclasses.dataclass(frozen=True, kw_only=True)
class Activity(BELModelElement):
    abundance: Abundance
    molecular_activity: MolecularActivity | None = None


@dataclasses.dataclass(frozen=True, kw_only=True)
class BiologicalProcess(BELModelElement):
    namespace: str
    identifier: str


@dataclasses.dataclass(frozen=True, kw_only=True)
class CellSecretion(BELModelElement):
    abundance: Abundance


@dataclasses.dataclass(frozen=True, kw_only=True)
class CellSurfaceExpression(BELModelElement):
    abundance: Abundance


@dataclasses.dataclass(frozen=True, kw_only=True)
class ComplexAbundance(Abundance):
    """A named or member-defined complex abundance.

    `members` is a `frozenset`, so a complex listing the same member twice
    (a homodimer such as `complex(p(HGNC:NFKB1), p(HGNC:NFKB1))`) collapses to
    a single member. BEL 2.1.3 defines no stoichiometry and says nothing about
    member uniqueness, so nothing distinguishes the collapsed form from the
    written one; the multiplicity is lost. Revisit only if BEL gains
    stoichiometry.
    """

    namespace: str | None = None
    identifier: str | None = None
    members: frozenset[Abundance] = dataclasses.field(
        default_factory=frozenset
    )


@dataclasses.dataclass(frozen=True, kw_only=True)
class CompositeAbundance(Abundance):
    """A named or member-defined composite abundance.

    `members` is a `frozenset`, so a composite listing the same member twice
    collapses to a single member; see `ComplexAbundance` for the rationale.
    """

    namespace: str | None = None
    identifier: str | None = None
    members: frozenset[Abundance] = dataclasses.field(
        default_factory=frozenset
    )


@dataclasses.dataclass(frozen=True, kw_only=True)
class Degradation(BELModelElement):
    abundance: Abundance


@dataclasses.dataclass(frozen=True, kw_only=True)
class Fragment(BELModelElement):
    start_stop: str
    descriptor: str | None


@dataclasses.dataclass(frozen=True, kw_only=True)
class Fusion(BELModelElement):
    namespace5: str
    identifier5: str
    range5: str
    namespace3: str
    identifier3: str
    range3: str


@dataclasses.dataclass(frozen=True, kw_only=True)
class Variant(BELModelElement):
    descriptor: str


@dataclasses.dataclass(frozen=True, kw_only=True)
class GeneAbundance(Abundance):
    namespace: str | None
    identifier: str | None
    fusion: Fusion | None
    variants: tuple[Variant] = dataclasses.field(default_factory=tuple)


@dataclasses.dataclass(frozen=True, kw_only=True)
class MicroRNAAbundance(Abundance):
    namespace: str | None
    identifier: str | None
    fusion: Fusion | None
    variants: tuple[Variant] = dataclasses.field(default_factory=tuple)


@dataclasses.dataclass(frozen=True, kw_only=True)
class Pathology(BELModelElement):
    namespace: str
    identifier: str


@dataclasses.dataclass(frozen=True, kw_only=True)
class PopulationAbundance(BELModelElement):
    namespace: str
    identifier: str
    location: Location | None = None


@dataclasses.dataclass(frozen=True, kw_only=True)
class ProteinModification(BELModelElement):
    namespace: str = ""
    identifier: str
    amino_acid: str | None = None
    residue: str | None = None


@dataclasses.dataclass(frozen=True, kw_only=True)
class ProteinAbundance(Abundance):
    namespace: str | None = None
    identifier: str | None = None
    fusion: Fusion | None = None
    variants: tuple[Variant] = dataclasses.field(default_factory=tuple)
    fragment: Fragment | None = None
    modifications: tuple[ProteinModification] = dataclasses.field(
        default_factory=tuple
    )


@dataclasses.dataclass(frozen=True, kw_only=True)
class Reaction(BELModelElement):
    reactants: frozenset[Abundance] = dataclasses.field(
        default_factory=frozenset
    )
    products: frozenset[Abundance] = dataclasses.field(
        default_factory=frozenset
    )


@dataclasses.dataclass(frozen=True, kw_only=True)
class RNAAbundance(Abundance):
    namespace: str | None
    identifier: str | None
    fusion: Fusion | None
    variants: tuple[Variant] = dataclasses.field(default_factory=tuple)


@dataclasses.dataclass(frozen=True, kw_only=True)
class Translocation(BELModelElement):
    abundance: Abundance
    from_namespace: str
    from_identifier: str
    to_namespace: str
    to_identifier: str


@dataclasses.dataclass(frozen=True, kw_only=True)
class Analogous(BELModelElement):
    source: BELModelElement
    target: BELModelElement


@dataclasses.dataclass(frozen=True, kw_only=True)
class Association(BELModelElement):
    source: BELModelElement
    target: BELModelElement


@dataclasses.dataclass(frozen=True, kw_only=True)
class BiomarkerFor(BELModelElement):
    source: BELModelElement
    target: BELModelElement


@dataclasses.dataclass(frozen=True, kw_only=True)
class CausesNoChange(BELModelElement):
    source: BELModelElement
    target: BELModelElement


@dataclasses.dataclass(frozen=True, kw_only=True)
class Decreases(BELModelElement):
    source: BELModelElement
    target: BELModelElement


@dataclasses.dataclass(frozen=True, kw_only=True)
class DirectlyDecreases(BELModelElement):
    source: BELModelElement
    target: BELModelElement


@dataclasses.dataclass(frozen=True, kw_only=True)
class DirectlyIncreases(BELModelElement):
    source: BELModelElement
    target: BELModelElement


@dataclasses.dataclass(frozen=True, kw_only=True)
class HasActivity(BELModelElement):
    source: BELModelElement
    target: BELModelElement


@dataclasses.dataclass(frozen=True, kw_only=True)
class HasComponent(BELModelElement):
    source: BELModelElement
    target: BELModelElement


@dataclasses.dataclass(frozen=True, kw_only=True)
class HasComponents(BELModelElement):
    source: BELModelElement
    target: BELModelElement


@dataclasses.dataclass(frozen=True, kw_only=True)
class HasMember(BELModelElement):
    source: BELModelElement
    target: BELModelElement


@dataclasses.dataclass(frozen=True, kw_only=True)
class HasMembers(BELModelElement):
    source: BELModelElement
    target: BELModelElement


@dataclasses.dataclass(frozen=True, kw_only=True)
class Increases(BELModelElement):
    source: BELModelElement
    target: BELModelElement


@dataclasses.dataclass(frozen=True, kw_only=True)
class IsA(BELModelElement):
    source: BELModelElement
    target: BELModelElement


@dataclasses.dataclass(frozen=True, kw_only=True)
class NegativeCorrelation(BELModelElement):
    source: BELModelElement
    target: BELModelElement


@dataclasses.dataclass(frozen=True, kw_only=True)
class Ortholgous(BELModelElement):
    source: BELModelElement
    target: BELModelElement


@dataclasses.dataclass(frozen=True, kw_only=True)
class PositiveCorrelation(BELModelElement):
    source: BELModelElement
    target: BELModelElement


@dataclasses.dataclass(frozen=True, kw_only=True)
class PrognosticBiomarkerFor(BELModelElement):
    source: BELModelElement
    target: BELModelElement


@dataclasses.dataclass(frozen=True, kw_only=True)
class RateLimitingStepFor(BELModelElement):
    source: BELModelElement
    target: BELModelElement


@dataclasses.dataclass(frozen=True, kw_only=True)
class Regulates(BELModelElement):
    source: BELModelElement
    target: BELModelElement


@dataclasses.dataclass(frozen=True, kw_only=True)
class SubProcessOf(BELModelElement):
    source: BELModelElement
    target: BELModelElement


@dataclasses.dataclass(frozen=True, kw_only=True)
class TranscribedTo(BELModelElement):
    source: BELModelElement
    target: BELModelElement


@dataclasses.dataclass(frozen=True, kw_only=True)
class TranslatedTo(BELModelElement):
    source: BELModelElement
    target: BELModelElement


@dataclasses.dataclass(frozen=True, kw_only=True)
class EquivalentTo(BELModelElement):
    source: BELModelElement
    target: BELModelElement


@dataclasses.dataclass(frozen=True, kw_only=True)
class NoCorrelation(BELModelElement):
    source: BELModelElement
    target: BELModelElement


@dataclasses.dataclass(frozen=True, kw_only=True)
class BELNamespaceDefinition(BELModelElement):
    name: str
    as_: str | tuple[str]


@dataclasses.dataclass(frozen=True, kw_only=True)
class BELGenericAnnotationDefinition:
    name: str
    as_: str | tuple[str]


@dataclasses.dataclass(frozen=True, kw_only=True)
class BELAnnotation:
    pass


@dataclasses.dataclass(frozen=True, kw_only=True)
class BELGenericAnnotation(BELAnnotation):
    """An annotation set on one or more statements.

    `name` is authoritative: four annotation names used by `SET` in real BEL
    documents (Citation, Support, STATEMENT_GROUP and DOCUMENT) have no
    `DEFINE ANNOTATION`, so `definition` is `None` for them. `definition` is a
    convenience back-pointer, populated when a definition exists; the complete
    set of definitions lives in `ReaderResult.annotation_definitions`.
    """

    name: str
    definition: BELGenericAnnotationDefinition | None = None
    args: tuple[str]


@dataclasses.dataclass(frozen=True, kw_only=True)
class BELDocumentAnnotation(BELAnnotation):
    name: str | None = None
    authors: str | None = None
    contact_info: str | None = None
    description: str | None = None
    licenses: str | None = None
    copyright: str | None = None
    version: str | None = None


@dataclasses.dataclass(frozen=True, kw_only=True)
class BELModel(momapy.core.Model):
    statements: frozenset[BELModelElement] = dataclasses.field(
        default_factory=frozenset
    )
    namespace_definitions: frozenset[BELNamespaceDefinition] = (
        dataclasses.field(default_factory=frozenset)
    )

    def is_submodel(self, other):
        return self.statements.issubset(
            other.statements
        ) and self.namespace_definitions.issubset(other.namespace_definitions)


@dataclasses.dataclass(frozen=True, kw_only=True)
class BELMap(momapy.core.map.Map):
    """Class for BEL maps.

    BEL support is model only: a `BELMap` carries a `BELModel` and no layout.
    The inherited `layout` and `layout_model_mapping` fields are always `None`.
    """

    model: BELModel | None = None


BELModelBuilder = momapy.builder.get_or_make_builder_cls(BELModel)
BELMapBuilder = momapy.builder.get_or_make_builder_cls(BELMap)


def _make_namespace_identifier_arg(namespace, identifier):
    if namespace:
        return f"{namespace}:{identifier}"
    else:
        return identifier


def _make_function_string(function_symbol, args):
    return f"{function_symbol}({', '.join(args)})"


def _make_relation_string(relation_symbol, source, target):
    return f"{source} {relation_symbol} {target}"
