import dataclasses

import momapy.core
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
    namespace: str | None = None
    identifier: str | None = None
    members: frozenset[Abundance] = dataclasses.field(
        default_factory=frozenset
    )


@dataclasses.dataclass(frozen=True, kw_only=True)
class CompositeAbundance(Abundance):
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
    variant: Variant | None


@dataclasses.dataclass(frozen=True, kw_only=True)
class MicroRNAAbundance(Abundance):
    namespace: str | None
    identifier: str | None
    fusion: Fusion | None
    variant: Variant | None


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
    variant: Variant | None = None
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
    variant: Variant | None


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
class BELNamespaceDefinition:
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

    def is_submodel(self, other):
        return self.statements.issubset(other.statements)


BELModelBuilder = momapy.builder.get_or_make_builder_cls(BELModel)


def _make_namespace_identifier_arg(namespace, identifier):
    if namespace:
        return f"{namespace}:{identifier}"
    else:
        return identifier


def _make_function_string(function_symbol, args):
    return f"{function_symbol}({', '.join(args)})"


def _make_relation_string(relation_symbol, source, target):
    return f"{source} {relation_symbol} {target}"
