"""BEL reader and writer.

Reads and writes BEL (Biological Expression Language) 2.1.3 documents.
`BELWriter` renders a `BELModel` back to BEL; `BELReader` parses a BEL
document into a `BELModel` plus the annotation contexts each statement was
asserted under.

`gmod()` is also read and written, as an extension: it is not part of BEL
2.1.3, but was introduced by PyBEL and is used by real documents.
"""

import collections
import dataclasses
import os
import typing

import frozendict
import pyparsing

import momapy.builder
import momapy.core
import momapy.core.map
import momapy.io.core
import momapy.utils
import momapy_bel.core


class BELWriter(momapy.io.core.Writer):
    _ELEMENT_CLS_TO_FUNC_NAME = {
        momapy_bel.core.List: "_list_to_string",
        momapy_bel.core.Location: "_location_to_string",
        momapy_bel.core.Abundance: "_abundance_to_string",
        momapy_bel.core.MolecularActivity: "_molecular_activity_to_string",
        momapy_bel.core.Activity: "_activity_to_string",
        momapy_bel.core.BiologicalProcess: "_biological_process_to_string",
        momapy_bel.core.CellSecretion: "_cell_secretion_to_string",
        momapy_bel.core.CellSurfaceExpression: "_cell_surface_expression_to_string",
        momapy_bel.core.ComplexAbundance: "_complex_abundance_to_string",
        momapy_bel.core.CompositeAbundance: "_composite_abundance_to_string",
        momapy_bel.core.Degradation: "_degradation_to_string",
        momapy_bel.core.Fragment: "_fragment_to_string",
        momapy_bel.core.Fusion: "_fusion_to_string",
        momapy_bel.core.Variant: "_variant_to_string",
        momapy_bel.core.GeneModification: "_gene_modification_to_string",
        momapy_bel.core.GeneAbundance: "_gene_abundance_to_string",
        momapy_bel.core.MicroRNAAbundance: "_microrna_abundance_to_string",
        momapy_bel.core.Pathology: "_pathology_to_string",
        momapy_bel.core.PopulationAbundance: "_population_abundance_to_string",
        momapy_bel.core.ProteinModification: "_protein_modification_to_string",
        momapy_bel.core.ProteinAbundance: "_protein_abundance_to_string",
        momapy_bel.core.Reaction: "_reaction_to_string",
        momapy_bel.core.RNAAbundance: "_rna_abundance_to_string",
        momapy_bel.core.Translocation: "_translocation_to_string",
        momapy_bel.core.Analogous: "_analogous_to_string",
        momapy_bel.core.Association: "_association_to_string",
        momapy_bel.core.BiomarkerFor: "_biomarker_for_to_string",
        momapy_bel.core.CausesNoChange: "_causes_no_change_to_string",
        momapy_bel.core.Decreases: "_decreases_to_string",
        momapy_bel.core.EquivalentTo: "_equivalent_to_to_string",
        momapy_bel.core.NoCorrelation: "_no_correlation_to_string",
        momapy_bel.core.DirectlyDecreases: "_directly_decreases_to_string",
        momapy_bel.core.DirectlyIncreases: "_directly_increases_to_string",
        momapy_bel.core.HasActivity: "_has_activity_to_string",
        momapy_bel.core.HasComponent: "_has_component_to_string",
        momapy_bel.core.HasComponents: "_has_components_to_string",
        momapy_bel.core.HasMember: "_has_member_to_string",
        momapy_bel.core.HasMembers: "_has_members_to_string",
        momapy_bel.core.Increases: "_increases_to_string",
        momapy_bel.core.IsA: "_is_a_to_string",
        momapy_bel.core.NegativeCorrelation: "_negative_correlation_to_string",
        momapy_bel.core.Ortholgous: "_orthologous_to_string",
        momapy_bel.core.PositiveCorrelation: "_positive_correlation_to_string",
        momapy_bel.core.PrognosticBiomarkerFor: "_prognostic_biomarker_for_to_string",
        momapy_bel.core.RateLimitingStepFor: "_rate_limiting_step_for_to_string",
        momapy_bel.core.Regulates: "_regulates_to_string",
        momapy_bel.core.SubProcessOf: "_subprocess_of_to_string",
        momapy_bel.core.TranscribedTo: "_transcribed_to_to_string",
        momapy_bel.core.TranslatedTo: "_translated_to_to_string",
        (
            momapy_bel.core.BELGenericAnnotation,
            "set",
        ): "_generic_annotation_to_set_string",
        (
            momapy_bel.core.BELGenericAnnotation,
            "unset",
        ): "_generic_annotation_to_unset_string",
        (
            momapy_bel.core.BELDocumentAnnotation,
            "set",
        ): "_document_annotation_to_set_string",
        momapy_bel.core.BELGenericAnnotationDefinition: "_annotation_definition_to_define_string",
        momapy_bel.core.BELNamespaceDefinition: "_namespace_definition_to_define_string",
    }

    @classmethod
    def write(
        cls,
        obj: "momapy_bel.core.BELModel | momapy_bel.core.BELMap",
        file_path: str | os.PathLike,
        element_to_annotations: dict | None = None,
        annotation_definitions: frozenset | None = None,
        with_abundances_as_statements: bool = False,
        with_biological_processes_as_statements: bool = False,
        with_reactions_as_statements: bool = True,
        with_degradations_as_statements: bool = True,
        with_translocations_as_statements: bool = True,
        **options: typing.Any,
    ) -> momapy.io.core.WriterResult:
        """Write a BEL model to a file.

        Args:
            obj: The `BELModel` to write, or a `BELMap` carrying one. The
                namespace definitions are taken from the model's
                `namespace_definitions` field.
            file_path: Path of the file to write to.
            element_to_annotations: Optional mapping from the model and its
                statements to the set of annotation contexts each was asserted
                under, as returned in `ReaderResult.element_to_annotations`.
                Each context is a `frozenset` of `BELAnnotation`s; a statement
                is emitted once per context.
            annotation_definitions: Optional set of
                `BELGenericAnnotationDefinition`s to emit as `DEFINE
                ANNOTATION` lines, as returned in
                `ReaderResult.annotation_definitions`.
            with_abundances_as_statements: Whether to emit bare abundances as
                statements. Defaults to `False`.
            with_biological_processes_as_statements: Whether to emit bare
                biological processes as statements. Defaults to `False`.
            with_reactions_as_statements: Whether to emit bare reactions as
                statements. Defaults to `True`.
            with_degradations_as_statements: Whether to emit bare degradations
                as statements. Defaults to `True`.
            with_translocations_as_statements: Whether to emit bare
                translocations as statements. Defaults to `True`.
            options: Additional writer options (accepted and ignored).

        Returns:
            A `WriterResult` holding `obj` and `file_path`.
        """
        if momapy.builder.isinstance_or_builder(obj, momapy.core.map.Map):
            bel_model = obj.model
        else:
            bel_model = obj
        if element_to_annotations is None:
            element_to_annotations = {}
        if annotation_definitions is None:
            annotation_definitions = frozenset()
        bel_string = cls._bel_model_to_string(
            bel_model,
            element_to_annotations,
            annotation_definitions,
            top_level_obj=obj,
            with_abundances_as_statements=with_abundances_as_statements,
            with_biological_processes_as_statements=with_biological_processes_as_statements,
            with_reactions_as_statements=with_reactions_as_statements,
            with_degradations_as_statements=with_degradations_as_statements,
            with_translocations_as_statements=with_translocations_as_statements,
        )
        with open(file_path, "w") as f:
            f.write(bel_string)
        return momapy.io.core.WriterResult(obj=obj, file_path=file_path)

    @classmethod
    def _list_to_string(cls, bel_list):
        args = [cls._bel_element_to_string(element) for element in bel_list.elements]
        return cls._make_function_string("list", args)

    @classmethod
    def _location_to_string(cls, location):
        args = [
            cls._make_namespace_identifier_arg(location.namespace, location.identifier)
        ]
        return cls._make_function_string("loc", args)

    @classmethod
    def _abundance_to_string(cls, abundance):
        args = [
            cls._make_namespace_identifier_arg(
                abundance.namespace, abundance.identifier
            )
        ]
        if abundance.location is not None:
            args.append(cls._bel_element_to_string(abundance.location))
        return cls._make_function_string("a", args)

    @classmethod
    def _molecular_activity_to_string(cls, molecular_activity):
        args = [
            cls._make_namespace_identifier_arg(
                molecular_activity.namespace, molecular_activity.identifier
            )
        ]
        return cls._make_function_string("ma", args)

    @classmethod
    def _activity_to_string(cls, activity):
        args = [cls._bel_element_to_string(activity.abundance)]
        if activity.molecular_activity is not None:
            args.append(cls._bel_element_to_string(activity.molecular_activity))
        return cls._make_function_string("act", args)

    @classmethod
    def _biological_process_to_string(cls, biological_process):
        args = [
            cls._make_namespace_identifier_arg(
                biological_process.namespace, biological_process.identifier
            )
        ]
        return cls._make_function_string("bp", args)

    @classmethod
    def _cell_secretion_to_string(cls, cell_secretion):
        args = [cls._bel_element_to_string(cell_secretion.abundance)]
        return cls._make_function_string("sec", args)

    @classmethod
    def _cell_surface_expression_to_string(cls, cell_surface_expression):
        args = [cls._bel_element_to_string(cell_surface_expression.abundance)]
        return cls._make_function_string("surf", args)

    @classmethod
    def _complex_abundance_to_string(cls, complex_abundance):
        if complex_abundance.members:
            args = sorted(
                [
                    cls._bel_element_to_string(member)
                    for member in complex_abundance.members
                ]
            )
        else:
            args = [
                cls._make_namespace_identifier_arg(
                    complex_abundance.namespace, complex_abundance.identifier
                )
            ]
        if complex_abundance.location is not None:
            args.append(cls._bel_element_to_string(complex_abundance.location))
        return cls._make_function_string("complex", args)

    @classmethod
    def _composite_abundance_to_string(cls, composite_abundance):
        if (
            composite_abundance.namespace is not None
            and composite_abundance.identifier is not None
        ):
            args = [
                cls._make_namespace_identifier_arg(
                    composite_abundance.namespace,
                    composite_abundance.identifier,
                )
            ]
        else:
            args = sorted(
                [
                    cls._bel_element_to_string(member)
                    for member in composite_abundance.members
                ]
            )
        return cls._make_function_string("composite", args)

    @classmethod
    def _degradation_to_string(cls, degradation):
        args = [cls._bel_element_to_string(degradation.abundance)]
        return cls._make_function_string("deg", args)

    @classmethod
    def _fragment_to_string(cls, fragment):
        # Always quoted: `?` is not a valid bareword and a start/stop such as
        # `672_713` must not be re-read as anything but one opaque token.
        args = [cls._quote(fragment.start_stop)]
        if fragment.descriptor is not None:
            args.append(cls._quote(fragment.descriptor))
        return cls._make_function_string("frag", args)

    @classmethod
    def _fusion_to_string(cls, fusion):
        args = [
            cls._make_namespace_identifier_arg(fusion.namespace5, fusion.identifier5),
            fusion.range5,
            cls._make_namespace_identifier_arg(fusion.namespace3, fusion.identifier3),
            fusion.range3,
        ]
        return cls._make_function_string("fus", args)

    @classmethod
    def _variant_to_string(cls, variant):
        # Always quoted: descriptors such as `P,S,396` contain commas and
        # would otherwise be re-read as three arguments.
        args = [cls._quote(variant.descriptor)]
        return cls._make_function_string("var", args)

    @classmethod
    def _gene_abundance_to_string(cls, gene_abundance):
        if (
            gene_abundance.namespace is not None
            and gene_abundance.identifier is not None
        ):
            args = [
                cls._make_namespace_identifier_arg(
                    gene_abundance.namespace, gene_abundance.identifier
                )
            ]
        else:
            args = [cls._bel_element_to_string(gene_abundance.fusion)]
        if gene_abundance.location is not None:
            args.append(cls._bel_element_to_string(gene_abundance.location))
        for variant in gene_abundance.variants:
            args.append(cls._bel_element_to_string(variant))
        for modification in gene_abundance.modifications:
            args.append(cls._bel_element_to_string(modification))
        return cls._make_function_string("g", args)

    @classmethod
    def _gene_modification_to_string(cls, gene_modification):
        args = [
            cls._make_namespace_identifier_arg(
                gene_modification.namespace, gene_modification.identifier
            )
        ]
        return cls._make_function_string("gmod", args)

    @classmethod
    def _microrna_abundance_to_string(cls, microrna):
        if microrna.namespace is not None and microrna.identifier is not None:
            args = [
                cls._make_namespace_identifier_arg(
                    microrna.namespace, microrna.identifier
                )
            ]
        else:
            args = [cls._bel_element_to_string(microrna.fusion)]
        if microrna.location is not None:
            args.append(cls._bel_element_to_string(microrna.location))
        for variant in microrna.variants:
            args.append(cls._bel_element_to_string(variant))
        for modification in microrna.modifications:
            args.append(cls._bel_element_to_string(modification))
        return cls._make_function_string("m", args)

    @classmethod
    def _pathology_to_string(cls, pathology):
        args = [
            cls._make_namespace_identifier_arg(
                pathology.namespace, pathology.identifier
            )
        ]
        return cls._make_function_string("path", args)

    @classmethod
    def _population_abundance_to_string(cls, population_abundance):
        args = [
            cls._make_namespace_identifier_arg(
                population_abundance.namespace, population_abundance.identifier
            )
        ]
        if population_abundance.location is not None:
            args.append(cls._bel_element_to_string(population_abundance.location))
        return cls._make_function_string("pop", args)

    @classmethod
    def _protein_modification_to_string(cls, protein_modification):
        args = [
            cls._make_namespace_identifier_arg(
                protein_modification.namespace, protein_modification.identifier
            )
        ]
        if protein_modification.amino_acid is not None:
            args.append(protein_modification.amino_acid)
        if protein_modification.residue is not None:
            args.append(protein_modification.residue)
        return cls._make_function_string("pmod", args)

    @classmethod
    def _protein_abundance_to_string(cls, protein_abundance):
        if (
            protein_abundance.namespace is not None
            and protein_abundance.identifier is not None
        ):
            args = [
                cls._make_namespace_identifier_arg(
                    protein_abundance.namespace, protein_abundance.identifier
                )
            ]
        else:
            args = [cls._bel_element_to_string(protein_abundance.fusion)]
        if protein_abundance.location is not None:
            args.append(cls._bel_element_to_string(protein_abundance.location))
        for variant in protein_abundance.variants:
            args.append(cls._bel_element_to_string(variant))
        if protein_abundance.fragment is not None:
            args.append(cls._bel_element_to_string(protein_abundance.fragment))
        if protein_abundance.modifications:
            for modification in protein_abundance.modifications:
                args.append(cls._bel_element_to_string(modification))
        return cls._make_function_string("p", args)

    @classmethod
    def _reaction_to_string(cls, reaction):
        args = [
            cls._make_function_string(
                "reactants",
                sorted(
                    [
                        cls._bel_element_to_string(reactant)
                        for reactant in reaction.reactants
                    ]
                ),
            ),
            cls._make_function_string(
                "products",
                sorted(
                    [
                        cls._bel_element_to_string(product)
                        for product in reaction.products
                    ]
                ),
            ),
        ]
        return cls._make_function_string("rxn", args)

    @classmethod
    def _rna_abundance_to_string(cls, rna_abundance):
        if rna_abundance.namespace is not None and rna_abundance.identifier is not None:
            args = [
                cls._make_namespace_identifier_arg(
                    rna_abundance.namespace, rna_abundance.identifier
                )
            ]
        else:
            args = [cls._bel_element_to_string(rna_abundance.fusion)]
        if rna_abundance.location is not None:
            args.append(cls._bel_element_to_string(rna_abundance.location))
        for variant in rna_abundance.variants:
            args.append(cls._bel_element_to_string(variant))
        for modification in rna_abundance.modifications:
            args.append(cls._bel_element_to_string(modification))
        return cls._make_function_string("r", args)

    @classmethod
    def _translocation_to_string(cls, translocation):
        args = [
            cls._bel_element_to_string(translocation.abundance),
            cls._make_function_string(
                "fromLoc",
                [
                    cls._make_namespace_identifier_arg(
                        translocation.from_namespace,
                        translocation.from_identifier,
                    )
                ],
            ),
            cls._make_function_string(
                "toLoc",
                [
                    cls._make_namespace_identifier_arg(
                        translocation.to_namespace, translocation.to_identifier
                    )
                ],
            ),
        ]
        return cls._make_function_string("tloc", args)

    @classmethod
    def _analogous_to_string(cls, analogous):
        return cls._make_relation_string(
            "analogous",
            cls._bel_element_to_string(analogous.source),
            cls._bel_element_to_string(analogous.target),
        )

    @classmethod
    def _association_to_string(cls, association):
        return cls._make_relation_string(
            "association",
            cls._bel_element_to_string(association.source),
            cls._bel_element_to_string(association.target),
        )

    @classmethod
    def _biomarker_for_to_string(cls, biomarker_for):
        return cls._make_relation_string(
            "biomarkerFor",
            cls._bel_element_to_string(biomarker_for.source),
            cls._bel_element_to_string(biomarker_for.target),
        )

    @classmethod
    def _causes_no_change_to_string(cls, causes_no_change):
        return cls._make_relation_string(
            "cnc",
            cls._bel_element_to_string(causes_no_change.source),
            cls._bel_element_to_string(causes_no_change.target),
        )

    @classmethod
    def _decreases_to_string(cls, decreases):
        return cls._make_relation_string(
            "-|",
            cls._bel_element_to_string(decreases.source),
            cls._bel_element_to_string(decreases.target),
        )

    @classmethod
    def _directly_decreases_to_string(cls, directly_decreases):
        return cls._make_relation_string(
            "=|",
            cls._bel_element_to_string(directly_decreases.source),
            cls._bel_element_to_string(directly_decreases.target),
        )

    @classmethod
    def _directly_increases_to_string(cls, directly_increases):
        return cls._make_relation_string(
            "=>",
            cls._bel_element_to_string(directly_increases.source),
            cls._bel_element_to_string(directly_increases.target),
        )

    @classmethod
    def _has_activity_to_string(cls, has_activity):
        return cls._make_relation_string(
            "hasActivity",
            cls._bel_element_to_string(has_activity.source),
            cls._bel_element_to_string(has_activity.target),
        )

    @classmethod
    def _has_component_to_string(cls, has_component):
        return cls._make_relation_string(
            "hasComponent",
            cls._bel_element_to_string(has_component.source),
            cls._bel_element_to_string(has_component.target),
        )

    @classmethod
    def _has_components_to_string(cls, has_components):
        return cls._make_relation_string(
            "hasComponents",
            cls._bel_element_to_string(has_components.source),
            cls._bel_element_to_string(has_components.target),
        )

    @classmethod
    def _has_member_to_string(cls, has_member):
        return cls._make_relation_string(
            "hasMember",
            cls._bel_element_to_string(has_member.source),
            cls._bel_element_to_string(has_member.target),
        )

    @classmethod
    def _has_members_to_string(cls, has_members):
        return cls._make_relation_string(
            "hasMembers",
            cls._bel_element_to_string(has_members.source),
            cls._bel_element_to_string(has_members.target),
        )

    @classmethod
    def _increases_to_string(cls, increases):
        return cls._make_relation_string(
            "->",
            cls._bel_element_to_string(increases.source),
            cls._bel_element_to_string(increases.target),
        )

    @classmethod
    def _is_a_to_string(cls, is_a):
        return cls._make_relation_string(
            "isA",
            cls._bel_element_to_string(is_a.source),
            cls._bel_element_to_string(is_a.target),
        )

    @classmethod
    def _negative_correlation_to_string(cls, negative_correlation):
        return cls._make_relation_string(
            "neg",
            cls._bel_element_to_string(negative_correlation.source),
            cls._bel_element_to_string(negative_correlation.target),
        )

    @classmethod
    def _orthologous_to_string(cls, orthologous):
        return cls._make_relation_string(
            "orthologous",
            cls._bel_element_to_string(orthologous.source),
            cls._bel_element_to_string(orthologous.target),
        )

    @classmethod
    def _positive_correlation_to_string(cls, positive_correlation):
        return cls._make_relation_string(
            "pos",
            cls._bel_element_to_string(positive_correlation.source),
            cls._bel_element_to_string(positive_correlation.target),
        )

    @classmethod
    def _prognostic_biomarker_for_to_string(cls, prognostic_biomarker_for):
        return cls._make_relation_string(
            "prognosticBiomarkerFor",
            cls._bel_element_to_string(prognostic_biomarker_for.source),
            cls._bel_element_to_string(prognostic_biomarker_for.target),
        )

    @classmethod
    def _rate_limiting_step_for_to_string(cls, rate_limiting_step_for):
        return cls._make_relation_string(
            "rateLimitingStepFor",
            cls._bel_element_to_string(rate_limiting_step_for.source),
            cls._bel_element_to_string(rate_limiting_step_for.target),
        )

    @classmethod
    def _regulates_to_string(cls, regulates):
        return cls._make_relation_string(
            "reg",
            cls._bel_element_to_string(regulates.source),
            cls._bel_element_to_string(regulates.target),
        )

    @classmethod
    def _subprocess_of_to_string(cls, subprocess_of):
        return cls._make_relation_string(
            "subProcessOf",
            cls._bel_element_to_string(subprocess_of.source),
            cls._bel_element_to_string(subprocess_of.target),
        )

    @classmethod
    def _transcribed_to_to_string(cls, transcribed_to):
        return cls._make_relation_string(
            ":>",
            cls._bel_element_to_string(transcribed_to.source),
            cls._bel_element_to_string(transcribed_to.target),
        )

    @classmethod
    def _translated_to_to_string(cls, translated_to):
        return cls._make_relation_string(
            ">>",
            cls._bel_element_to_string(translated_to.source),
            cls._bel_element_to_string(translated_to.target),
        )

    @classmethod
    def _equivalent_to_to_string(cls, equivalent_to):
        return cls._make_relation_string(
            "eq",
            cls._bel_element_to_string(equivalent_to.source),
            cls._bel_element_to_string(equivalent_to.target),
        )

    @classmethod
    def _no_correlation_to_string(cls, no_correlation):
        return cls._make_relation_string(
            "noCorrelation",
            cls._bel_element_to_string(no_correlation.source),
            cls._bel_element_to_string(no_correlation.target),
        )

    @classmethod
    def _generic_annotation_to_set_string(cls, annotation):
        # `annotation.name`, not `annotation.definition.name`: Citation,
        # Support, STATEMENT_GROUP and DOCUMENT are never `DEFINE`d, so their
        # back-pointer is `None`.
        return cls._make_set_string(annotation.name, annotation.args)

    @classmethod
    def _generic_annotation_to_unset_string(cls, annotation):
        return cls._make_unset_string(annotation.name)

    @classmethod
    def _annotation_definition_to_define_string(cls, annotation_definition):
        return cls._make_define_string(
            "ANNOTATION",
            annotation_definition.name,
            annotation_definition.as_type.value,
            annotation_definition.as_,
        )

    @classmethod
    def _namespace_definition_to_define_string(cls, namespace_definition):
        return cls._make_define_string(
            "NAMESPACE",
            namespace_definition.name,
            namespace_definition.as_type.value,
            namespace_definition.as_,
        )

    @classmethod
    def _document_annotation_to_set_string(cls, document_annotation):
        output_strings = []
        for attribute_name, document_key in (
            ("name", "Name"),
            ("description", "Description"),
            ("version", "Version"),
            ("authors", "Authors"),
            ("licenses", "Licenses"),
            ("copyright", "Copyright"),
            ("contact_info", "ContactInfo"),
        ):
            value = getattr(document_annotation, attribute_name)
            if value is not None:
                output_strings.append(
                    cls._make_set_string(f"DOCUMENT {document_key}", [value])
                )
        return "\n".join(output_strings)

    @classmethod
    def _annotation_sort_key(cls, bel_annotation):
        if isinstance(bel_annotation, momapy_bel.core.BELGenericAnnotation):
            return (bel_annotation.name, tuple(bel_annotation.args))
        return ("", ())

    @classmethod
    def _context_sort_key(cls, bel_context):
        return tuple(
            sorted(
                cls._annotation_sort_key(bel_annotation)
                for bel_annotation in bel_context
            )
        )

    @classmethod
    def _bel_model_to_string(
        cls,
        bel_model,
        bel_element_to_annotations,
        bel_annotation_definitions,
        top_level_obj=None,
        with_abundances_as_statements=False,
        with_biological_processes_as_statements=False,
        with_reactions_as_statements=True,
        with_degradations_as_statements=True,
        with_translocations_as_statements=True,
    ):
        output_strings = []
        # The document-level annotations may be keyed on the model or, when the
        # reader returned a map, on the map.
        bel_model_contexts = bel_element_to_annotations.get(bel_model)
        if bel_model_contexts is None and top_level_obj is not None:
            bel_model_contexts = bel_element_to_annotations.get(top_level_obj)
        if bel_model_contexts is not None:
            for bel_context in sorted(bel_model_contexts, key=cls._context_sort_key):
                for bel_annotation in sorted(
                    bel_context, key=cls._annotation_sort_key
                ):
                    output_strings.append(
                        cls._bel_annotation_to_string(bel_annotation)
                    )
        for bel_namespace_definition in sorted(
            bel_model.namespace_definitions, key=lambda definition: definition.name
        ):
            output_strings.append(
                cls._namespace_definition_to_define_string(bel_namespace_definition)
            )
        for bel_annotation_definition in sorted(
            bel_annotation_definitions, key=lambda definition: definition.name
        ):
            output_strings.append(
                cls._annotation_definition_to_define_string(bel_annotation_definition)
            )
        bel_statements = [
            bel_statement
            for bel_statement in bel_model.statements
            if (
                (
                    with_abundances_as_statements
                    or not isinstance(bel_statement, momapy_bel.core.Abundance)
                )
                and (
                    with_biological_processes_as_statements
                    or not isinstance(
                        bel_statement, momapy_bel.core.BiologicalProcess
                    )
                )
                and (
                    with_reactions_as_statements
                    or not isinstance(bel_statement, momapy_bel.core.Reaction)
                )
                and (
                    with_degradations_as_statements
                    or not isinstance(bel_statement, momapy_bel.core.Degradation)
                )
                and (
                    with_translocations_as_statements
                    or not isinstance(bel_statement, momapy_bel.core.Translocation)
                )
            )
        ]
        # `statements` is a frozenset, so emission order is not meaningful;
        # sorting makes the output deterministic.
        for statement_string, bel_statement in sorted(
            (
                (cls._bel_element_to_string(bel_statement), bel_statement)
                for bel_statement in bel_statements
            ),
            key=lambda pair: pair[0],
        ):
            bel_contexts = bel_element_to_annotations.get(bel_statement)
            if not bel_contexts:
                output_strings.append(statement_string)
                continue
            # One assertion of the statement per context it was made under.
            for bel_context in sorted(bel_contexts, key=cls._context_sort_key):
                unset_strings = []
                for bel_annotation in sorted(
                    bel_context, key=cls._annotation_sort_key
                ):
                    output_strings.append(
                        cls._bel_annotation_to_string(
                            bel_annotation, set_or_unset="set"
                        )
                    )
                    unset_strings.append(
                        cls._bel_annotation_to_string(
                            bel_annotation, set_or_unset="unset"
                        )
                    )
                output_strings.append(statement_string)
                output_strings.extend(unset_strings)
        bel_string = "\n".join(output_strings)
        return bel_string

    @classmethod
    def _quote(cls, value):
        # `\"` is the only escape BEL defines, so it is the only one to emit.
        escaped = value.replace('"', '\\"')
        return f'"{escaped}"'

    @classmethod
    def _make_namespace_identifier_arg(cls, namespace, identifier):
        if not identifier.isalnum():
            identifier = cls._quote(identifier)
        if namespace:
            return f"{namespace}:{identifier}"
        else:
            return identifier

    @classmethod
    def _make_function_string(cls, function_symbol, args):
        return f"{function_symbol}({', '.join(args)})"

    @classmethod
    def _make_relation_string(cls, relation_symbol, source, target):
        return f"{source} {relation_symbol} {target}"

    @classmethod
    def _make_set_or_define_args_string(cls, args):
        args = [cls._quote(arg) for arg in args]
        if len(args) > 1:
            args_string = f"{{{', '.join(args)}}}"
        else:
            args_string = str(args[0])
        return args_string

    @classmethod
    def _make_set_string(cls, annotation_name, args):
        args_string = cls._make_set_or_define_args_string(args)
        return f"SET {annotation_name} = {args_string}"

    @classmethod
    def _make_unset_string(cls, annotation_name):
        return f"UNSET {annotation_name}"

    @classmethod
    def _make_define_string(cls, define_type, definition_name, as_type, as_):
        args_string = cls._make_set_or_define_args_string(as_)
        return f"DEFINE {define_type} {definition_name} AS {as_type} {args_string}"

    @classmethod
    def _get_transformation_func(cls, key):
        transformation_func = cls._ELEMENT_CLS_TO_FUNC_NAME[key]
        return transformation_func

    @classmethod
    def _bel_element_to_string(cls, bel_element):
        key = type(bel_element)
        transformation_func_name = cls._get_transformation_func(key)
        transformation_func = getattr(cls, transformation_func_name)
        bel_string = transformation_func(bel_element)
        return bel_string

    @classmethod
    def _bel_annotation_to_string(cls, bel_annotation, set_or_unset="set"):
        key = (
            type(bel_annotation),
            set_or_unset,
        )
        transformation_func_name = cls._get_transformation_func(key)
        transformation_func = getattr(cls, transformation_func_name)
        bel_string = transformation_func(bel_annotation)
        return bel_string




@dataclasses.dataclass(kw_only=True)
class BELReaderResult(momapy.io.core.ReaderResult):
    """A `ReaderResult` carrying the BEL annotation vocabularies.

    BEL declares the vocabulary of its annotations in the document itself
    (`DEFINE ANNOTATION`), which no other momapy format does, so the
    definitions live on this subclass rather than on `ReaderResult`.

    This is the authoritative, complete set: unlike
    `BELGenericAnnotation.definition`, it also carries vocabularies that are
    declared but never `SET`.
    """

    annotation_definitions: frozenset | None = dataclasses.field(
        default=None,
        metadata={
            "description": (
                "The `BELGenericAnnotationDefinition`s declared by the "
                "document, including those never used by a `SET` line. `None` "
                "when the document was read with `with_annotations=False`."
            )
        },
    )


class BELReadingError(Exception):
    """Raised when a BEL document contains lines that cannot be parsed.

    Every failing line is reported: a third-party corpus is better served by
    one error listing all its bad lines than by dying on the first.
    """

    def __init__(self, file_path, failures):
        self.file_path = file_path
        self.failures = failures
        line_number, line, message = failures[0]
        super().__init__(
            f"{len(failures)} unparsable line(s) in '{file_path}'; "
            f"first at line {line_number}: {message} in {line!r}"
        )


class BELReader(momapy.io.core.Reader):
    """Reader for BEL 2.1.3 documents."""

    @classmethod
    def check_file(cls, file_path: str | os.PathLike) -> bool:
        """Return `True` if the file looks like a BEL document.

        The check is made on bytes rather than on decoded text, so it does not
        depend on the ambient locale encoding.

        Args:
            file_path: Path of the file to check.

        Returns:
            `True` if the file contains a `SET DOCUMENT` or a
            `DEFINE NAMESPACE` line, `False` otherwise.
        """
        try:
            with open(file_path, "rb") as f:
                for line in f:
                    if b"SET DOCUMENT" in line or b"DEFINE NAMESPACE" in line:
                        return True
            return False
        except Exception:
            return False

    @classmethod
    def read(
        cls,
        file_path: str | os.PathLike,
        return_type: typing.Literal["map", "model", "layout"] = "map",
        with_model: bool = True,
        with_layout: bool = True,
        with_annotations: bool = True,
        with_notes: bool = True,
        **options: typing.Any,
    ) -> momapy.io.core.ReaderResult:
        """Read a BEL document and return a reader result object.

        Args:
            file_path: Path of the BEL file to read.
            return_type: Shape of `result.obj`: `"map"` (default) returns a
                `BELMap`, `"model"` returns the bare `BELModel`. BEL has no
                layout, so `"layout"` raises `NotImplementedError`.
            with_model: Whether to build the model. When `False` and
                `return_type="map"`, the map's `model` is `None` and no
                statement is parsed. Defaults to `True`.
            with_layout: Accepted for signature parity with the other readers;
                BEL has no layout, so it has no effect. Defaults to `True`.
            with_annotations: Whether to read the `SET`/`UNSET` annotation
                contexts and the `DEFINE ANNOTATION` vocabularies. Defaults to
                `True`.
            with_notes: Accepted for signature parity; BEL has no notes, so it
                has no effect. Defaults to `True`.
            options: Additional reader-specific options (ignored).

        Returns:
            A `BELReaderResult` whose `obj` is a `BELMap` or a `BELModel`, whose
            `element_to_annotations` maps the object and each statement to the
            **set of annotation contexts** it was asserted under (each context
            a `frozenset` of `BELAnnotation`s), and whose
            `annotation_definitions` holds every `DEFINE ANNOTATION`
            vocabulary, including those never `SET`.

        Raises:
            NotImplementedError: If `return_type="layout"` (BEL has no layout).
            BELReadingError: If any line of the document fails to parse.
        """
        if return_type not in ("map", "model", "layout"):
            raise ValueError(
                f"invalid return_type {return_type!r}: expected 'map', "
                "'model' or 'layout'"
            )
        if return_type == "layout":
            raise NotImplementedError(
                "BEL has no layout; return_type='layout' is not supported"
            )
        momapy.utils.check_file_exists(file_path)
        with open(file_path, encoding="utf-8") as f:
            text = f.read()
        with_model = with_model or return_type == "model"
        if with_model:
            model_builder = momapy_bel.core.BELModelBuilder()
        else:
            model_builder = None
        failures = []
        # One `BELGenericAnnotation` per active annotation name: `SET`
        # replaces and `UNSET` deletes, so the one-value-per-name invariant
        # holds by construction and a context is never built by union.
        name_to_annotation = {}
        document_values = {}
        name_to_annotation_definition = {}
        element_to_contexts = collections.defaultdict(set)
        for line_number, line in _join_bel_continuation_lines(text):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            try:
                if _bel_set_keyword.matches(line, parse_all=False):
                    if not with_annotations:
                        continue
                    cls._handle_set_line(
                        line, name_to_annotation, document_values,
                        name_to_annotation_definition,
                    )
                elif _bel_unset_keyword.matches(line, parse_all=False):
                    if not with_annotations:
                        continue
                    cls._handle_unset_line(line, name_to_annotation)
                elif _bel_define_keyword.matches(line, parse_all=False):
                    cls._handle_define_line(
                        line,
                        model_builder,
                        name_to_annotation_definition,
                        with_annotations=with_annotations,
                    )
                elif model_builder is not None:
                    statement = parse_bel_statement(line)
                    model_builder.statements.add(statement)
                    if with_annotations:
                        element_to_contexts[statement].add(
                            frozenset(name_to_annotation.values())
                        )
            except pyparsing.ParseBaseException as exception:
                failures.append((line_number, line, str(exception)))
            except ValueError as exception:
                failures.append((line_number, line, str(exception)))
        if failures:
            raise BELReadingError(file_path, failures)
        if model_builder is not None:
            model = momapy.builder.object_from_builder(model_builder)
        else:
            model = None
        if return_type == "model":
            obj = model
        else:
            map_builder = momapy_bel.core.BELMapBuilder()
            map_builder.model = model
            obj = momapy.builder.object_from_builder(map_builder)
        if with_annotations:
            element_to_annotations = {
                element: frozenset(contexts)
                for element, contexts in element_to_contexts.items()
            }
            if document_values:
                element_to_annotations[obj] = frozenset(
                    [
                        frozenset(
                            [momapy_bel.core.BELDocumentAnnotation(**document_values)]
                        )
                    ]
                )
            element_to_annotations = frozendict.frozendict(element_to_annotations)
            annotation_definitions = frozenset(
                name_to_annotation_definition.values()
            )
        else:
            element_to_annotations = None
            annotation_definitions = None
        return BELReaderResult(
            obj=obj,
            element_to_annotations=element_to_annotations,
            element_to_notes=None,
            annotation_definitions=annotation_definitions,
            id_to_element=_build_bel_id_to_element(obj, model),
            source_id_to_model_element=None,
            source_id_to_layout_element=None,
            source_id_to_annotations=None,
            source_id_to_notes=None,
            file_path=file_path,
        )

    @classmethod
    def _handle_set_line(
        cls,
        line: str,
        name_to_annotation: dict,
        document_values: dict,
        name_to_annotation_definition: dict,
    ) -> None:
        results = _bel_set_line.parse_string(line, parse_all=True)
        if results[0] == "DOCUMENT":
            document_key = results[1]
            attribute_name = _BEL_DOCUMENT_KEY_TO_ATTRIBUTE_NAME.get(document_key)
            if attribute_name is None:
                raise ValueError(f"unknown SET DOCUMENT key '{document_key}'")
            document_values[attribute_name] = results[2]
            return
        name = results[0]
        args = tuple(results[1])
        annotation = momapy_bel.core.BELGenericAnnotation(
            name=name,
            definition=name_to_annotation_definition.get(name),
            args=args,
        )
        name_to_annotation[_bel_annotation_slot(name)] = annotation

    @classmethod
    def _handle_unset_line(cls, line: str, name_to_annotation: dict) -> None:
        results = _bel_unset_line.parse_string(line, parse_all=True)
        if results[0] == "ALL":
            name_to_annotation.clear()
            return
        for name in results:
            # Unsetting an inactive name is not an error: real documents unset
            # names they never set and leave statement groups unclosed.
            name_to_annotation.pop(_bel_annotation_slot(name), None)

    @classmethod
    def _handle_define_line(
        cls,
        line: str,
        model_builder: typing.Any,
        name_to_annotation_definition: dict,
        with_annotations: bool = True,
    ) -> None:
        results = _bel_define_line.parse_string(line, parse_all=True)
        define_type, name = results[0], results[1]
        as_type = momapy_bel.core.BELAsDefinitionType(results[2])
        as_ = tuple(results[3])
        if (
            as_type is not momapy_bel.core.BELAsDefinitionType.LIST
            and len(as_) != 1
        ):
            raise ValueError(
                f"DEFINE {define_type} {name} AS {as_type.value} takes one "
                f"value, got {len(as_)}"
            )
        if define_type == "NAMESPACE":
            if model_builder is not None:
                model_builder.namespace_definitions.add(
                    momapy_bel.core.BELNamespaceDefinition(
                        name=name, as_type=as_type, as_=as_
                    )
                )
        elif with_annotations:
            name_to_annotation_definition[name] = (
                momapy_bel.core.BELGenericAnnotationDefinition(
                    name=name, as_type=as_type, as_=as_
                )
            )


def parse_bel_statement(line: str) -> momapy_bel.core.BELModelElement:
    """Parse one BEL statement line into a model element.

    A line with no relation (a term-only line) yields the term itself; a line
    with a relation yields the corresponding relation element.

    Args:
        line: One logical BEL line, with continuations already joined.

    Returns:
        The `BELModelElement` the line denotes.

    Raises:
        pyparsing.ParseBaseException: If the line is not a BEL statement.
        ValueError: If a function is unknown or is given an argument it cannot
            hold.
    """
    return _bel_statement.parse_string(line, parse_all=True)[0]


def count_bel_call_sites(element: momapy_bel.core.BELModelElement) -> int:
    """Count the BEL function call sites a parsed element accounts for.

    Every element maps to one written function call, except `Reaction` and
    `Translocation`, which absorb the `reactants`/`products` and
    `fromLoc`/`toLoc` calls of their source text into their own fields.

    Args:
        element: A parsed BEL element.

    Returns:
        The number of function call sites in the text `element` was parsed
        from, counted with multiplicity.
    """
    if not isinstance(element, momapy_bel.core.BELModelElement):
        return 0
    if isinstance(element, tuple(_BEL_RELATION_TO_CLASS.values())):
        # A relation is written as an infix operator, not as a function call.
        count = 0
    elif isinstance(
        element, (momapy_bel.core.Reaction, momapy_bel.core.Translocation)
    ):
        # `rxn()` absorbs `reactants()`/`products()`, `tloc()` absorbs
        # `fromLoc()`/`toLoc()`.
        count = 3
    else:
        count = 1
    for field in dataclasses.fields(element):
        value = getattr(element, field.name)
        if isinstance(value, (frozenset, tuple)):
            for item in value:
                count += count_bel_call_sites(item)
        else:
            count += count_bel_call_sites(value)
    return count


_BEL_WORD_RELATION_TO_CLASS = {
    "analogous": momapy_bel.core.Analogous,
    "association": momapy_bel.core.Association,
    "biomarkerFor": momapy_bel.core.BiomarkerFor,
    "causesNoChange": momapy_bel.core.CausesNoChange,
    "cnc": momapy_bel.core.CausesNoChange,
    "decreases": momapy_bel.core.Decreases,
    "directlyDecreases": momapy_bel.core.DirectlyDecreases,
    "directlyIncreases": momapy_bel.core.DirectlyIncreases,
    "equivalentTo": momapy_bel.core.EquivalentTo,
    "eq": momapy_bel.core.EquivalentTo,
    "hasActivity": momapy_bel.core.HasActivity,
    "hasComponent": momapy_bel.core.HasComponent,
    "hasComponents": momapy_bel.core.HasComponents,
    "hasMember": momapy_bel.core.HasMember,
    "hasMembers": momapy_bel.core.HasMembers,
    "increases": momapy_bel.core.Increases,
    "isA": momapy_bel.core.IsA,
    "negativeCorrelation": momapy_bel.core.NegativeCorrelation,
    "neg": momapy_bel.core.NegativeCorrelation,
    "noCorrelation": momapy_bel.core.NoCorrelation,
    "orthologous": momapy_bel.core.Ortholgous,
    "positiveCorrelation": momapy_bel.core.PositiveCorrelation,
    "pos": momapy_bel.core.PositiveCorrelation,
    "prognosticBiomarkerFor": momapy_bel.core.PrognosticBiomarkerFor,
    "rateLimitingStepOf": momapy_bel.core.RateLimitingStepFor,
    # Not a BEL 2.1.3 spelling, but the one `BELWriter` emits, so that the
    # writer's own output reparses.
    "rateLimitingStepFor": momapy_bel.core.RateLimitingStepFor,
    "regulates": momapy_bel.core.Regulates,
    "reg": momapy_bel.core.Regulates,
    "subProcessOf": momapy_bel.core.SubProcessOf,
    "transcribedTo": momapy_bel.core.TranscribedTo,
    "translatedTo": momapy_bel.core.TranslatedTo,
}
_BEL_SYMBOL_RELATION_TO_CLASS = {
    "->": momapy_bel.core.Increases,
    "-|": momapy_bel.core.Decreases,
    "=>": momapy_bel.core.DirectlyIncreases,
    "=|": momapy_bel.core.DirectlyDecreases,
    "--": momapy_bel.core.Association,
    ":>": momapy_bel.core.TranscribedTo,
    ">>": momapy_bel.core.TranslatedTo,
}
_BEL_RELATION_TO_CLASS = (
    _BEL_WORD_RELATION_TO_CLASS | _BEL_SYMBOL_RELATION_TO_CLASS
)
_BEL_DOCUMENT_KEY_TO_ATTRIBUTE_NAME = {
    "Name": "name",
    "Description": "description",
    "Version": "version",
    "Authors": "authors",
    "Licenses": "licenses",
    "Copyright": "copyright",
    "ContactInfo": "contact_info",
}
# `Support` and `Evidence` are the same slot: documents `SET Support` and
# `UNSET Evidence` interchangeably.
_BEL_ANNOTATION_NAME_TO_SLOT = {"Evidence": "Support"}


_bel_quoted_string = pyparsing.quoted_string.copy()
_bel_bareword = pyparsing.Word(pyparsing.alphanums + "_.")
_bel_namespace = pyparsing.Word(pyparsing.alphas, pyparsing.alphanums + "_")
_bel_namespace_identifier = (
    _bel_namespace
    + pyparsing.Suppress(":")
    + (_bel_quoted_string | _bel_bareword)
)
_bel_function_name = pyparsing.Word(pyparsing.alphas, pyparsing.alphanums + "_")
# Two `one_of` calls: `as_keyword=True` wraps each token in `\b` and so cannot
# match `->`, which has no word character at either end.
_bel_relation = pyparsing.one_of(
    list(_BEL_WORD_RELATION_TO_CLASS), as_keyword=True
) | pyparsing.one_of(list(_BEL_SYMBOL_RELATION_TO_CLASS))
_bel_annotation_name = pyparsing.Word(
    pyparsing.alphas + "_", pyparsing.alphanums + "_"
)
_bel_annotation_value = _bel_quoted_string | _bel_bareword
_bel_annotation_value_list = (
    pyparsing.Suppress("{")
    + pyparsing.DelimitedList(_bel_annotation_value, ",")
    + pyparsing.Suppress("}")
)
_bel_set_keyword = pyparsing.Keyword("SET")
_bel_unset_keyword = pyparsing.Keyword("UNSET")
_bel_define_keyword = pyparsing.Keyword("DEFINE")

_bel_term = pyparsing.Forward()
_bel_statement = pyparsing.Forward()
_bel_argument = (
    _bel_term
    | pyparsing.Group(_bel_namespace_identifier)
    | _bel_quoted_string
    | _bel_bareword
)
_bel_term <<= (
    _bel_function_name
    + pyparsing.Suppress("(")
    + pyparsing.DelimitedList(_bel_argument, ",")
    + pyparsing.Suppress(")")
)
_bel_object = _bel_term | (
    pyparsing.Suppress("(") + _bel_statement + pyparsing.Suppress(")")
)
_bel_statement <<= _bel_term + pyparsing.Opt(_bel_relation + _bel_object)
_bel_document_set_line = (
    pyparsing.Suppress(_bel_set_keyword)
    + pyparsing.Keyword("DOCUMENT")
    + _bel_annotation_name
    + pyparsing.Suppress("=")
    + _bel_quoted_string
)
_bel_generic_set_line = (
    pyparsing.Suppress(_bel_set_keyword)
    + _bel_annotation_name
    + pyparsing.Suppress("=")
    + pyparsing.Group(_bel_annotation_value_list | _bel_annotation_value)
)
_bel_set_line = _bel_document_set_line | _bel_generic_set_line
_bel_unset_line = pyparsing.Suppress(_bel_unset_keyword) + (
    pyparsing.Keyword("ALL")
    | (
        pyparsing.Suppress("{")
        + pyparsing.DelimitedList(_bel_annotation_name, ",")
        + pyparsing.Suppress("}")
    )
    | _bel_annotation_name
)
_bel_define_line = (
    pyparsing.Suppress(_bel_define_keyword)
    + (pyparsing.Keyword("NAMESPACE") | pyparsing.Keyword("ANNOTATION"))
    + _bel_annotation_name
    + pyparsing.Suppress(pyparsing.Keyword("AS"))
    + (
        pyparsing.Keyword("URL")
        | pyparsing.Keyword("PATTERN")
        | pyparsing.Keyword("LIST")
    )
    + pyparsing.Group(_bel_annotation_value_list | _bel_annotation_value)
)


@_bel_quoted_string.set_parse_action
def _resolve_bel_quoted_string(results: pyparsing.ParseResults) -> str:
    # `\"` is the only real escape in BEL; every other backslash is literal,
    # which rules out both `QuotedString(esc_char="\\")` and a regex that
    # treats `\d` as an escape.
    return results[0][1:-1].replace('\\"', '"')


@_bel_term.set_parse_action
def _resolve_bel_term(
    results: pyparsing.ParseResults,
) -> momapy_bel.core.BELModelElement:
    # The single place the one-element-one-parse-action convention does not
    # apply: `p`/`path`/`pop`/`pmod`/`products` and friends prefix-collide, so
    # the function name is resolved by dict lookup rather than by alternation.
    function_name = results[0]
    resolver = _BEL_FUNCTION_NAME_TO_RESOLVER.get(function_name)
    if resolver is None:
        raise ValueError(f"unknown BEL function '{function_name}'")
    return resolver(list(results)[1:])


@_bel_statement.set_parse_action
def _resolve_bel_statement(
    results: pyparsing.ParseResults,
) -> momapy_bel.core.BELModelElement:
    if len(results) == 1:
        return results[0]
    source, relation, target = results
    return _BEL_RELATION_TO_CLASS[relation](source=source, target=target)


@dataclasses.dataclass(frozen=True)
class _BELReactants:
    """The `reactants()` modifier of a `rxn()` term, absorbed by `Reaction`."""

    members: tuple


@dataclasses.dataclass(frozen=True)
class _BELProducts:
    """The `products()` modifier of a `rxn()` term, absorbed by `Reaction`."""

    members: tuple


@dataclasses.dataclass(frozen=True)
class _BELFromLocation:
    """The `fromLoc()` modifier of a `tloc()` term."""

    namespace: str
    identifier: str


@dataclasses.dataclass(frozen=True)
class _BELToLocation:
    """The `toLoc()` modifier of a `tloc()` term."""

    namespace: str
    identifier: str


def _join_bel_continuation_lines(text: str) -> list[tuple[int, str]]:
    logical_lines = []
    buffer = None
    start_line_number = 0
    for line_number, line in enumerate(text.split("\n"), start=1):
        if buffer is None:
            buffer = ""
            start_line_number = line_number
        if line.endswith("\\"):
            buffer += line[:-1]
            continue
        logical_lines.append((start_line_number, buffer + line))
        buffer = None
    if buffer is not None:
        logical_lines.append((start_line_number, buffer))
    return logical_lines


def _build_bel_id_to_element(
    obj: typing.Any, model: "momapy_bel.core.BELModel | None"
) -> frozendict.frozendict:
    # Deliberately not `momapy.io._utils.build_id_to_element`: that module is
    # private to momapy, and BEL has no layout, so the walk is three lines.
    id_to_element = {}
    if model is not None:
        id_to_element[model.id_] = model
        for element in model.descendants():
            id_to_element[element.id_] = element
    id_to_element[obj.id_] = obj
    return frozendict.frozendict(id_to_element)


def _bel_annotation_slot(name: str) -> str:
    return _BEL_ANNOTATION_NAME_TO_SLOT.get(name, name)


def _bel_namespace_identifier_from_argument(
    argument: typing.Any, function_name: str
) -> tuple[str, str]:
    if isinstance(argument, pyparsing.ParseResults) and len(argument) == 2:
        return argument[0], argument[1]
    raise ValueError(
        f"expected a namespace:identifier argument in {function_name}(), "
        f"got {argument!r}"
    )


def _resolve_bel_location(arguments: list) -> momapy_bel.core.Location:
    if len(arguments) != 1:
        raise ValueError(f"loc() takes one argument, got {len(arguments)}")
    namespace, identifier = _bel_namespace_identifier_from_argument(
        arguments[0], "loc"
    )
    return momapy_bel.core.Location(namespace=namespace, identifier=identifier)


def _resolve_bel_from_location(arguments: list) -> _BELFromLocation:
    if len(arguments) != 1:
        raise ValueError(f"fromLoc() takes one argument, got {len(arguments)}")
    namespace, identifier = _bel_namespace_identifier_from_argument(
        arguments[0], "fromLoc"
    )
    return _BELFromLocation(namespace=namespace, identifier=identifier)


def _resolve_bel_to_location(arguments: list) -> _BELToLocation:
    if len(arguments) != 1:
        raise ValueError(f"toLoc() takes one argument, got {len(arguments)}")
    namespace, identifier = _bel_namespace_identifier_from_argument(
        arguments[0], "toLoc"
    )
    return _BELToLocation(namespace=namespace, identifier=identifier)


def _resolve_bel_molecular_activity(
    arguments: list,
) -> momapy_bel.core.MolecularActivity:
    if len(arguments) != 1:
        raise ValueError(f"ma() takes one argument, got {len(arguments)}")
    argument = arguments[0]
    if isinstance(argument, str):
        # A default-namespace activity such as `ma(kin)`.
        return momapy_bel.core.MolecularActivity(namespace="", identifier=argument)
    namespace, identifier = _bel_namespace_identifier_from_argument(argument, "ma")
    return momapy_bel.core.MolecularActivity(
        namespace=namespace, identifier=identifier
    )


def _resolve_bel_gene_modification(
    arguments: list,
) -> momapy_bel.core.GeneModification:
    if len(arguments) != 1:
        raise ValueError(f"gmod() takes one argument, got {len(arguments)}")
    argument = arguments[0]
    if isinstance(argument, str):
        # A default-namespace modification such as `gmod(Me)`.
        return momapy_bel.core.GeneModification(
            namespace="", identifier=argument
        )
    namespace, identifier = _bel_namespace_identifier_from_argument(
        argument, "gmod"
    )
    return momapy_bel.core.GeneModification(
        namespace=namespace, identifier=identifier
    )


def _resolve_bel_variant(arguments: list) -> momapy_bel.core.Variant:
    if len(arguments) != 1:
        raise ValueError(f"var() takes one argument, got {len(arguments)}")
    if not isinstance(arguments[0], str):
        raise ValueError(f"unexpected argument {arguments[0]!r} in var()")
    # Stored verbatim: real documents violate the spec's `/\S+/` with values
    # such as `G, 93, A`, and BEL1-to-HGVS normalisation is a separate pass.
    return momapy_bel.core.Variant(descriptor=arguments[0])


def _resolve_bel_fragment(arguments: list) -> momapy_bel.core.Fragment:
    if not 1 <= len(arguments) <= 2:
        raise ValueError(
            f"frag() takes one or two arguments, got {len(arguments)}"
        )
    for argument in arguments:
        if not isinstance(argument, str):
            raise ValueError(f"unexpected argument {argument!r} in frag()")
    descriptor = arguments[1] if len(arguments) == 2 else None
    return momapy_bel.core.Fragment(
        start_stop=arguments[0], descriptor=descriptor
    )


def _resolve_bel_fusion(arguments: list) -> momapy_bel.core.Fusion:
    if len(arguments) not in (2, 4):
        raise ValueError(
            f"fus() takes two or four arguments, got {len(arguments)}"
        )
    namespace5, identifier5 = _bel_namespace_identifier_from_argument(
        arguments[0], "fus"
    )
    if len(arguments) == 2:
        range5 = "?"
        three_prime_argument = arguments[1]
        range3 = "?"
    else:
        range5 = arguments[1]
        three_prime_argument = arguments[2]
        range3 = arguments[3]
    namespace3, identifier3 = _bel_namespace_identifier_from_argument(
        three_prime_argument, "fus"
    )
    return momapy_bel.core.Fusion(
        namespace5=namespace5,
        identifier5=identifier5,
        range5=range5,
        namespace3=namespace3,
        identifier3=identifier3,
        range3=range3,
    )


def _resolve_bel_protein_modification(
    arguments: list,
) -> momapy_bel.core.ProteinModification:
    if not 1 <= len(arguments) <= 3:
        raise ValueError(
            f"pmod() takes one to three arguments, got {len(arguments)}"
        )
    argument = arguments[0]
    if isinstance(argument, str):
        namespace, identifier = "", argument
    else:
        namespace, identifier = _bel_namespace_identifier_from_argument(
            argument, "pmod"
        )
    # The remaining arguments are positional, not looked up: eight one-letter
    # codes are both BEL1 modification types and amino-acid abbreviations, so
    # `pmod(Ph, P, 213)` only resolves by position.
    amino_acid = None
    residue = None
    if len(arguments) > 1:
        if not isinstance(arguments[1], str):
            raise ValueError(f"unexpected argument {arguments[1]!r} in pmod()")
        amino_acid = arguments[1]
    if len(arguments) > 2:
        if not isinstance(arguments[2], str):
            raise ValueError(f"unexpected argument {arguments[2]!r} in pmod()")
        residue = arguments[2]
    return momapy_bel.core.ProteinModification(
        namespace=namespace,
        identifier=identifier,
        amino_acid=amino_acid,
        residue=residue,
    )


def _resolve_bel_abundance(arguments: list) -> momapy_bel.core.Abundance:
    namespace, identifier = _bel_namespace_identifier_from_argument(
        arguments[0], "a"
    )
    location = None
    for argument in arguments[1:]:
        if isinstance(argument, momapy_bel.core.Location):
            location = argument
        else:
            raise ValueError(f"unexpected argument {argument!r} in a()")
    return momapy_bel.core.Abundance(
        namespace=namespace, identifier=identifier, location=location
    )


def _resolve_bel_population_abundance(
    arguments: list,
) -> momapy_bel.core.PopulationAbundance:
    namespace, identifier = _bel_namespace_identifier_from_argument(
        arguments[0], "pop"
    )
    location = None
    for argument in arguments[1:]:
        if isinstance(argument, momapy_bel.core.Location):
            location = argument
        else:
            raise ValueError(f"unexpected argument {argument!r} in pop()")
    return momapy_bel.core.PopulationAbundance(
        namespace=namespace, identifier=identifier, location=location
    )


def _resolve_bel_biological_process(
    arguments: list,
) -> momapy_bel.core.BiologicalProcess:
    if len(arguments) != 1:
        raise ValueError(f"bp() takes one argument, got {len(arguments)}")
    namespace, identifier = _bel_namespace_identifier_from_argument(
        arguments[0], "bp"
    )
    return momapy_bel.core.BiologicalProcess(
        namespace=namespace, identifier=identifier
    )


def _resolve_bel_pathology(arguments: list) -> momapy_bel.core.Pathology:
    if len(arguments) != 1:
        raise ValueError(f"path() takes one argument, got {len(arguments)}")
    namespace, identifier = _bel_namespace_identifier_from_argument(
        arguments[0], "path"
    )
    return momapy_bel.core.Pathology(namespace=namespace, identifier=identifier)


def _resolve_bel_activity(arguments: list) -> momapy_bel.core.Activity:
    if not arguments:
        raise ValueError("act() takes at least one argument, got none")
    abundance = arguments[0]
    if not isinstance(abundance, momapy_bel.core.BELModelElement):
        raise ValueError(f"unexpected argument {abundance!r} in act()")
    molecular_activity = None
    for argument in arguments[1:]:
        if isinstance(argument, momapy_bel.core.MolecularActivity):
            molecular_activity = argument
        else:
            raise ValueError(f"unexpected argument {argument!r} in act()")
    return momapy_bel.core.Activity(
        abundance=abundance, molecular_activity=molecular_activity
    )


def _make_bel_single_abundance_resolver(
    cls: type, function_name: str
) -> typing.Callable[[list], momapy_bel.core.BELModelElement]:
    def _resolve(arguments: list) -> momapy_bel.core.BELModelElement:
        if len(arguments) != 1:
            raise ValueError(
                f"{function_name}() takes one argument, got {len(arguments)}"
            )
        abundance = arguments[0]
        if not isinstance(abundance, momapy_bel.core.BELModelElement):
            raise ValueError(
                f"unexpected argument {abundance!r} in {function_name}()"
            )
        return cls(abundance=abundance)

    return _resolve


_resolve_bel_cell_secretion = _make_bel_single_abundance_resolver(
    momapy_bel.core.CellSecretion, "sec"
)
_resolve_bel_cell_surface_expression = _make_bel_single_abundance_resolver(
    momapy_bel.core.CellSurfaceExpression, "surf"
)
_resolve_bel_degradation = _make_bel_single_abundance_resolver(
    momapy_bel.core.Degradation, "deg"
)


def _make_bel_nucleic_abundance_resolver(
    cls: type, function_name: str
) -> typing.Callable[[list], momapy_bel.core.BELModelElement]:
    def _resolve(arguments: list) -> momapy_bel.core.BELModelElement:
        if not arguments:
            raise ValueError(
                f"{function_name}() takes at least one argument, got none"
            )
        namespace = None
        identifier = None
        fusion = None
        if isinstance(arguments[0], momapy_bel.core.Fusion):
            fusion = arguments[0]
        else:
            namespace, identifier = _bel_namespace_identifier_from_argument(
                arguments[0], function_name
            )
        location = None
        variants = []
        modifications = []
        for argument in arguments[1:]:
            if isinstance(argument, momapy_bel.core.Location):
                location = argument
            elif isinstance(argument, momapy_bel.core.Variant):
                variants.append(argument)
            elif isinstance(argument, momapy_bel.core.GeneModification):
                modifications.append(argument)
            else:
                raise ValueError(
                    f"unexpected argument {argument!r} in {function_name}()"
                )
        return cls(
            namespace=namespace,
            identifier=identifier,
            fusion=fusion,
            location=location,
            variants=tuple(variants),
            modifications=tuple(modifications),
        )

    return _resolve


_resolve_bel_gene_abundance = _make_bel_nucleic_abundance_resolver(
    momapy_bel.core.GeneAbundance, "g"
)
_resolve_bel_rna_abundance = _make_bel_nucleic_abundance_resolver(
    momapy_bel.core.RNAAbundance, "r"
)
_resolve_bel_microrna_abundance = _make_bel_nucleic_abundance_resolver(
    momapy_bel.core.MicroRNAAbundance, "m"
)


def _resolve_bel_protein_abundance(
    arguments: list,
) -> momapy_bel.core.ProteinAbundance:
    if not arguments:
        raise ValueError("p() takes at least one argument, got none")
    namespace = None
    identifier = None
    fusion = None
    if isinstance(arguments[0], momapy_bel.core.Fusion):
        fusion = arguments[0]
    else:
        namespace, identifier = _bel_namespace_identifier_from_argument(
            arguments[0], "p"
        )
    location = None
    fragment = None
    variants = []
    modifications = []
    # A loop with an `else` that raises, not a `next(...)` over a predicate: a
    # protein carrying two `var()`s must keep both, and an unknown modifier
    # must be reported rather than silently dropped.
    for argument in arguments[1:]:
        if isinstance(argument, momapy_bel.core.Location):
            location = argument
        elif isinstance(argument, momapy_bel.core.Fragment):
            fragment = argument
        elif isinstance(argument, momapy_bel.core.Variant):
            variants.append(argument)
        elif isinstance(argument, momapy_bel.core.ProteinModification):
            modifications.append(argument)
        else:
            raise ValueError(f"unexpected argument {argument!r} in p()")
    return momapy_bel.core.ProteinAbundance(
        namespace=namespace,
        identifier=identifier,
        location=location,
        fusion=fusion,
        variants=tuple(variants),
        fragment=fragment,
        modifications=tuple(modifications),
    )


def _make_bel_members_abundance_resolver(
    cls: type, function_name: str
) -> typing.Callable[[list], momapy_bel.core.BELModelElement]:
    def _resolve(arguments: list) -> momapy_bel.core.BELModelElement:
        if not arguments:
            raise ValueError(
                f"{function_name}() takes at least one argument, got none"
            )
        namespace = None
        identifier = None
        members = []
        location = None
        rest = arguments
        if not isinstance(arguments[0], momapy_bel.core.BELModelElement):
            namespace, identifier = _bel_namespace_identifier_from_argument(
                arguments[0], function_name
            )
            rest = arguments[1:]
        for argument in rest:
            if isinstance(argument, momapy_bel.core.Location):
                location = argument
            elif isinstance(argument, momapy_bel.core.BELModelElement):
                members.append(argument)
            else:
                raise ValueError(
                    f"unexpected argument {argument!r} in {function_name}()"
                )
        return cls(
            namespace=namespace,
            identifier=identifier,
            members=frozenset(members),
            location=location,
        )

    return _resolve


_resolve_bel_complex_abundance = _make_bel_members_abundance_resolver(
    momapy_bel.core.ComplexAbundance, "complex"
)
_resolve_bel_composite_abundance = _make_bel_members_abundance_resolver(
    momapy_bel.core.CompositeAbundance, "composite"
)


def _resolve_bel_list(arguments: list) -> momapy_bel.core.List:
    for argument in arguments:
        if not isinstance(argument, momapy_bel.core.BELModelElement):
            raise ValueError(f"unexpected argument {argument!r} in list()")
    return momapy_bel.core.List(elements=tuple(arguments))


def _resolve_bel_reactants(arguments: list) -> _BELReactants:
    for argument in arguments:
        if not isinstance(argument, momapy_bel.core.BELModelElement):
            raise ValueError(f"unexpected argument {argument!r} in reactants()")
    return _BELReactants(members=tuple(arguments))


def _resolve_bel_products(arguments: list) -> _BELProducts:
    for argument in arguments:
        if not isinstance(argument, momapy_bel.core.BELModelElement):
            raise ValueError(f"unexpected argument {argument!r} in products()")
    return _BELProducts(members=tuple(arguments))


def _resolve_bel_reaction(arguments: list) -> momapy_bel.core.Reaction:
    reactants = None
    products = None
    for argument in arguments:
        if isinstance(argument, _BELReactants):
            reactants = argument
        elif isinstance(argument, _BELProducts):
            products = argument
        else:
            raise ValueError(f"unexpected argument {argument!r} in rxn()")
    if reactants is None or products is None:
        raise ValueError("rxn() takes a reactants() and a products() argument")
    return momapy_bel.core.Reaction(
        reactants=frozenset(reactants.members),
        products=frozenset(products.members),
    )


def _resolve_bel_translocation(arguments: list) -> momapy_bel.core.Translocation:
    if not arguments:
        raise ValueError("tloc() takes at least one argument, got none")
    abundance = arguments[0]
    if not isinstance(abundance, momapy_bel.core.BELModelElement):
        raise ValueError(f"unexpected argument {abundance!r} in tloc()")
    from_location = None
    to_location = None
    # The one term whose resolver reshapes rather than maps: `fromLoc()` and
    # `toLoc()` are flattened into `Translocation`'s four scalar fields.
    for argument in arguments[1:]:
        if isinstance(argument, _BELFromLocation):
            from_location = argument
        elif isinstance(argument, _BELToLocation):
            to_location = argument
        else:
            raise ValueError(f"unexpected argument {argument!r} in tloc()")
    if from_location is None or to_location is None:
        raise ValueError("tloc() takes a fromLoc() and a toLoc() argument")
    return momapy_bel.core.Translocation(
        abundance=abundance,
        from_namespace=from_location.namespace,
        from_identifier=from_location.identifier,
        to_namespace=to_location.namespace,
        to_identifier=to_location.identifier,
    )


# Hardcoded from the BEL 2.1.3 spec rather than derived from any corpus: every
# function has both its short and its long spelling, and both map to the same
# resolver.
_BEL_FUNCTION_NAME_TO_RESOLVER = {
    "a": _resolve_bel_abundance,
    "abundance": _resolve_bel_abundance,
    "act": _resolve_bel_activity,
    "activity": _resolve_bel_activity,
    "bp": _resolve_bel_biological_process,
    "biologicalProcess": _resolve_bel_biological_process,
    "sec": _resolve_bel_cell_secretion,
    "cellSecretion": _resolve_bel_cell_secretion,
    "surf": _resolve_bel_cell_surface_expression,
    "cellSurfaceExpression": _resolve_bel_cell_surface_expression,
    "complex": _resolve_bel_complex_abundance,
    "complexAbundance": _resolve_bel_complex_abundance,
    "composite": _resolve_bel_composite_abundance,
    "compositeAbundance": _resolve_bel_composite_abundance,
    "deg": _resolve_bel_degradation,
    "degradation": _resolve_bel_degradation,
    "g": _resolve_bel_gene_abundance,
    "geneAbundance": _resolve_bel_gene_abundance,
    "gmod": _resolve_bel_gene_modification,
    "geneModification": _resolve_bel_gene_modification,
    "list": _resolve_bel_list,
    "m": _resolve_bel_microrna_abundance,
    "microRNAAbundance": _resolve_bel_microrna_abundance,
    "path": _resolve_bel_pathology,
    "pathology": _resolve_bel_pathology,
    "pop": _resolve_bel_population_abundance,
    "populationAbundance": _resolve_bel_population_abundance,
    "p": _resolve_bel_protein_abundance,
    "proteinAbundance": _resolve_bel_protein_abundance,
    "rxn": _resolve_bel_reaction,
    "reaction": _resolve_bel_reaction,
    "r": _resolve_bel_rna_abundance,
    "rnaAbundance": _resolve_bel_rna_abundance,
    "tloc": _resolve_bel_translocation,
    "translocation": _resolve_bel_translocation,
    "frag": _resolve_bel_fragment,
    "fragment": _resolve_bel_fragment,
    "fromLoc": _resolve_bel_from_location,
    "fus": _resolve_bel_fusion,
    "fusion": _resolve_bel_fusion,
    "loc": _resolve_bel_location,
    "location": _resolve_bel_location,
    "ma": _resolve_bel_molecular_activity,
    "molecularActivity": _resolve_bel_molecular_activity,
    "products": _resolve_bel_products,
    "pmod": _resolve_bel_protein_modification,
    "proteinModification": _resolve_bel_protein_modification,
    "reactants": _resolve_bel_reactants,
    "toLoc": _resolve_bel_to_location,
    "var": _resolve_bel_variant,
    "variant": _resolve_bel_variant,
}


momapy.io.core.register_reader("bel", BELReader)
momapy.io.core.register_writer("bel", BELWriter)
