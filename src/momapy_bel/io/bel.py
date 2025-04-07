import urllib.parse

import momapy.io
import momapy_bel.core


class BELWriter(momapy.io.Writer):
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
        momapy_bel.core.TranslatedTo: "_translated_to",
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
        obj: momapy_bel.core.BELModel,
        file_path: str,
        namespace_definitions,
        annotation_definitions,
        annotations,
        with_abundances_as_statements=False,
        with_biological_processes_as_statements=False,
        with_reactions_as_statements=True,
        with_degradations_as_statements=True,
    ):
        bel_string = cls._bel_model_to_string(
            obj,
            namespace_definitions,
            annotation_definitions,
            annotations,
            with_abundances_as_statements=with_abundances_as_statements,
            with_biological_processes_as_statements=with_biological_processes_as_statements,
            with_reactions_as_statements=with_reactions_as_statements,
            with_degradations_as_statements=with_degradations_as_statements,
        )
        with open(file_path, "w") as f:
            f.write(bel_string)

    @classmethod
    def _list_to_string(cls, bel_list):
        args = [
            cls._bel_element_to_string(element)
            for element in bel_list.elements
        ]
        return cls._make_function_string("list", args)

    @classmethod
    def _location_to_string(cls, location):
        args = [
            cls._make_namespace_identifier_arg(
                location.namespace, location.identifier
            )
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
            args.append(
                cls._bel_element_to_string(activity.molecular_activity)
            )
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
        args = [fragment.start_stop]
        if fragment.descriptor is not None:
            args.append(fragment.descriptor)
        return cls._make_function_string("frag", args)

    @classmethod
    def _fusion_to_string(cls, fusion):
        args = [
            cls._make_namespace_identifier_arg(
                fusion.namespace5, fusion.identifier5
            ),
            fusion.range5,
            cls._make_namespace_identifier_arg(
                fusion.namespace3, fusion.identifier3
            ),
            fusion.range3,
        ]
        return cls._make_namespace_identifier_arg("fus", args)

    @classmethod
    def _variant_to_string(cls, variant):
        args = [variant.descriptor]
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
        if gene_abundance.variant is not None:
            args.append(cls._bel_element_to_string(gene_abundance.variant))
        return cls._make_function_string("g", args)

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
        if microrna.variant is not None:
            args.append(cls._bel_element_to_string(microrna.variant))
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
            args.append(
                cls._bel_element_to_string(population_abundance.location)
            )
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
        args = [
            cls._make_namespace_identifier_arg(
                protein_abundance.namespace, protein_abundance.identifier
            )
        ]
        if protein_abundance.location is not None:
            args.append(cls._bel_element_to_string(protein_abundance.location))
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
        if (
            rna_abundance.namespace is not None
            and rna_abundance.identifier is not None
        ):
            args = [
                cls._make_namespace_identifier_arg(
                    rna_abundance.namespace, rna_abundance.identifier
                )
            ]
        else:
            args = [cls._bel_element_to_string(rna_abundance.fusion)]
        if rna_abundance.location is not None:
            args.append(cls._bel_element_to_string(rna_abundance.location))
        if rna_abundance.variant is not None:
            args.append(cls._bel_element_to_string(rna_abundance.variant))
        return cls._make_function_string("r", args)

    @classmethod
    def _translocation_to_string(cls, translocation):
        args = [
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
    def _generic_annotation_to_set_string(cls, annotation):
        return cls._make_set_string(
            annotation.definition.name, annotation.args
        )

    @classmethod
    def _generic_annotation_to_unset_string(cls, annotation):
        return cls._make_unset_string(annotation.definition.name)

    @classmethod
    def _is_url(cls, s):
        result = urllib.parse.urlparse(s)
        if result.scheme != "":
            return True
        return False

    @classmethod
    def _get_as_type_from_as(cls, as_):
        if isinstance(as_, list):
            as_type = "LIST"
        elif cls._is_url(as_):
            as_type = "URL"
        else:
            as_type = "PATTERN"
        return as_type

    @classmethod
    def _annotation_definition_to_define_string(cls, annotation_definition):
        as_type = cls._get_as_type_from_as(annotation_definition.as_)
        if isinstance(annotation_definition.as_, str):
            as_ = [annotation_definition.as_]
        else:
            as_ = annotation_definition.as_
        return cls._make_define_string(
            "ANNOTATION", annotation_definition.name, as_type, as_
        )

    @classmethod
    def _namespace_definition_to_define_string(cls, namespace_definition):
        as_type = cls._get_as_type_from_as(namespace_definition.as_)
        if isinstance(namespace_definition.as_, str):
            as_ = [namespace_definition.as_]
        else:
            as_ = namespace_definition.as_
        return cls._make_define_string(
            "NAMESPACE", namespace_definition.name, as_type, as_
        )

    @classmethod
    def _document_annotation_to_set_string(cls, document_annotation):
        output_strings = []
        if document_annotation.name is not None:
            output_strings.append(
                cls._make_set_string(
                    "DOCUMENT Name", [document_annotation.name]
                )
            )
        if document_annotation.description is not None:
            output_strings.append(
                cls._make_set_string(
                    "DOCUMENT Description", [document_annotation.description]
                )
            )
        if document_annotation.authors is not None:
            output_strings.append(
                cls._make_set_string(
                    "DOCUMENT Authors", [document_annotation.authors]
                )
            )
        return "\n".join(output_strings)

    @classmethod
    def _bel_model_to_string(
        cls,
        bel_model,
        bel_namespace_definitions,
        bel_annotation_definitions,
        bel_annotations,
        with_abundances_as_statements=False,
        with_biological_processes_as_statements=False,
        with_reactions_as_statements=True,
        with_degradations_as_statements=True,
    ):
        output_strings = []
        bel_model_annotations = bel_annotations.get(bel_model)
        if bel_model_annotations is not None:
            for bel_model_annotation in bel_model_annotations:
                output_strings.append(
                    cls._bel_annotation_to_string(bel_model_annotation)
                )
        for bel_namespace_definition in bel_namespace_definitions:
            define_string = cls._namespace_definition_to_define_string(
                bel_namespace_definition
            )
            output_strings.append(define_string)
        for bel_annotation_definition in bel_annotation_definitions:
            define_string = cls._annotation_definition_to_define_string(
                bel_annotation_definition
            )
            output_strings.append(define_string)
        for bel_statement in bel_model.statements:
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
                    or not isinstance(
                        bel_statement, momapy_bel.core.Degradation
                    )
                )
            ):
                unset_strings = []
                bel_statement_annotations = bel_annotations.get(bel_statement)
                if bel_statement_annotations is not None:
                    for bel_annotation in bel_statement_annotations:
                        set_string = cls._bel_annotation_to_string(
                            bel_annotation, set_or_unset="set"
                        )
                        output_strings.append(set_string)
                        unset_string = cls._bel_annotation_to_string(
                            bel_annotation, set_or_unset="unset"
                        )
                        unset_strings.append(unset_string)
                output_strings.append(
                    cls._bel_element_to_string(bel_statement)
                )
                for unset_string in unset_strings:
                    output_strings.append(unset_string)
        bel_string = "\n".join(output_strings)
        return bel_string

    @classmethod
    def _make_namespace_identifier_arg(cls, namespace, identifier):
        if not identifier.isalnum():
            identifier = f'"{identifier}"'
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
        args = [f'"{arg}"' for arg in args]
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


momapy.io.register_writer("bel", BELWriter)
