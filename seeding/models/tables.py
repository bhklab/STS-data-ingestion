from sqlalchemy import (
    Boolean,
    DateTime,
    Float,
    ForeignKey,
    ForeignKeyConstraint,
    Index,
    Integer,
    String,
    Text,
    UniqueConstraint,
)
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column, relationship


class Base(DeclarativeBase):
    pass


"""
Clinical tables: tables related to the storage/organization of data from human
in vivo datasets (patients in clinic).
"""


class ClinicalDataset(Base):
    __tablename__ = "clinical_dataset"

    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    name: Mapped[str] = mapped_column(Text())


class ClinicalSample(Base):
    __tablename__ = "clinical_sample"

    id: Mapped[str] = mapped_column(String(100), primary_key=True)
    dataset_id: Mapped[int] = mapped_column(
        Integer,
        ForeignKey("datasets.id"),
    )
    race: Mapped[str | None] = mapped_column(String(30), nullable=True) # race_standardized
    ethnicity: Mapped[str | None] = mapped_column(String(30), nullable=True) # ethnicity_standardized
    sex: Mapped[str | None] = mapped_column(String(1), nullable=True) # sex_curated
    age: Mapped[int | None] = mapped_column(Integer, nullable=True) # age_curated
    histology: Mapped[str | None] = mapped_column(String(50), nullable=True) # histo
    tissue: Mapped[str | None] = mapped_column(String(50), nullable=True) # cancer_type
    tissue_origin: Mapped[str | None] = mapped_column(String(50), nullable=True) # tissue_or_organ_of_origin
    

class ClinicalAntigen(Base):
    __tablename__ = "clinical_antigen"

    id: Mapped[str] = mapped_column(
        String(255),
        primary_key=True,
    )
    catalogue_number: Mapped[str | None] = mapped_column(
        String(255),
        nullable=True,
    )
    peptide_target: Mapped[str | None] = mapped_column(
        String(255),
        nullable=True,
    )
    peptide_target_gene: Mapped[str | None] = mapped_column(
        String(255),
        ForeignKey("pre_clinical_gene.id"),
        nullable=True,
    )


class ClinicalProbe(Base):
    __tablename__ = "clinical_probe"

    id: Mapped[str] = mapped_column(
        String(255),
        primary_key=True,
    )
    name: Mapped[str | None] = mapped_column(
        String(255),
        nullable=True,
    )


class ClinicalRNA(Base):
    __tablename__ = "clinical_rna"

    gene_id: Mapped[str] = mapped_column(
        String(255),
        ForeignKey("pre_clinical_gene.id"),
        primary_key=True,
    )
    sample_id: Mapped[str] = mapped_column(
        String(255),
        ForeignKey("clinical_sample.id", ondelete="CASCADE"),
        primary_key=True,
    )
    value: Mapped[float | None] = mapped_column(
        Float,
        nullable=True,
    ) # tpm_unstrand matrix


class ClinicalMutation(Base):
    __tablename__ = "clinical_mutation"

    gene_id: Mapped[str] = mapped_column(
        String(255),
        ForeignKey("pre_clinical_gene.id"),
        primary_key=True,
    )
    sample_id: Mapped[str] = mapped_column(
        String(255),
        ForeignKey("clinical_sample.id", ondelete="CASCADE"),
        primary_key=True,
    )
    mutation: Mapped[Boolean | None] = mapped_column(
        Boolean,
        nullable=True,
    ) # mutation matrix
    oncoprint: Mapped[str | None] = mapped_column(
        String(255),
        nullable=True,
    ) # oncoprint matrix
        
    


class ClinicalCNV(Base):
    __tablename__ = "clinical_cnv"

    gene_id: Mapped[str] = mapped_column(
        String(255),
        ForeignKey("pre_clinical_gene.id"),
        primary_key=True,
    )
    sample_id: Mapped[str] = mapped_column(
        String(255),
        ForeignKey("clinical_sample.id", ondelete="CASCADE"),
        primary_key=True,
    )
    value: Mapped[int | None] = mapped_column(
        Integer,
        nullable=True,
    ) # copy_number


class ClinicalRPPA(Base):
    __tablename__ = "clinical_rppa"

    antigen_id: Mapped[str] = mapped_column(
        String(255),
        ForeignKey("clinical_antigen.id"),
        primary_key=True,
    )
    sample_id: Mapped[str] = mapped_column(
        String(255),
        ForeignKey("clinical_sample.id", ondelete="CASCADE"),
        primary_key=True,
    )
    value: Mapped[float | None] = mapped_column(
        Float,
        nullable=True,
    ) # expression matrix


class ClinicalMiRNA(Base):
    __tablename__ = "clinical_mirna"

    id: Mapped[str] = mapped_column(
        String(255),
        primary_key=True,
    )
    sample_id: Mapped[str] = mapped_column(
        String(255),
        ForeignKey("clinical_sample.id", ondelete="CASCADE"),
        primary_key=True,
    )
    gene_id: Mapped[str | None] = mapped_column(
        String(255),
        ForeignKey("pre_clinical_gene.id"),
        nullable=True,
    )
    value: Mapped[float | None] = mapped_column(
        Float,
        nullable=True,
    ) #rpm matrix


class ClinicalMethylation(Base):
    __tablename__ = "clinical_methylation"

    probe_id: Mapped[str] = mapped_column(
        String(255),
        ForeignKey("clinical_probe.id"),
        primary_key=True,
    )
    sample_id: Mapped[str] = mapped_column(
        String(255),
        ForeignKey("clinical_sample.id", ondelete="CASCADE"),
        primary_key=True,
    )
    value: Mapped[float | None] = mapped_column(
        Float,
        nullable=True,
    ) #listData [[1]] matrix


"""
Pre-clinical tables: tables related to the storage/organization of data from
in vitro datasets.

Expected extracted CSV shapes:

pre_clinical_cell_line.csv:
    cell_line_name,tissueid,mod_tissueid,accession,category,sex,age
    dataset_id is added during seeding from the dataset being loaded.

pre_clinical_sample.csv:
    id,cell_line_name
    id is the final dataset-prefixed sample ID, e.g. CCLE_A204..., gcsi_..., PRISM_...
    The composite foreign key (cell_line_name, dataset_id) references
    pre_clinical_cell_line(cell_line_name, dataset_id).

pre_clinical_treatment_response.csv:
    cell_line_name,treatment_id,cid,ic50_recomputed,acc_recomputed,mechanism_of_action
    dataset_id is added during seeding from the dataset being loaded.
    cid is the compound identifier from the PharmacoSet treatment slot.

Molecular CSVs:
    pre_clinical_rna_seq.csv: sample_id,gene_id,value
    pre_clinical_microarray.csv: sample_id,gene_id,value
    pre_clinical_copy_number_variation.csv: sample_id,gene_id,value
    pre_clinical_mutation.csv: sample_id,gene_id,value

pre_clinical_gene.csv:
    id,name
    Ensembl version suffixes such as .20 or .20_PAR_Y should be stripped by the
    extraction/seeding layer before insertion.
"""


class Dataset(Base):
    __tablename__ = "datasets"

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    name: Mapped[str] = mapped_column(String(255), nullable=False, unique=True)
    version: Mapped[str | None] = mapped_column(String(50), nullable=True)
    software: Mapped[str | None] = mapped_column(String(255), nullable=True)
    link: Mapped[str | None] = mapped_column(Text, nullable=True)
    publication: Mapped[str | None] = mapped_column(String(255), nullable=True)
    PMID: Mapped[str | None] = mapped_column(String(255), nullable=True)
    description: Mapped[str | None] = mapped_column(Text, nullable=True)
    key_study_findings: Mapped[str | None] = mapped_column(Text, nullable=True)
    clinical: Mapped[bool | None] = mapped_column(Boolean, nullable=True)

    cell_lines: Mapped[list["PreClinicalCellLine"]] = relationship(
        back_populates="dataset",
        cascade="all, delete-orphan",
    )
    samples: Mapped[list["PreClinicalSample"]] = relationship(
        back_populates="dataset",
        cascade="all, delete-orphan",
        overlaps="cell_line,samples",
    )
    treatment_responses: Mapped[list["PreClinicalTreatmentResponse"]] = relationship(
        back_populates="dataset",
        cascade="all, delete-orphan",
        overlaps="cell_line,treatment_responses",
    )


class PreClinicalCellLine(Base):
    __tablename__ = "pre_clinical_cell_line"

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    dataset_id: Mapped[int] = mapped_column(
        Integer,
        ForeignKey("datasets.id", ondelete="CASCADE"),
        nullable=False,
    )

    cell_line_name: Mapped[str] = mapped_column(String(255), nullable=False)
    tissueid: Mapped[str | None] = mapped_column(String(255), nullable=True)
    mod_tissueid: Mapped[str | None] = mapped_column(String(255), nullable=True)
    accession: Mapped[str | None] = mapped_column(String(255), nullable=True)
    category: Mapped[str | None] = mapped_column(String(255), nullable=True)
    sex: Mapped[str | None] = mapped_column(String(30), nullable=True)
    age: Mapped[int | None] = mapped_column(Integer, nullable=True)

    __table_args__ = (
        UniqueConstraint(
            "cell_line_name",
            "dataset_id",
            name="uq_pc_cell_line_dataset",
        ),
        Index("ix_pc_cell_line_accession", "accession"),
    )

    dataset: Mapped["Dataset"] = relationship(back_populates="cell_lines")
    samples: Mapped[list["PreClinicalSample"]] = relationship(
        back_populates="cell_line",
        cascade="all, delete-orphan",
        overlaps="dataset,samples",
    )
    treatment_responses: Mapped[list["PreClinicalTreatmentResponse"]] = relationship(
        back_populates="cell_line",
        cascade="all, delete-orphan",
        overlaps="dataset,treatment_responses",
    )


class PreClinicalCellLineInfo(Base):
    __tablename__ = "pre_clinical_cell_line_info"

    accession: Mapped[str] = mapped_column(
        String(255),
        ForeignKey("pre_clinical_cell_line.accession"),
        primary_key=True,
    )
    cell_line_name: Mapped[str | None] = mapped_column(String(255), nullable=True)
    category: Mapped[str | None] = mapped_column(String(255), nullable=True)
    date: Mapped[str | None] = mapped_column(String(255), nullable=True)
    age_at_sampling: Mapped[str | None] = mapped_column(String(255), nullable=True)
    sex_of_cell: Mapped[str | None] = mapped_column(String(255), nullable=True)
    hierarchy: Mapped[str | None] = mapped_column(Text, nullable=True)
    cell_type: Mapped[str | None] = mapped_column(Text, nullable=True)
    derived_from_site: Mapped[str | None] = mapped_column(Text, nullable=True)
    donor_information: Mapped[str | None] = mapped_column(Text, nullable=True)
    doubling_time: Mapped[str | None] = mapped_column(Text, nullable=True)
    genome_ancestry: Mapped[str | None] = mapped_column(Text, nullable=True)
    hla_typing: Mapped[str | None] = mapped_column(Text, nullable=True)
    microsatellite_instability: Mapped[str | None] = mapped_column(Text, nullable=True)
    omics: Mapped[str | None] = mapped_column(Text, nullable=True)
    part_of: Mapped[str | None] = mapped_column(Text, nullable=True)
    population: Mapped[str | None] = mapped_column(Text, nullable=True)
    sequence_variation: Mapped[str | None] = mapped_column(Text, nullable=True)
    anecdotal: Mapped[str | None] = mapped_column(Text, nullable=True)
    biotechnology: Mapped[str | None] = mapped_column(Text, nullable=True)
    discontinued: Mapped[str | None] = mapped_column(Text, nullable=True)
    group_col: Mapped[str | None] = mapped_column(Text, nullable=True)
    misspelling: Mapped[str | None] = mapped_column(Text, nullable=True)
    registration: Mapped[str | None] = mapped_column(Text, nullable=True)
    virology: Mapped[str | None] = mapped_column(Text, nullable=True)
    caution: Mapped[str | None] = mapped_column(Text, nullable=True)
    characteristics: Mapped[str | None] = mapped_column(Text, nullable=True)
    karyotypic_information: Mapped[str | None] = mapped_column(Text, nullable=True)
    problematic_cell_line: Mapped[str | None] = mapped_column(Text, nullable=True)
    transformant: Mapped[str | None] = mapped_column(Text, nullable=True)
    miscellaneous: Mapped[str | None] = mapped_column(Text, nullable=True)
    from_col: Mapped[str | None] = mapped_column(Text, nullable=True)
    genetic_integration: Mapped[str | None] = mapped_column(Text, nullable=True)
    knockout_cell: Mapped[str | None] = mapped_column(Text, nullable=True)
    selected_for_resistance_to: Mapped[str | None] = mapped_column(Text, nullable=True)

    # Flattened from the AnnotationDB diseases list.
    disease_ids: Mapped[str | None] = mapped_column(Text, nullable=True)
    disease_descriptions: Mapped[str | None] = mapped_column(Text, nullable=True)

    # OncoTree-derived fields. `first_level` is the level-2 parent code with
    # underscores replaced by spaces, e.g. SOFT_TISSUE -> SOFT TISSUE.
    # `second_level` is the OncoTree level-2 tumor-type name.
    first_level: Mapped[str | None] = mapped_column(String(255), nullable=True)
    second_level: Mapped[str | None] = mapped_column(String(255), nullable=True)


class PreClinicalSample(Base):
    __tablename__ = "pre_clinical_sample"

    id: Mapped[str] = mapped_column(String(255), primary_key=True)
    dataset_id: Mapped[int] = mapped_column(
        Integer,
        ForeignKey("datasets.id", ondelete="CASCADE"),
        nullable=False,
    )
    cell_line_name: Mapped[str] = mapped_column(String(255), nullable=False)

    __table_args__ = (
        ForeignKeyConstraint(
            ["cell_line_name", "dataset_id"],
            [
                "pre_clinical_cell_line.cell_line_name",
                "pre_clinical_cell_line.dataset_id",
            ],
            ondelete="CASCADE",
            name="fk_pc_sample_cell_line_dataset",
        ),
        UniqueConstraint(
            "id",
            "dataset_id",
            name="uq_pc_sample_dataset",
        ),
    )

    dataset: Mapped["Dataset"] = relationship(
        back_populates="samples",
        overlaps="cell_line,samples",
    )
    cell_line: Mapped["PreClinicalCellLine"] = relationship(
        back_populates="samples",
        overlaps="dataset,samples,cell_lines",
    )

    rna_seq_data: Mapped[list["PreClinicalRnaSeq"]] = relationship(
        back_populates="sample",
        cascade="all, delete-orphan",
    )
    microarray_data: Mapped[list["PreClinicalMicroarray"]] = relationship(
        back_populates="sample",
        cascade="all, delete-orphan",
    )
    copy_number_variations: Mapped[list["PreClinicalCopyNumberVariation"]] = relationship(
        back_populates="sample",
        cascade="all, delete-orphan",
    )
    mutations: Mapped[list["PreClinicalMutation"]] = relationship(
        back_populates="sample",
        cascade="all, delete-orphan",
    )


class PreClinicalTreatmentResponse(Base):
    __tablename__ = "pre_clinical_treatment_response"

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    dataset_id: Mapped[int] = mapped_column(
        Integer,
        ForeignKey("datasets.id", ondelete="CASCADE"),
        nullable=False,
    )
    cell_line_name: Mapped[str] = mapped_column(String(255), nullable=False)
    treatment_id: Mapped[str] = mapped_column(String(255), nullable=False)
    cid: Mapped[str | None] = mapped_column(String(255), nullable=True)
    ic50_recomputed: Mapped[float | None] = mapped_column(Float, nullable=True)
    acc_recomputed: Mapped[float | None] = mapped_column(Float, nullable=True)
    mechanism_of_action: Mapped[str | None] = mapped_column(String(255), nullable=True)

    __table_args__ = (
        ForeignKeyConstraint(
            ["cell_line_name", "dataset_id"],
            [
                "pre_clinical_cell_line.cell_line_name",
                "pre_clinical_cell_line.dataset_id",
            ],
            ondelete="CASCADE",
            name="fk_pc_tr_cell_line_dataset",
        ),
        UniqueConstraint(
            "dataset_id",
            "cell_line_name",
            "treatment_id",
            name="uq_pc_tr_dataset_cell_line_treatment",
        ),
        Index("ix_pc_tr_cid", "cid"),
    )

    dataset: Mapped["Dataset"] = relationship(
        back_populates="treatment_responses",
        overlaps="cell_line,treatment_responses",
    )
    cell_line: Mapped["PreClinicalCellLine"] = relationship(
        back_populates="treatment_responses",
        overlaps="dataset,treatment_responses",
    )


class PreClinicalDrug(Base):
    __tablename__ = "drugs"

    # PubChem CID from pre_clinical_treatment_response.cid.
    # This FK intentionally follows the requested direction: drugs.cid references
    # treatment-response cid values that already exist in the loaded assay data.
    cid: Mapped[str] = mapped_column(
        String(255),
        ForeignKey("pre_clinical_treatment_response.cid", ondelete="CASCADE"),
        primary_key=True,
    )

    title: Mapped[str | None] = mapped_column(String(255), nullable=True)
    mapped_name: Mapped[str | None] = mapped_column(Text, nullable=True)
    molecule_chembl_id: Mapped[str | None] = mapped_column(String(255), nullable=True)
    molecule_chembl_id_from_synonyms: Mapped[bool | None] = mapped_column(Boolean, nullable=True)
    molecular_formula: Mapped[str | None] = mapped_column(String(255), nullable=True)
    molecular_weight: Mapped[str | None] = mapped_column(String(255), nullable=True)
    smiles: Mapped[str | None] = mapped_column(Text, nullable=True)
    connectivity_smiles: Mapped[str | None] = mapped_column(Text, nullable=True)
    inchi: Mapped[str | None] = mapped_column(Text, nullable=True)
    inchikey: Mapped[str | None] = mapped_column(String(255), nullable=True)
    iupac_name: Mapped[str | None] = mapped_column(Text, nullable=True)
    xlogp: Mapped[float | None] = mapped_column(Float, nullable=True)
    exact_mass: Mapped[str | None] = mapped_column(String(255), nullable=True)
    monoisotopic_mass: Mapped[str | None] = mapped_column(String(255), nullable=True)
    tpsa: Mapped[float | None] = mapped_column(Float, nullable=True)
    complexity: Mapped[float | None] = mapped_column(Float, nullable=True)
    charge: Mapped[int | None] = mapped_column(Integer, nullable=True)
    h_bond_donor_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    h_bond_acceptor_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    rotatable_bond_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    heavy_atom_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    isotope_atom_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    atom_stereo_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    defined_atom_stereo_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    undefined_atom_stereo_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    bond_stereo_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    defined_bond_stereo_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    undefined_bond_stereo_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    covalent_unit_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    volume_3d: Mapped[float | None] = mapped_column(Float, nullable=True)
    x_steric_quadrupole_3d: Mapped[float | None] = mapped_column(Float, nullable=True)
    y_steric_quadrupole_3d: Mapped[float | None] = mapped_column(Float, nullable=True)
    z_steric_quadrupole_3d: Mapped[float | None] = mapped_column(Float, nullable=True)
    feature_count_3d: Mapped[int | None] = mapped_column(Integer, nullable=True)
    feature_acceptor_count_3d: Mapped[int | None] = mapped_column(Integer, nullable=True)
    feature_donor_count_3d: Mapped[int | None] = mapped_column(Integer, nullable=True)
    feature_anion_count_3d: Mapped[int | None] = mapped_column(Integer, nullable=True)
    feature_cation_count_3d: Mapped[int | None] = mapped_column(Integer, nullable=True)
    feature_ring_count_3d: Mapped[int | None] = mapped_column(Integer, nullable=True)
    feature_hydrophobe_count_3d: Mapped[int | None] = mapped_column(Integer, nullable=True)
    conformer_model_rmsd_3d: Mapped[float | None] = mapped_column(Float, nullable=True)
    effective_rotor_count_3d: Mapped[float | None] = mapped_column(Float, nullable=True)
    conformer_count_3d: Mapped[int | None] = mapped_column(Integer, nullable=True)
    fingerprint_2d: Mapped[str | None] = mapped_column(Text, nullable=True)
    patent_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    patent_family_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    literature_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    annotation_types: Mapped[str | None] = mapped_column(Text, nullable=True)
    annotation_type_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    chembl_max_phase: Mapped[int | None] = mapped_column(Integer, nullable=True)
    drug_like: Mapped[bool | None] = mapped_column(Boolean, nullable=True)
    fda_approval: Mapped[bool | None] = mapped_column(Boolean, nullable=True)
    date_added: Mapped[object | None] = mapped_column(DateTime, nullable=True)
    atc_code: Mapped[str | None] = mapped_column(String(255), nullable=True)

    # First object from the AnnotationDB mechanisms list, if present.
    mechanism_molecule_chembl_id: Mapped[str | None] = mapped_column(String(255), nullable=True)
    mechanism_parent_molecule_chembl_id: Mapped[str | None] = mapped_column(String(255), nullable=True)
    mechanism_action_type: Mapped[str | None] = mapped_column(String(255), nullable=True)
    mechanism_of_action: Mapped[str | None] = mapped_column(Text, nullable=True)


class PreClinicalGene(Base):
    __tablename__ = "pre_clinical_gene"

    id: Mapped[str] = mapped_column(String(255), primary_key=True)
    name: Mapped[str | None] = mapped_column(String(255), nullable=True)


class PreClinicalRnaSeq(Base):
    __tablename__ = "pre_clinical_rna_seq"

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    sample_id: Mapped[str] = mapped_column(
        String(255),
        ForeignKey("pre_clinical_sample.id", ondelete="CASCADE"),
        nullable=False,
    )
    gene_id: Mapped[str] = mapped_column(String(255), nullable=False)
    value: Mapped[float] = mapped_column(Float, nullable=False)

    __table_args__ = (
        UniqueConstraint(
            "sample_id",
            "gene_id",
            name="uq_pc_rna_sample_gene",
        ),
    )

    sample: Mapped["PreClinicalSample"] = relationship(back_populates="rna_seq_data")


class PreClinicalMicroarray(Base):
    __tablename__ = "pre_clinical_microarray"

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    sample_id: Mapped[str] = mapped_column(
        String(255),
        ForeignKey("pre_clinical_sample.id", ondelete="CASCADE"),
        nullable=False,
    )
    gene_id: Mapped[str] = mapped_column(String(255), nullable=False)
    value: Mapped[float] = mapped_column(Float, nullable=False)

    __table_args__ = (
        UniqueConstraint(
            "sample_id",
            "gene_id",
            name="uq_pc_microarray_sample_gene",
        ),
    )

    sample: Mapped["PreClinicalSample"] = relationship(back_populates="microarray_data")


class PreClinicalCopyNumberVariation(Base):
    __tablename__ = "pre_clinical_copy_number_variation"

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    sample_id: Mapped[str] = mapped_column(
        String(255),
        ForeignKey("pre_clinical_sample.id", ondelete="CASCADE"),
        nullable=False,
    )
    gene_id: Mapped[str] = mapped_column(String(255), nullable=False)
    value: Mapped[float] = mapped_column(Float, nullable=False)

    __table_args__ = (
        UniqueConstraint(
            "sample_id",
            "gene_id",
            name="uq_pc_cnv_sample_gene",
        ),
    )

    sample: Mapped["PreClinicalSample"] = relationship(
        back_populates="copy_number_variations"
    )


class PreClinicalMutation(Base):
    __tablename__ = "pre_clinical_mutation"

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    sample_id: Mapped[str] = mapped_column(
        String(255),
        ForeignKey("pre_clinical_sample.id", ondelete="CASCADE"),
        nullable=False,
    )
    gene_id: Mapped[str] = mapped_column(String(255), nullable=False)
    value: Mapped[str] = mapped_column(String(255), nullable=False)

    __table_args__ = (
        UniqueConstraint(
            "sample_id",
            "gene_id",
            name="uq_pc_mutation_sample_gene",
        ),
    )

    sample: Mapped["PreClinicalSample"] = relationship(back_populates="mutations")