from sqlalchemy import (
    Float,
    ForeignKey,
    ForeignKeyConstraint,
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


class ClinicalPatient(Base):
    __tablename__ = "clinical_patient"

    id: Mapped[str] = mapped_column(String(100), primary_key=True)
    dataset_id: Mapped[int] = mapped_column(
        Integer,
        ForeignKey("clinical_dataset.id"),
    )
    race: Mapped[str | None] = mapped_column(String(30), nullable=True)
    ethnicity: Mapped[str | None] = mapped_column(String(30), nullable=True)
    gender: Mapped[str | None] = mapped_column(String(1), nullable=True)
    sex_at_birth: Mapped[str | None] = mapped_column(String(1), nullable=True)


"""
Pre-clinical tables: tables related to the storage/organization of data from
in vitro datasets.

Expected extracted CSV shapes:

pre_clinical_cell_line.csv:
    cell_line_name,tissueid,mod_tissueid,accession,category,sex,age
    dataset_id is added during seeding from the dataset being loaded.

pre_clinical_sample.csv:
    sampleid,dataset_id,cell_line_name
    If an extraction script still writes id/cell_line_name, the seeding layer
    should map id -> sampleid and store cell_line_name directly.
    The composite foreign key (cell_line_name, dataset_id) references
    pre_clinical_cell_line(cell_line_name, dataset_id).

pre_clinical_treatment_response.csv:
    cell_line_name,treatment_id,ic50_recomputed,acc_recomputed,mechanism_of_action
    dataset_id is added during seeding from the dataset being loaded.
    cell_line_name is stored directly and enforced as a composite foreign key
    with dataset_id back to pre_clinical_cell_line(cell_line_name, dataset_id).

Molecular CSVs:
    pre_clinical_rna_seq.csv: sample_id,gene_id,expression_value
    pre_clinical_microarray.csv: sample_id,gene_id,expression_value
    pre_clinical_copy_number_variation.csv: sample_id,gene_id,value
    pre_clinical_mutation.csv: sample_id,gene_id,value

pre_clinical_gene.csv:
    id,name
"""


class PreClinicalDataset(Base):
    __tablename__ = "pre_clinical_dataset"

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    name: Mapped[str] = mapped_column(String(255), nullable=False, unique=True)

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
        ForeignKey("pre_clinical_dataset.id", ondelete="CASCADE"),
        nullable=False,
    )

    # Fields from pre_clinical_cell_line.csv.
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
    )

    dataset: Mapped["PreClinicalDataset"] = relationship(back_populates="cell_lines")
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


class PreClinicalSample(Base):
    __tablename__ = "pre_clinical_sample"

    # This is the final dataset-prefixed sample ID used by the extracted CSVs,
    # e.g. CCLE_<sampleid>, gcsi_<sampleid>, PRISM_<sampleid>.
    sampleid: Mapped[str] = mapped_column(String(255), primary_key=True)

    dataset_id: Mapped[int] = mapped_column(
        Integer,
        ForeignKey("pre_clinical_dataset.id", ondelete="CASCADE"),
        nullable=False,
    )

    # Store the human-readable cell-line key directly. The composite FK below
    # enforces that (cell_line_name, dataset_id) exists in pre_clinical_cell_line.
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
            "sampleid",
            "dataset_id",
            name="uq_pc_sample_dataset",
        ),
    )

    dataset: Mapped["PreClinicalDataset"] = relationship(
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
        ForeignKey("pre_clinical_dataset.id", ondelete="CASCADE"),
        nullable=False,
    )

    # Store the human-readable cell-line key directly. The composite FK below
    # enforces that (cell_line_name, dataset_id) exists in pre_clinical_cell_line.
    cell_line_name: Mapped[str] = mapped_column(String(255), nullable=False)

    treatment_id: Mapped[str] = mapped_column(String(255), nullable=False)
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
    )

    dataset: Mapped["PreClinicalDataset"] = relationship(
        back_populates="treatment_responses",
        overlaps="cell_line,treatment_responses",
    )
    cell_line: Mapped["PreClinicalCellLine"] = relationship(
        back_populates="treatment_responses",
        overlaps="dataset,treatment_responses",
    )


class PreClinicalGene(Base):
    __tablename__ = "pre_clinical_gene"

    id: Mapped[str] = mapped_column(String(255), primary_key=True)
    name: Mapped[str | None] = mapped_column(String(255), nullable=True)


class PreClinicalRnaSeq(Base):
    __tablename__ = "pre_clinical_rna_seq"

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    sample_id: Mapped[str] = mapped_column(
        String(255),
        ForeignKey("pre_clinical_sample.sampleid", ondelete="CASCADE"),
        nullable=False,
    )
    gene_id: Mapped[str] = mapped_column(String(255), nullable=False)
    expression_value: Mapped[float] = mapped_column(Float, nullable=False)

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
        ForeignKey("pre_clinical_sample.sampleid", ondelete="CASCADE"),
        nullable=False,
    )
    gene_id: Mapped[str] = mapped_column(String(255), nullable=False)
    expression_value: Mapped[float] = mapped_column(Float, nullable=False)

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
        ForeignKey("pre_clinical_sample.sampleid", ondelete="CASCADE"),
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
        ForeignKey("pre_clinical_sample.sampleid", ondelete="CASCADE"),
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