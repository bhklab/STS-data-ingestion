from sqlalchemy import String, Float, Integer, ForeignKey, Text, Boolean, DateTime, func, UniqueConstraint
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column, relationship
from datetime import datetime


# Base table to be referenced
class Base(DeclarativeBase):
    pass


"""
    Clinical tables: tables related to the storage/organization of data from in human in vivo datasets (patients in clinic)
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
    race: Mapped[str] = mapped_column(String(30), nullable=True)
    ethnicity: Mapped[str] = mapped_column(String(30), nullable=True)
    gender: Mapped[str] = mapped_column(String(1), nullable=True)
    sex_at_birth: Mapped[str] = mapped_column(String(1), nullable=True)


"""
    Pre clinical tables: tables related to the storage/organization of data from in vitro datasets
"""


class PreClinicalDataset(Base):
    __tablename__ = "pre_clinical_dataset"

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    name: Mapped[str] = mapped_column(String(255), nullable=False)

    cell_lines: Mapped[list["PreClinicalCellLine"]] = relationship(
        back_populates="dataset",
        cascade="all, delete-orphan",
    )


class PreClinicalCellLine(Base):
    __tablename__ = "pre_clinical_cell_line"

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    cell_line_name: Mapped[str] = mapped_column(String(255), nullable=False)
    dataset_id: Mapped[int] = mapped_column(
        Integer,
        ForeignKey("pre_clinical_dataset.id", ondelete="CASCADE"),
        nullable=False,
    )

    accession: Mapped[str | None] = mapped_column(String(255), nullable=True)
    category: Mapped[str | None] = mapped_column(String(255), nullable=True)
    sex: Mapped[str | None] = mapped_column(String(30), nullable=True)
    age: Mapped[int | None] = mapped_column(Integer, nullable=True)

    __table_args__ = (
        UniqueConstraint(
            "cell_line_name",
            "dataset_id",
            name="uq_cell_line_name_dataset_id",
        ),
    )

    dataset: Mapped["PreClinicalDataset"] = relationship(back_populates="cell_lines")
    treatment_responses: Mapped[list["PreClinicalTreatmentResponse"]] = relationship(
        back_populates="cell_line",
        cascade="all, delete-orphan",
    )
    samples: Mapped[list["Sample"]] = relationship(
        back_populates="cell_line",
        cascade="all, delete-orphan",
    )


class PreClinicalTreatmentResponse(Base):
    __tablename__ = "pre_clinical_treatment_response"

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    cell_line_id: Mapped[int] = mapped_column(
        Integer,
        ForeignKey("pre_clinical_cell_line.id", ondelete="CASCADE"),
        nullable=False,
    )

    treatment_id: Mapped[str] = mapped_column(String(255), nullable=False)
    ic50_recomputed: Mapped[float | None] = mapped_column(Float, nullable=True)
    acc_recomputed: Mapped[float | None] = mapped_column(Float, nullable=True)
    mechanism_of_action: Mapped[str | None] = mapped_column(String(255), nullable=True)

    __table_args__ = (
        UniqueConstraint(
            "cell_line_id",
            "treatment_id",
            name="uq_treatment_response_cell_line_treatment",
        ),
    )

    cell_line: Mapped["PreClinicalCellLine"] = relationship(back_populates="treatment_responses")


class PreClinicalSample(Base):
    __tablename__ = "pre_clinical_sample"

    id: Mapped[str] = mapped_column(String(255), primary_key=True)

    cell_line_id: Mapped[int | None] = mapped_column(
        Integer,
        ForeignKey("pre_clinical_cell_line.id", ondelete="CASCADE"),
        nullable=True,
    )

    site_primary: Mapped[str | None] = mapped_column(String(255), nullable=True)
    site_subtype1: Mapped[str | None] = mapped_column(String(255), nullable=True)
    site_subtype2: Mapped[str | None] = mapped_column(String(255), nullable=True)
    site_subtype3: Mapped[str | None] = mapped_column(String(255), nullable=True)

    histology: Mapped[str | None] = mapped_column(String(255), nullable=True)
    histology_subtype1: Mapped[str | None] = mapped_column(String(255), nullable=True)
    histology_subtype2: Mapped[str | None] = mapped_column(String(255), nullable=True)
    histology_subtype3: Mapped[str | None] = mapped_column(String(255), nullable=True)

    gender: Mapped[str | None] = mapped_column(String(30), nullable=True)
    age: Mapped[int | None] = mapped_column(Integer, nullable=True)
    race: Mapped[str | None] = mapped_column(String(30), nullable=True)
    diseases: Mapped[str | None] = mapped_column(Text, nullable=True)
    disease_type: Mapped[str | None] = mapped_column(String(128), nullable=True)

    cell_line: Mapped["PreClinicalCellLine"] = relationship(back_populates="samples")

class PreClinicalRnaSeq(Base):
    __tablename__ = "pre_clinical_rna_seq"

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    sample_id: Mapped[str] = mapped_column(
        String(255),
        ForeignKey("pre_clinical_sample.id", ondelete="CASCADE"),
        nullable=False,
    )

    gene_id: Mapped[str] = mapped_column(String(255), nullable=False)
    expression_value: Mapped[float] = mapped_column(Float, nullable=False)

    __table_args__ = (
        UniqueConstraint(
            "sample_id",
            "gene_id",
            name="uq_rna_seq_sample_gene",
        ),
    )

    sample: Mapped["PreClinicalSample"] = relationship(back_populates="rna_seq_data")

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
            name="uq_mutation_sample_gene",
        ),
    )

    sample: Mapped["PreClinicalSample"] = relationship(back_populates="mutations")

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
            name="uq_copy_number_variation_sample_gene",
        ),
    )

    sample: Mapped["PreClinicalSample"] = relationship(back_populates="copy_number_variations")

class 