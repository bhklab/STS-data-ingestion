from sqlalchemy import String, Float, Integer, ForeignKey, Text, Boolean, DateTime, func
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

# Patient/Sample/Outcomes Tables

class ClinicalPatient(Base):
    __tablename__ = "clinical_patient"

    id: Mapped[str] = mapped_column(String(100), primary_key=True)
    dataset_id: Mapped[int] = mapped_column(
        Integer,
        ForeignKey("clinical_dataset.id"),
    )
    race: Mapped[str | None] = mapped_column(String(30), nullable=True)
    ethnicity: Mapped[str | None] = mapped_column(String(30), nullable=True)
    gender: Mapped[str | None] = mapped_column(String(1), nullable=True) # If the clinical gives full male/female, map to M/F
    sex_at_birth: Mapped[str | None] = mapped_column(String(1), nullable=True) # If the clinical gives full male/female, map to M/F

class ClinicalSample(Base):
	__tablename__ = "clinical_sample"

	id: Mapped[str] = mapped_column(String(100), primary_key=True)
	patient_id: Mapped[str] = mapped_column(String(100), ForeignKey("clinical_patient.id"), nullable=False)
	age_at_index: Mapped[int | None] = mapped_column(Integer, nullable=True)
	days_to_death: Mapped[int | None] = mapped_column(Integer, nullable=True)
	disease_type: Mapped[str | None] = mapped_column(String(100), nullable=True)
	primary_diagnosis_l1: Mapped[str | None] = mapped_column(String(100), nullable=True) # Primary diagnosis L1 to be retrieved from Oncotree (primary_diagnosis is already a base field in TCGA)
	primary_diagnosis_l2: Mapped[str | None] = mapped_column(String(100), nullable=True) # Primary diagnosis L2 to be retrieved from Oncotree (primary_diagnosis is already a base field in TCGA)
	tumor_definition: Mapped[str | None] = mapped_column(String(100), nullable=True)
	classification_of_tumor: Mapped[str | None] = mapped_column(String(100), nullable=True)

class ClinicalOutcome(Base):
	__tablename__ = "clinical_outcome"

	patient_id: Mapped[str] = mapped_column(String(100), ForeignKey("clinical_patient.id"), primary_key=True, nullable=False)
	recist: Mapped[str | None] = mapped_column(String(20), nullable=True)
	response: Mapped[str | None] = mapped_column(String(20), nullable=True)
	# Survival times (in days)
	survival_time_os: Mapped[int | None] = mapped_column(Integer, nullable=True)
	survival_time_pfs: Mapped[int | None] = mapped_column(Integer, nullable=True)
	survival_time_rfs: Mapped[int | None] = mapped_column(Integer, nullable=True)
	# Events (binary)
	event_occured_os: Mapped[bool | None] = mapped_column(Boolean, nullable=True)
	event_occured_pfs: Mapped[bool | None] = mapped_column(Boolean, nullable=True)
	event_occured_rfs: Mapped[bool | None] = mapped_column(Boolean, nullable=True)




# Assay Layers

class ClinicalRNA(Base):
	__tablename__ = "clinical_rna"

	gene_id: Mapped[str] = mapped_column(String(100), ForeignKey("gene.id"), primary_key=True)
	sample_id: Mapped[str] = mapped_column(String(100), ForeignKey("clinical_sample.id"), primary_key=True)
	value: Mapped[float] = mapped_column(Float(), nullable=False) # 0.0+

class ClinicalMutation(Base):
	__tablename__ = "clinical_mutation"

	gene_id: Mapped[str] = mapped_column(String(100), ForeignKey("gene.id"), primary_key=True)
	sample_id: Mapped[str] = mapped_column(String(100), ForeignKey("clinical_sample.id"), primary_key=True)
	value: Mapped[int] = mapped_column(Integer, nullable=False) # 0 or 1

class ClinicalCNV(Base):
	__tablename__ = "clinical_cnv"

	gene_id: Mapped[str] = mapped_column(String(100), ForeignKey("gene.id"), primary_key=True)
	sample_id: Mapped[str] = mapped_column(String(100), ForeignKey("clinical_sample.id"), primary_key=True)
	value: Mapped[int] = mapped_column(Integer, nullable=False) # 1+

class ClinicalRPPA(Base):
	__tablename__ = "clinical_rppa"

	antigen_id: Mapped[str] = mapped_column(String(100), ForeignKey("antigen.id"), primary_key=True) # Antigen identifier, Ex. AGID00100
	sample_id: Mapped[str] = mapped_column(String(100), ForeignKey("clinical_sample.id"), primary_key=True)
	value: Mapped[float] = mapped_column(Float(), nullable=False) # negative or positive float

class ClinicalMirna(Base):
	__tablename__ = "clinical_mirna"

	mirna_id: Mapped[str] = mapped_column(String(100), ForeignKey("mirna.id"), primary_key=True) # MiRBase ID, Ex. hsa-let-7a-1 (found as a synonym to an ENSEMBL gene)
	sample_id: Mapped[str] = mapped_column(String(100), ForeignKey("clinical_sample.id"), primary_key=True)
	gene_id: Mapped[str] = mapped_column(String(100), ForeignKey("gene.id"), nullable=False) # extracted from id (Ex. id --> ENSEMBL))
	value: Mapped[float] = mapped_column(Float(), nullable=False) # 0.0+

class ClinicalMethylation(Base):
	__tablename__ = "clinical_mirna"

	probe_id: Mapped[str] = mapped_column(String(100), ForeignKey("mirna.id"), primary_key=True) # MiRBase ID, Ex. hsa-let-7a-1 (found as a synonym to an ENSEMBL gene)
	sample_id: Mapped[str] = mapped_column(String(100), ForeignKey("clinical_sample.id"), primary_key=True)
	value: Mapped[float] = mapped_column(Float(), nullable=False) # 0.0 - 1.0


# Gene bridge tables

class Antigen(Base):
	__tablename__ = "antigen"

	id: Mapped[str] = mapped_column(String(100), primary_key=True) # Ex. AGID00100
	catalogue_number: Mapped[int] = mapped_column(Integer, nullable=False)
	protein_target: Mapped[str] = mapped_column(String(100), nullable=False) # Used to map via uniprot to ENSEMBL, Ex. 1433BETA --> uniprot --> ENSG00000166913
	protein_target_gene_id: Mapped[str] = mapped_column(String(100), ForeignKey("gene.id"), nullable=False) # ENSEMBL ID of the protein target

# Gene table
class gene(Base):
	__tablename__ = "gene"

	id: Mapped[str] = mapped_column(String(100), primary_key=True) # ENSEMBL ID
	name: Mapped[str] = mapped_column(String(100), nullable=False)




"""
    Pre clinical tables: tables related to the storage/organization of data from in vitro datasets
"""


class PreClinicalDataset(Base):
    __tablename__ = "pre_clinical_dataset"

    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    name: Mapped[str] = mapped_column(Text())
