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

    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    name: Mapped[str] = mapped_column(Text())
