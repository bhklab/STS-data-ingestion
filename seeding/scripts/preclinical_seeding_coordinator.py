from .seeding_coordinator_engine import alchemy_engine

from ..models.tables import (
    Base,
    Dataset,
    PreClinicalCellLine,
    PreClinicalTreatmentResponse,
    PreClinicalSample,
    PreClinicalRnaSeq,
    PreClinicalMutation,
    PreClinicalCopyNumberVariation,
    PreClinicalGene,
)


def main():
    Base.metadata.create_all(
        bind=alchemy_engine(),
        tables=[
            Dataset.__table__,
            PreClinicalCellLine.__table__,
            PreClinicalTreatmentResponse.__table__,
            PreClinicalSample.__table__,
            PreClinicalRnaSeq.__table__,
            PreClinicalMutation.__table__,
            PreClinicalCopyNumberVariation.__table__,
            PreClinicalGene.__table__,
        ],
    )


if __name__ == "__main__":
    main()