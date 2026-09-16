from .seeding_coordinator_engine import alchemy_engine

from ..models.tables import (
    Base,
    Dataset,
    ClinicalSample,
    ClinicalAntigen,
    ClinicalProbe,
    ClinicalRNA,
    ClinicalMutation,
    ClinicalCNV,
    ClinicalRPPA,
    ClinicalMiRNA,
    ClinicalMethylation,
    ClinicalSlide,
    ClinicalTile,
    ClinicalEmbedding,
    PreClinicalGene,
)


def main():
    Base.metadata.create_all(
        bind=alchemy_engine(),
        tables=[
            Dataset.__table__,
            PreClinicalGene.__table__,
            ClinicalSample.__table__,
            ClinicalAntigen.__table__,
            ClinicalProbe.__table__,
            ClinicalRNA.__table__,
            ClinicalMutation.__table__,
            ClinicalCNV.__table__,
            ClinicalRPPA.__table__,
            ClinicalMiRNA.__table__,
            ClinicalMethylation.__table__,
            ClinicalSlide.__table__,
            ClinicalTile.__table__,
            ClinicalEmbedding.__table__,
        ],
    )


if __name__ == "__main__":
    main()
