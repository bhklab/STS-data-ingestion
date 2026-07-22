from datetime import datetime
import pandas as pd
from .seeding_coordinator_engine import alchemy_engine

from ..models.tables import Base, ClinicalDataset, ClinicalPatient

Base.metadata.create_all(
    bind=alchemy_engine(), tables=[ClinicalPatient.__table__]
)  # Defined table creation
