from sqlalchemy import create_engine
import os
from datetime import datetime
from dotenv import load_dotenv
from urllib.parse import quote_plus
import pandas as pd

load_dotenv(override=True)


def alchemy_engine():
    password_cleaned = quote_plus(os.getenv("DATABASE_PASS"))
    engine = create_engine(
        f"mysql+pymysql://{os.getenv('DATABASE_USER')}:{password_cleaned}"
        f"@{os.getenv('DATABASE_IP')}:{os.getenv('PORT')}/{os.getenv('SELECTED_DB')}",
        echo=True,
    )
    return engine
