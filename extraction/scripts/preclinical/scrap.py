import pandas as pd

for name, path in {
    "GCSI": "extraction/data/proc/preclinical/GCSI/pre_clinical_copy_number_variation.csv",
    "CCLE": "extraction/data/proc/preclinical/CCLE/pre_clinical_copy_number_variation.csv",
}.items():
    df = pd.read_csv(path)
    print("\n", name)
    print(df["value"].describe())
    print("median:", df["value"].median())
    print("negative rows:", (df["value"] < 0).sum())