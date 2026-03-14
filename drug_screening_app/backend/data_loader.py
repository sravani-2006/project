import pandas as pd


REQUIRED_COLUMNS = {"ligand_id", "smiles"}


def load_ligand_dataset(csv_path: str) -> pd.DataFrame:
    """Load ligand dataset and validate required columns."""
    df = pd.read_csv(csv_path)
    missing = REQUIRED_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(f"Ligand CSV missing required columns: {sorted(missing)}")
    return df[["ligand_id", "smiles"]].dropna().reset_index(drop=True)


def sample_demo_dataset(df: pd.DataFrame, sample_size: int = 100, random_state: int = 42) -> pd.DataFrame:
    """Random sample for lightweight demo runs."""
    n = min(sample_size, len(df))
    return df.sample(n=n, random_state=random_state).reset_index(drop=True)


def select_top_molecules_by_score(docking_df: pd.DataFrame, top_n: int = 5000) -> pd.DataFrame:
    """Sort by docking score ascending and keep top molecules."""
    return docking_df.sort_values("docking_score", ascending=True).head(top_n).reset_index(drop=True)
