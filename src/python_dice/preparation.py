from typing import Optional
import pandas as pd

def filter_dea_table(dea_df: pd.DataFrame, pval_threshold: float = 0.05, fc_threshold: Optional[float] = None) -> pd.DataFrame:
    """Filter a DEA table that is a result of a GEO2R analysis with the given p-value and Fold Change thresholds

    Parameters
    ----------
    dea_df : pd.DataFrame
        The resulting DEA table from the GEO2R analysis, in a Dataframe format
    pval_threshold : float
        The threshold for the adjusted p-value - will not keep genes with padj ABOVE this given threshold
    fc_threshold : float
        The threshold for the log2FoldChange - will not keep genes with log2FC BELOW this given threshold

    Returns
    -------
    pd.DataFrame
        the resulting filtered DataFrame
    """
    dea_df = dea_df[dea_df["padj"] < pval_threshold]
    if fc_threshold is not None:    
        dea_df = dea_df[dea_df["log2FoldChange"] > fc_threshold]

    return dea_df