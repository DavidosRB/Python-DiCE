from typing import Optional, Any
import pandas as pd
import os 
import requests

# TODO: Implement Data Normalization function (normalize gene expression using z-scores) - not done here since my data is already normalised

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




def get_gene_id_mappings(force_download: bool = False)-> dict[str, str]:
    url = "https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_pub_eg_id&col=gd_pub_ensembl_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_hgnc_id&format=text&submit=submit"
    gene_id_mappings_file_path = os.path.join("..", "data", "gene_id_mappings.txt")

    if not os.path.exists(gene_id_mappings_file_path) or force_download:
        # Get the file from the genenames custom download
        response = requests.get(url)
        # Save the response text to a local .txt file
        with open(gene_id_mappings_file_path, "w", encoding="utf-8") as file:
            file.write(response.text)
        print(f"File saved to {gene_id_mappings_file_path}")
    else:
        print(f"File {gene_id_mappings_file_path} already exists. If you want to download it again, set force_download to True.")
    
    # Now, we turn this downloaded file into a dataframe and return the mapping from each of the columns to HGNC ID as a dictionary
    # Note that the NCBI ID will be cast into a string to prevent it from becoming a float 
    gene_id_df = pd.read_csv(os.path.join("..", "data", "gene_id_mappings.txt"), sep="\t", dtype={"NCBI Gene ID": str})

    gene_id_dict = {}

    cols = list(gene_id_df.columns)

    for _, row in gene_id_df.iterrows():
        for col in cols:
            col_val = row[col]
            if not pd.isna(row[col]):
                gene_id_dict[col_val] = row["HGNC ID"]

    return gene_id_dict

def convert_series_to_hgnc(series: pd.Series|Any, gene_id_mappings: dict) -> pd.Series:
    if not isinstance(series, pd.Series):
        try:
            series = pd.Series(series)
        except:
            raise Exception("Couldn't convert input object to a pandas Series object.")
    

    series = series.map(lambda x: gene_id_mappings.get(str(x)))

    return series