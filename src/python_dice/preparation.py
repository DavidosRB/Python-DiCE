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
    url = "https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_pub_eg_id&col=gd_pub_ensembl_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_hgnc_id&format=text&submit=submit"
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
    
    # Now, we turn this downloaded file into a dataframe
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

def filter_candidate_genes(raw_df: pd.DataFrame, dea_df: pd.DataFrame, gene_id_column: Optional[str] = None, candidate_gene_column: str = "GeneID") -> pd.DataFrame:
    """Filter the raw DataFrame and keep only genes present in the candidate genes

    Parameters
    ----------
    raw_df : pd.DataFrame
        The raw gene expression data to filter using candidate genes
    gene_id_column : str
        The column that contains the gene IDs in the raw data
    dea_df: pd.DataFrame 
        The filtered DEA data that contains only the candidate genes
    candidate_gene_column: str
        The column where to find the gene IDs in the DEA filtered data. From GEO2R results, this should usually be "GeneID".

    Returns
    -------
    pd.DataFrame
        The filtered DataFrame that only contains candidate genes
    """
    # Turn candidate genes to a set for easier comparison
    set_candidate_genes = set(dea_df[candidate_gene_column])

    if gene_id_column is None:
        filtered_raw_df = raw_df.loc[raw_df.index.isin(set_candidate_genes)]
    else:
        filtered_raw_df = raw_df.loc[raw_df[gene_id_column].isin(set_candidate_genes)]


    return filtered_raw_df

def get_string_db(force_download: bool = False)-> tuple[pd.DataFrame, pd.DataFrame]:
    """Get the two necessary files from the StringDB and return them as Pandas DataFrames

    Parameters
    ----------
    force_download : bool, optional
        Whether or not to force the download of the files even if they already exist, by default False

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame]
        The two Dataframes we extract from the StringDB - the links between proteins and the (additional) information we need
    """
    links_url = "https://stringdb-downloads.org/download/stream/protein.links.v12.0/9606.protein.links.v12.0.min400.onlyAB.txt.gz"
    info_url = "https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz"
    

    protein_links_file_path = os.path.join("..", "data", "9606.protein.links.v12.0.min400.onlyAB.txt.gz")
    protein_info_file_path = os.path.join("..", "data", "9606.protein.info.v12.0.txt.gz")

    # First, download protein links
    if not os.path.exists(protein_links_file_path) or force_download:
        # Get the file from the stringdb downloads site
        response = requests.get(links_url)
        response.raw.decode_content = False  # Prevent auto-decompression
        # Save the response text to a local .txt file
        with open(protein_links_file_path, "wb") as file: # Use "wb" to write the content as bytes since it's a gzipped file
            file.write(response.content) # Also write content instead of text since it's a gzipped file
        print(f"File saved to {protein_links_file_path}")
    else:
        print(f"File {protein_links_file_path} already exists. If you want to download it again, set force_download to True.")
    
    # Then, download protein info
    if not os.path.exists(protein_info_file_path) or force_download:
        # Get the file from the stringdb downloads site
        response = requests.get(info_url)
        response.raw.decode_content = False  # Prevent auto-decompression
        # Save the response text to a local .txt file
        with open(protein_info_file_path, "wb") as file:
            file.write(response.content)
        print(f"File saved to {protein_info_file_path}")
    else:
        print(f"File {protein_info_file_path} already exists. If you want to download it again, set force_download to True.")


    # We can then turn these files into dataframes using Pandas and return them

    protein_links_df = pd.read_csv(protein_links_file_path, sep=" ", compression="gzip") # Note that the file is gzipped, so we need to specify this in the read_csv function
    protein_info_df = pd.read_csv(protein_info_file_path, sep="\t", compression="gzip")

    return protein_links_df, protein_info_df


def map_protein_links_symbols(protein_links_df: pd.DataFrame, protein_info_df: pd.DataFrame) -> pd.DataFrame:
    """Map the Ensembl Protein IDs from the StringDB's Protein Links file/table to use HGNC symbols instead

    Parameters
    ----------
    protein_links_df : pd.DataFrame
        The StringDB's links and relationships between proteins, in a Pandas DataFrame
    protein_info_df : pd.DataFrame
        The StringDB's (additional) information for every protein, including the preferred HGNC symbol for eveery Ensembl Protein ID, in a Pandas DataFrame

    Returns
    -------
    pd.DataFrame
        The original protein_links_df, but now with HGNC symbols instead of Ensembl Protein IDs in the 2 protein columns
    """

    # Get the mapping from Ensemble Protein ID to HGNC Symbol using the protein Info DataFrame
    # Zipping the two columns and turning them into a dictionary turns them into a dictionary with ENSP IDs as Keys and HGNC symbols as value
    ensp_symbol_mapping_dict = dict(zip(protein_info_df["#string_protein_id"], protein_info_df["preferred_name"]))
    # Using this dictionary, we can then map the column of the protein links DataFrame
    mapped_protein_links_df = protein_links_df.copy()
    # Map both relevant columns to use HGNC Symbols instead of ENSP IDs
    mapped_protein_links_df["protein1"] = mapped_protein_links_df["protein1"].map(ensp_symbol_mapping_dict)
    mapped_protein_links_df["protein2"] = mapped_protein_links_df["protein2"].map(ensp_symbol_mapping_dict)

    return mapped_protein_links_df

