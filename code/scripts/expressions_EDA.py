import io

import boto3
import pandas as pd


def load_expressions_from_s3(
    samples,
    output_path="expression_matrix.csv",
    bucket="bostongene-deconv-oss-rnaseq",
    verbose=True,
):
    """
    Download expression matrices from S3 for a given list of samples.

    Args:
        samples (list of str): List of sample IDs.
        output_path (str): Path to save the final concatenated expression matrix.
        bucket (str): Name of the S3 bucket.
        verbose (bool): If True, prints status messages.

    Returns:
        pd.DataFrame: Combined expression matrix (genes x samples).
    """
    s3 = boto3.client("s3")
    expr_list = []

    for sample in samples:
        df = None
        try:
            key = f"database/{sample}/kallisto/expressions.tsv"
            response = s3.get_object(Bucket=bucket, Key=key)
            df = pd.read_csv(
                io.BytesIO(response["Body"].read()), sep="\t", index_col=0
            )
            if verbose:
                print(f"Loaded: {key}")
        except s3.exceptions.NoSuchKey:
            try:
                key = f"database/{sample}/kallisto/{sample}-kallisto-Xena-gene-TPM_without_noncoding.tsv"
                response = s3.get_object(Bucket=bucket, Key=key)
                df = pd.read_csv(
                    io.BytesIO(response["Body"].read()), sep="\t", index_col=0
                )
                if verbose:
                    print(f"Loaded (old format): {key}")
            except s3.exceptions.NoSuchKey:
                if verbose:
                    print(f"No file found for {sample}")

        if df is not None:
            df.columns = [sample]
            expr_list.append(df)

    if expr_list:
        expr = pd.concat(expr_list, axis=1)
        expr.to_csv(output_path)
        if verbose:
            print(
                f"\n Final expression matrix saved to '{output_path}', shape: {expr.shape}"
            )
        return expr
    else:
        if verbose:
            print("No expression files were successfully loaded.")
        return pd.DataFrame()


def load_est_counts_from_s3(
    samples, bucket="bostongene-deconv-oss-rnaseq", verbose=True
) -> pd.DataFrame:
    """
    Download estimated count matrices (without noncoding genes) from S3 for given samples.

    Parameters:
        samples (list of str): List of sample IDs.
        bucket (str): S3 bucket name (default: 'bostongene-deconv-oss-rnaseq').
        verbose (bool): If True, prints status messages.

    Returns:
        pd.DataFrame: Combined expression matrix (genes x samples).
    """
    s3 = boto3.client("s3")
    expr_list = []

    for sample in samples:
        key = f"database/{sample}/kallisto/{sample}-kallisto-Xena-gene-est_counts_without_noncoding.tsv"

        try:
            response = s3.get_object(Bucket=bucket, Key=key)
            df = pd.read_csv(
                io.BytesIO(response["Body"].read()), sep="\t", index_col=0
            )
            df.columns = [sample]
            expr_list.append(df)
            if verbose:
                print(f"Loaded: {key}")
        except s3.exceptions.NoSuchKey:
            if verbose:
                print(f"File not found: {key}")
            continue

    if expr_list:
        expr_matrix = pd.concat(expr_list, axis=1)
        if verbose:
            print(f"\nFinal expression matrix shape: {expr_matrix.shape}")
        return expr_matrix
    else:
        if verbose:
            print("No expression files were successfully loaded.")
        return pd.DataFrame()


def renormalize_tpm(
    expr: pd.DataFrame, genes_in_expression: list[str]
) -> pd.DataFrame:
    """Renornalize expression values so that expressions of one sample sum up to 10^6.

    Args:
        expr (pd.DataFrame): Expression values: genes as index.
        genes_in_expression (List[str]): Genes to take into account.

    Returns:
        pd.DataFrame: Renormalized expression values with genes_in_expression as index
        and the same columns.
    """

    expr = expr.loc[list(genes_in_expression)]
    expr = (expr / expr.sum()) * 10**6

    return expr
