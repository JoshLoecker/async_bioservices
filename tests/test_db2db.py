import pytest
import pandas as pd

from async_bioservices import db2db, InputDatabase, OutputDatabase, TaxonID


# Cache tests
# 1: This uses async_cache, but it will be a cache miss. Retrieve results from server and cache them
# 2: This uses async_cache, but it will be a cache hit. Retrieve results from cache
# 3: This uses biodbnet_cache, but it will be a cache miss. Retrieve results from server and cache them
# 4: This uses biodbnet_cache, but it will be a cache hit. Retrieve results from cache
# 5: This does not use any cache. Retrieve results from server and do not cache them
@pytest.mark.parametrize("async_cache, biodbnet_cache",
                         [(True, False), (True, False), (False, True), (False, True), (False, False)])
@pytest.mark.parametrize("taxon_id", [TaxonID.HOMO_SAPIENS, TaxonID.MUS_MUSCULUS])
@pytest.mark.asyncio
async def test_fetch_gene_info(async_cache, biodbnet_cache, taxon_id):
    if taxon_id == TaxonID.HOMO_SAPIENS:
        input_values = ["1", "2", "3", "4"]
    elif taxon_id == TaxonID.MUS_MUSCULUS:
        input_values = ["14910", "22059", "11816", "21898"]
    
    result: pd.DataFrame = await db2db(
        input_values=input_values,
        input_db=InputDatabase.GENE_ID,
        output_db=OutputDatabase.GENE_SYMBOL,
        taxon_id=taxon_id,
        async_cache=async_cache,
        biodbnet_cache=biodbnet_cache,
        quiet=True
    )
    
    print(result)
    
    assert isinstance(result, pd.DataFrame)
    assert len(result) == len(input_values)
    assert 'Gene Symbol' in result.columns
    assert 'Gene ID' in result.columns
    
    if taxon_id == TaxonID.HOMO_SAPIENS:
        assert result[result["Gene ID"] == "1"]["Gene Symbol"].values[0] == "A1BG"
        assert result[result["Gene ID"] == "2"]["Gene Symbol"].values[0] == "A2M"
        assert result[result["Gene ID"] == "3"]["Gene Symbol"].values[0] == "A2MP1"
        assert result[result["Gene ID"] == "4"]["Gene Symbol"].values[0] == "-"
    elif taxon_id == TaxonID.MUS_MUSCULUS:
        assert result[result["Gene ID"] == "14910"]["Gene Symbol"].values[0] == "Gt(ROSA)26Sor"
        assert result[result["Gene ID"] == "22059"]["Gene Symbol"].values[0] == "Trp53"
        assert result[result["Gene ID"] == "11816"]["Gene Symbol"].values[0] == "Apoe"
        assert result[result["Gene ID"] == "21898"]["Gene Symbol"].values[0] == "Tlr4"


@pytest.mark.asyncio
async def test_fetch_gene_info_invalid_input_db_in_output_db():
    # Cannot have the same input_db and output_db values
    with pytest.raises(ValueError):
        await db2db(
            input_values=["123"],
            input_db=InputDatabase.GENE_ID,
            output_db=OutputDatabase.GENE_ID,
            taxon_id=TaxonID.HOMO_SAPIENS,
            async_cache=True,
            quiet=True
        )


@pytest.mark.asyncio
async def test_fetch_gene_info_invalid_concurrency():
    # Cannot have more than 25 concurrent workers
    with pytest.raises(ValueError):
        await db2db(
            input_values=["123"],
            input_db=InputDatabase.GENE_ID,
            output_db=OutputDatabase.GENE_SYMBOL,
            quiet=True,
            concurrency=25
        )


@pytest.mark.asyncio
async def test_fetch_gene_info_invalid_batch_length():
    # Cannot have more than 500 items per batch length
    with pytest.raises(ValueError):
        await db2db(
            input_values=["123"],
            input_db=InputDatabase.GENE_ID,
            output_db=OutputDatabase.GENE_SYMBOL,
            quiet=True,
            batch_length=600
        )
