import pytest
import pandas as pd

from async_bioservices import db2db, InputDatabase, OutputDatabase, TaxonID


@pytest.mark.parametrize("cache", [True, False])
@pytest.mark.asyncio
async def test_fetch_gene_info(cache):
    input_values = [str(i) for i in range(100)]
    input_db = InputDatabase.GENE_ID
    output_db = OutputDatabase.GENE_SYMBOL
    taxon_id = TaxonID.HOMO_SAPIENS
    quiet = True
    
    result: pd.DataFrame = await db2db(
        input_values=input_values,
        input_db=input_db,
        output_db=output_db,
        taxon_id=taxon_id,
        cache=cache,
        quiet=quiet
    )
    
    assert isinstance(result, pd.DataFrame)
    assert len(result) == len(input_values)
    
    assert 'Gene Symbol' in result.columns
    assert 'Gene ID' in result.columns


@pytest.mark.asyncio
async def test_fetch_gene_info_invalid_input_db_in_output_db():
    with pytest.raises(ValueError):
        await db2db(
            input_values=["123"],
            input_db=InputDatabase.GENE_ID,
            output_db=OutputDatabase.GENE_ID,
            taxon_id=TaxonID.HOMO_SAPIENS,
            cache=True,
            quiet=True
        )


@pytest.mark.asyncio
async def test_fetch_gene_info_invalid_concurrency():
    with pytest.raises(ValueError):
        await db2db(
            input_values=["123"],
            input_db=InputDatabase.GENE_ID,
            output_db=OutputDatabase.GENE_SYMBOL,
            taxon_id=TaxonID.HOMO_SAPIENS,
            cache=True,
            quiet=True,
            concurrency=25
        )


@pytest.mark.asyncio
async def test_fetch_gene_info_invalid_batch_length():
    with pytest.raises(ValueError):
        await db2db(
            input_values=["123"],
            input_db=InputDatabase.GENE_ID,
            output_db=OutputDatabase.GENE_SYMBOL,
            taxon_id=TaxonID.HOMO_SAPIENS,
            cache=True,
            quiet=True,
            batch_length=600
        )
