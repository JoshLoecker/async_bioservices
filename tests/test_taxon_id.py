from async_bioservices.taxon_id import TaxonID


def test_taxon_ids_values():
    assert TaxonID.HOMO_SAPIENS.value == 9606
    assert TaxonID.MUS_MUSCULUS.value == 10090


def test_taxon_ids_names():
    assert TaxonID.HOMO_SAPIENS.name == 'HOMO_SAPIENS'
    assert TaxonID.MUS_MUSCULUS.name == 'MUS_MUSCULUS'


def test_taxon_ids_type():
    assert isinstance(TaxonID.HOMO_SAPIENS, TaxonID)
    assert isinstance(TaxonID.MUS_MUSCULUS, TaxonID)
