from async_bioservices.taxon_ids import TaxonIDs


def test_taxon_ids_values():
    assert TaxonIDs.HOMO_SAPIENS.value == 9606
    assert TaxonIDs.MUS_MUSCULUS.value == 10090


def test_taxon_ids_names():
    assert TaxonIDs.HOMO_SAPIENS.name == 'HOMO_SAPIENS'
    assert TaxonIDs.MUS_MUSCULUS.name == 'MUS_MUSCULUS'


def test_taxon_ids_type():
    assert isinstance(TaxonIDs.HOMO_SAPIENS, TaxonIDs)
    assert isinstance(TaxonIDs.MUS_MUSCULUS, TaxonIDs)
