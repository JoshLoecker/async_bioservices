from async_bioservices.input_database import InputDatabase


def test_name():
    assert InputDatabase.AFFY_GENECHIP_ARRAY.name == "AFFY_GENECHIP_ARRAY"
    assert InputDatabase.AFFY_ID.name == "AFFY_ID"
    assert InputDatabase.AFFY_TRANSCRIPT_CLUSTER_ID.name == "AFFY_TRANSCRIPT_CLUSTER_ID"
    assert InputDatabase.AGILENT_ID.name == "AGILENT_ID"
    assert InputDatabase.BIOCARTA_PATHWAY_NAME.name == "BIOCARTA_PATHWAY_NAME"
    assert InputDatabase.CODELINK_ID.name == "CODELINK_ID"
    assert InputDatabase.DBSNP_ID.name == "DBSNP_ID"
    assert InputDatabase.DRUGBANK_DRUG_ID.name == "DRUGBANK_DRUG_ID"
    assert InputDatabase.DRUGBANK_DRUG_NAME.name == "DRUGBANK_DRUG_NAME"
    assert InputDatabase.EC_NUMBER.name == "EC_NUMBER"
    assert InputDatabase.ENSEMBL_GENE_ID.name == "ENSEMBL_GENE_ID"
    assert InputDatabase.ENSEMBL_PROTEIN_ID.name == "ENSEMBL_PROTEIN_ID"
    assert InputDatabase.ENSEMBL_TRANSCRIPT_ID.name == "ENSEMBL_TRANSCRIPT_ID"
    assert InputDatabase.EST_ACCESSION.name == "EST_ACCESSION"
    assert InputDatabase.FLYBASE_GENE_ID.name == "FLYBASE_GENE_ID"
    assert InputDatabase.GENBANK_NUCLEOTIDE_ACCESSION.name == "GENBANK_NUCLEOTIDE_ACCESSION"
    assert InputDatabase.GENBANK_PROTEIN_ACCESSION.name == "GENBANK_PROTEIN_ACCESSION"
    assert InputDatabase.GENE_ID.name == "GENE_ID"
    assert InputDatabase.GENE_SYMBOL.name == "GENE_SYMBOL"
    assert InputDatabase.GENE_SYMBOL_AND_SYNONYMS.name == "GENE_SYMBOL_AND_SYNONYMS"
    assert InputDatabase.GENE_SYMBOL_ORDERED_LOCUS.name == "GENE_SYMBOL_ORDERED_LOCUS"
    assert InputDatabase.GENE_SYMBOL_ORF.name == "GENE_SYMBOL_ORF"
    assert InputDatabase.GI_NUMBER.name == "GI_NUMBER"
    assert InputDatabase.GO_ID.name == "GO_ID"
    assert InputDatabase.GSEA_STANDARD_NAME.name == "GSEA_STANDARD_NAME"
    assert InputDatabase.H_INV_LOCUS_ID.name == "H_INV_LOCUS_ID"
    assert InputDatabase.H_INV_PROTEIN_ID.name == "H_INV_PROTEIN_ID"
    assert InputDatabase.H_INV_TRANSCRIPT_ID.name == "H_INV_TRANSCRIPT_ID"
    assert InputDatabase.HGNC_ID.name == "HGNC_ID"
    assert InputDatabase.HMDB_METABOLITE.name == "HMDB_METABOLITE"
    assert InputDatabase.HOMOLOGENE_ID.name == "HOMOLOGENE_ID"
    assert InputDatabase.ILLUMINA_ID.name == "ILLUMINA_ID"
    assert InputDatabase.INTERPRO_ID.name == "INTERPRO_ID"
    assert InputDatabase.IPI_ID.name == "IPI_ID"
    assert InputDatabase.KEGG_COMPOUND_ID.name == "KEGG_COMPOUND_ID"
    assert InputDatabase.KEGG_COMPOUND_NAME.name == "KEGG_COMPOUND_NAME"
    assert InputDatabase.KEGG_DISEASE_ID.name == "KEGG_DISEASE_ID"
    assert InputDatabase.KEGG_DRUG_ID.name == "KEGG_DRUG_ID"
    assert InputDatabase.KEGG_DRUG_NAME.name == "KEGG_DRUG_NAME"
    assert InputDatabase.KEGG_GENE_ID.name == "KEGG_GENE_ID"
    assert InputDatabase.KEGG_PATHWAY_ID.name == "KEGG_PATHWAY_ID"
    assert InputDatabase.MAIZEGDB_ID.name == "MAIZEGDB_ID"
    assert InputDatabase.MGI_ID.name == "MGI_ID"
    assert InputDatabase.MIM_ID.name == "MIM_ID"
    assert InputDatabase.MIRBASE_ID.name == "MIRBASE_ID"
    assert InputDatabase.MIRBASE_MATURE_MIRNA_ACC.name == "MIRBASE_MATURE_MIRNA_ACC"
    assert InputDatabase.NCIPID_PATHWAY_NAME.name == "NCIPID_PATHWAY_NAME"
    assert InputDatabase.ORGANISM_SCIENTIFIC_NAME.name == "ORGANISM_SCIENTIFIC_NAME"
    assert InputDatabase.PDB_ID.name == "PDB_ID"
    assert InputDatabase.PFAM_ID.name == "PFAM_ID"
    assert InputDatabase.PHARMGKB_ID.name == "PHARMGKB_ID"
    assert InputDatabase.PUBCHEM_ID.name == "PUBCHEM_ID"
    assert InputDatabase.REACTOME_PATHWAY_NAME.name == "REACTOME_PATHWAY_NAME"
    assert InputDatabase.REFSEQ_GENOMIC_ACCESSION.name == "REFSEQ_GENOMIC_ACCESSION"
    assert InputDatabase.REFSEQ_MRNA_ACCESSION.name == "REFSEQ_MRNA_ACCESSION"
    assert InputDatabase.REFSEQ_PROTEIN_ACCESSION.name == "REFSEQ_PROTEIN_ACCESSION"
    assert InputDatabase.SGD_ID.name == "SGD_ID"
    assert InputDatabase.TAIR_ID.name == "TAIR_ID"
    assert InputDatabase.TAXON_ID.name == "TAXON_ID"
    assert InputDatabase.UNIGENE_ID.name == "UNIGENE_ID"
    assert InputDatabase.UNIPROT_ACCESSION.name == "UNIPROT_ACCESSION"
    assert InputDatabase.UNIPROT_ENTRY_NAME.name == "UNIPROT_ENTRY_NAME"
    assert InputDatabase.UNIPROT_PROTEIN_NAME.name == "UNIPROT_PROTEIN_NAME"
    assert InputDatabase.UNISTS_ID.name == "UNISTS_ID"


def test_type():
    assert isinstance(InputDatabase.AFFY_GENECHIP_ARRAY, InputDatabase)
    assert isinstance(InputDatabase.AFFY_ID, InputDatabase)
    assert isinstance(InputDatabase.AFFY_TRANSCRIPT_CLUSTER_ID, InputDatabase)
    assert isinstance(InputDatabase.AGILENT_ID, InputDatabase)
    assert isinstance(InputDatabase.BIOCARTA_PATHWAY_NAME, InputDatabase)
    assert isinstance(InputDatabase.CODELINK_ID, InputDatabase)
    assert isinstance(InputDatabase.DBSNP_ID, InputDatabase)
    assert isinstance(InputDatabase.DRUGBANK_DRUG_ID, InputDatabase)
    assert isinstance(InputDatabase.DRUGBANK_DRUG_NAME, InputDatabase)
    assert isinstance(InputDatabase.EC_NUMBER, InputDatabase)
    assert isinstance(InputDatabase.ENSEMBL_GENE_ID, InputDatabase)
    assert isinstance(InputDatabase.ENSEMBL_PROTEIN_ID, InputDatabase)
    assert isinstance(InputDatabase.ENSEMBL_TRANSCRIPT_ID, InputDatabase)
    assert isinstance(InputDatabase.EST_ACCESSION, InputDatabase)
    assert isinstance(InputDatabase.FLYBASE_GENE_ID, InputDatabase)
    assert isinstance(InputDatabase.GENBANK_NUCLEOTIDE_ACCESSION, InputDatabase)
    assert isinstance(InputDatabase.GENBANK_PROTEIN_ACCESSION, InputDatabase)
    assert isinstance(InputDatabase.GENE_ID, InputDatabase)
    assert isinstance(InputDatabase.GENE_SYMBOL, InputDatabase)
    assert isinstance(InputDatabase.GENE_SYMBOL_AND_SYNONYMS, InputDatabase)
    assert isinstance(InputDatabase.GENE_SYMBOL_ORDERED_LOCUS, InputDatabase)
    assert isinstance(InputDatabase.GENE_SYMBOL_ORF, InputDatabase)
    assert isinstance(InputDatabase.GI_NUMBER, InputDatabase)
    assert isinstance(InputDatabase.GO_ID, InputDatabase)
    assert isinstance(InputDatabase.GSEA_STANDARD_NAME, InputDatabase)
    assert isinstance(InputDatabase.H_INV_LOCUS_ID, InputDatabase)
    assert isinstance(InputDatabase.H_INV_PROTEIN_ID, InputDatabase)
    assert isinstance(InputDatabase.H_INV_TRANSCRIPT_ID, InputDatabase)
    assert isinstance(InputDatabase.HGNC_ID, InputDatabase)
    assert isinstance(InputDatabase.HMDB_METABOLITE, InputDatabase)
    assert isinstance(InputDatabase.HOMOLOGENE_ID, InputDatabase)
    assert isinstance(InputDatabase.ILLUMINA_ID, InputDatabase)
    assert isinstance(InputDatabase.INTERPRO_ID, InputDatabase)
    assert isinstance(InputDatabase.IPI_ID, InputDatabase)
    assert isinstance(InputDatabase.KEGG_COMPOUND_ID, InputDatabase)
    assert isinstance(InputDatabase.KEGG_COMPOUND_NAME, InputDatabase)
    assert isinstance(InputDatabase.KEGG_DISEASE_ID, InputDatabase)
    assert isinstance(InputDatabase.KEGG_DRUG_ID, InputDatabase)
    assert isinstance(InputDatabase.KEGG_DRUG_NAME, InputDatabase)
    assert isinstance(InputDatabase.KEGG_GENE_ID, InputDatabase)
    assert isinstance(InputDatabase.KEGG_PATHWAY_ID, InputDatabase)
    assert isinstance(InputDatabase.MAIZEGDB_ID, InputDatabase)
    assert isinstance(InputDatabase.MGI_ID, InputDatabase)
    assert isinstance(InputDatabase.MIM_ID, InputDatabase)
    assert isinstance(InputDatabase.MIRBASE_ID, InputDatabase)
    assert isinstance(InputDatabase.MIRBASE_MATURE_MIRNA_ACC, InputDatabase)
    assert isinstance(InputDatabase.NCIPID_PATHWAY_NAME, InputDatabase)
    assert isinstance(InputDatabase.ORGANISM_SCIENTIFIC_NAME, InputDatabase)
    assert isinstance(InputDatabase.PDB_ID, InputDatabase)
    assert isinstance(InputDatabase.PFAM_ID, InputDatabase)
    assert isinstance(InputDatabase.PHARMGKB_ID, InputDatabase)
    assert isinstance(InputDatabase.PUBCHEM_ID, InputDatabase)
    assert isinstance(InputDatabase.REACTOME_PATHWAY_NAME, InputDatabase)
    assert isinstance(InputDatabase.REFSEQ_GENOMIC_ACCESSION, InputDatabase)
    assert isinstance(InputDatabase.REFSEQ_MRNA_ACCESSION, InputDatabase)
    assert isinstance(InputDatabase.REFSEQ_PROTEIN_ACCESSION, InputDatabase)
    assert isinstance(InputDatabase.SGD_ID, InputDatabase)
    assert isinstance(InputDatabase.TAIR_ID, InputDatabase)
    assert isinstance(InputDatabase.TAXON_ID, InputDatabase)
    assert isinstance(InputDatabase.UNIGENE_ID, InputDatabase)
    assert isinstance(InputDatabase.UNIPROT_ACCESSION, InputDatabase)
    assert isinstance(InputDatabase.UNIPROT_ENTRY_NAME, InputDatabase)
    assert isinstance(InputDatabase.UNIPROT_PROTEIN_NAME, InputDatabase)
    assert isinstance(InputDatabase.UNISTS_ID, InputDatabase)


def test_value():
    assert InputDatabase.AFFY_GENECHIP_ARRAY.value == "Affy GeneChip Array"
    assert InputDatabase.AFFY_ID.value == "Affy ID"
    assert InputDatabase.AFFY_TRANSCRIPT_CLUSTER_ID.value == "Affy Transcript Cluster ID"
    assert InputDatabase.AGILENT_ID.value == "Agilent ID"
    assert InputDatabase.BIOCARTA_PATHWAY_NAME.value == "Biocarta Pathway Name"
    assert InputDatabase.CODELINK_ID.value == "CodeLink ID"
    assert InputDatabase.DBSNP_ID.value == "dbSNP ID"
    assert InputDatabase.DRUGBANK_DRUG_ID.value == "DrugBank Drug ID"
    assert InputDatabase.DRUGBANK_DRUG_NAME.value == "DrugBank Drug Name"
    assert InputDatabase.EC_NUMBER.value == "EC Number"
    assert InputDatabase.ENSEMBL_GENE_ID.value == "Ensembl Gene ID"
    assert InputDatabase.ENSEMBL_PROTEIN_ID.value == "Ensembl Protein ID"
    assert InputDatabase.ENSEMBL_TRANSCRIPT_ID.value == "Ensembl Transcript ID"
    assert InputDatabase.EST_ACCESSION.value == "EST Accession"
    assert InputDatabase.FLYBASE_GENE_ID.value == "FlyBase Gene ID"
    assert InputDatabase.GENBANK_NUCLEOTIDE_ACCESSION.value == "GenBank Nucleotide Accession"
    assert InputDatabase.GENBANK_PROTEIN_ACCESSION.value == "GenBank Protein Accession"
    assert InputDatabase.GENE_ID.value == "Gene ID"
    assert InputDatabase.GENE_SYMBOL.value == "Gene Symbol"
    assert InputDatabase.GENE_SYMBOL_AND_SYNONYMS.value == "Gene Symbol and Synonyms"
    assert InputDatabase.GENE_SYMBOL_ORDERED_LOCUS.value == "Gene Symbol Ordered Locus"
    assert InputDatabase.GENE_SYMBOL_ORF.value == "Gene Symbol ORF"
    assert InputDatabase.GI_NUMBER.value == "GI Number"
    assert InputDatabase.GO_ID.value == "GO ID"
    assert InputDatabase.GSEA_STANDARD_NAME.value == "GSEA Standard Name"
    assert InputDatabase.H_INV_LOCUS_ID.value == "H-Inv Locus ID"
    assert InputDatabase.H_INV_PROTEIN_ID.value == "H-Inv Protein ID"
    assert InputDatabase.H_INV_TRANSCRIPT_ID.value == "H-Inv Transcript ID"
    assert InputDatabase.HGNC_ID.value == "HGNC ID"
    assert InputDatabase.HMDB_METABOLITE.value == "HMDB Metabolite"
    assert InputDatabase.HOMOLOGENE_ID.value == "HomoloGene ID"
    assert InputDatabase.ILLUMINA_ID.value == "Illumina ID"
    assert InputDatabase.INTERPRO_ID.value == "InterPro ID"
    assert InputDatabase.IPI_ID.value == "IPI ID"
    assert InputDatabase.KEGG_COMPOUND_ID.value == "KEGG Compound ID"
    assert InputDatabase.KEGG_COMPOUND_NAME.value == "KEGG Compound Name"
    assert InputDatabase.KEGG_DISEASE_ID.value == "KEGG Disease ID"
    assert InputDatabase.KEGG_DRUG_ID.value == "KEGG Drug ID"
    assert InputDatabase.KEGG_DRUG_NAME.value == "KEGG Drug Name"
    assert InputDatabase.KEGG_GENE_ID.value == "KEGG Gene ID"
    assert InputDatabase.KEGG_PATHWAY_ID.value == "KEGG Pathway ID"
    assert InputDatabase.MAIZEGDB_ID.value == "MaizeGDB ID"
    assert InputDatabase.MGI_ID.value == "MGI ID"
    assert InputDatabase.MIM_ID.value == "MIM ID"
    assert InputDatabase.MIRBASE_ID.value == "miRBase ID"
    assert InputDatabase.MIRBASE_MATURE_MIRNA_ACC.value == "miRBase Mature miRNA Acc"
    assert InputDatabase.NCIPID_PATHWAY_NAME.value == "NCIPID Pathway Name"
    assert InputDatabase.ORGANISM_SCIENTIFIC_NAME.value == "Organism Scientific Name"
    assert InputDatabase.PDB_ID.value == "PDB ID"
    assert InputDatabase.PFAM_ID.value == "Pfam ID"
    assert InputDatabase.PHARMGKB_ID.value == "PharmGKB ID"
    assert InputDatabase.PUBCHEM_ID.value == "PubChem ID"
    assert InputDatabase.REACTOME_PATHWAY_NAME.value == "Reactome Pathway Name"
    assert InputDatabase.REFSEQ_GENOMIC_ACCESSION.value == "RefSeq Genomic Accession"
    assert InputDatabase.REFSEQ_MRNA_ACCESSION.value == "RefSeq mRNA Accession"
    assert InputDatabase.REFSEQ_PROTEIN_ACCESSION.value == "RefSeq Protein Accession"
    assert InputDatabase.SGD_ID.value == "SGD ID"
    assert InputDatabase.TAIR_ID.value == "TAIR ID"
    assert InputDatabase.TAXON_ID.value == "Taxon ID"
    assert InputDatabase.UNIGENE_ID.value == "UniGene ID"
    assert InputDatabase.UNIPROT_ACCESSION.value == "UniProt Accession"
    assert InputDatabase.UNIPROT_ENTRY_NAME.value == "UniProt Entry Name"
    assert InputDatabase.UNIPROT_PROTEIN_NAME.value == "UniProt Protein Name"
    assert InputDatabase.UNISTS_ID.value == "UniSTS ID"
