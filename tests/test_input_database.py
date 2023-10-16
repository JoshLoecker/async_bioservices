from async_bioservices.input_database import InputDatabase


def test_name():
    assert InputDatabase.AFFY_GENECHIP_ARRAY.name == ""
    assert InputDatabase.AFFY_ID.name == ""
    assert InputDatabase.AFFY_TRANSCRIPT_CLUSTER_ID.name == "AFFY_GENECHIP_ARRAY"
    assert InputDatabase.AGILENT_ID.name == "AFFY_ID"
    assert InputDatabase.BIOCARTA_PATHWAY_NAME.name == "AFFY_TRANSCRIPT_CLUSTER_ID"
    assert InputDatabase.CODELINK_ID.name == "AGILENT_ID"
    assert InputDatabase.DBSNP_ID.name == "BIOCARTA_PATHWAY_NAME"
    assert InputDatabase.DRUGBANK_DRUG_ID.name == "CODELINK_ID"
    assert InputDatabase.DRUGBANK_DRUG_NAME.name == "DBSNP_ID"
    assert InputDatabase.EC_NUMBER.name == "DRUGBANK_DRUG_ID"
    assert InputDatabase.ENSEMBL_GENE_ID.name == "DRUGBANK_DRUG_NAME"
    assert InputDatabase.ENSEMBL_PROTEIN_ID.name == "EC_NUMBER"
    assert InputDatabase.ENSEMBL_TRANSCRIPT_ID.name == "ENSEMBL_GENE_ID"
    assert InputDatabase.EST_ACCESSION.name == "ENSEMBL_PROTEIN_ID"
    assert InputDatabase.FLYBASE_GENE_ID.name == "ENSEMBL_TRANSCRIPT_ID"
    assert InputDatabase.GENBANK_NUCLEOTIDE_ACCESSION.name == "EST_ACCESSION"
    assert InputDatabase.GENBANK_PROTEIN_ACCESSION.name == "FLYBASE_GENE_ID"
    assert InputDatabase.GENE_ID.name == "GENBANK_NUCLEOTIDE_ACCESSION"
    assert InputDatabase.GENE_SYMBOL.name == "GENBANK_PROTEIN_ACCESSION"
    assert InputDatabase.GENE_SYMBOL_AND_SYNONYMS.name == "GENE_ID"
    assert InputDatabase.GENE_SYMBOL_ORDERED_LOCUS.name == "GENE_SYMBOL"
    assert InputDatabase.GENE_SYMBOL_ORF.name == "GENE_SYMBOL_AND_SYNONYMS"
    assert InputDatabase.GI_NUMBER.name == "GENE_SYMBOL_ORDERED_LOCUS"
    assert InputDatabase.GO_ID.name == "GENE_SYMBOL_ORF"
    assert InputDatabase.GSEA_STANDARD_NAME.name == "GI_NUMBER"
    assert InputDatabase.H_INV_LOCUS_ID.name == "GO_ID"
    assert InputDatabase.H_INV_PROTEIN_ID.name == "GSEA_STANDARD_NAME"
    assert InputDatabase.H_INV_TRANSCRIPT_ID.name == "H_INV_LOCUS_ID"
    assert InputDatabase.HGNC_ID.name == "H_INV_PROTEIN_ID"
    assert InputDatabase.HMDB_METABOLITE.name == "H_INV_TRANSCRIPT_ID"
    assert InputDatabase.HOMOLOGENE_ID.name == "HGNC_ID"
    assert InputDatabase.ILLUMINA_ID.name == "HMDB_METABOLITE"
    assert InputDatabase.INTERPRO_ID.name == "HOMOLOGENE_ID"
    assert InputDatabase.IPI_ID.name == "ILLUMINA_ID"
    assert InputDatabase.KEGG_COMPOUND_ID.name == "INTERPRO_ID"
    assert InputDatabase.KEGG_COMPOUND_NAME.name == "IPI_ID"
    assert InputDatabase.KEGG_DISEASE_ID.name == "KEGG_COMPOUND_ID"
    assert InputDatabase.KEGG_DRUG_ID.name == "KEGG_COMPOUND_NAME"
    assert InputDatabase.KEGG_DRUG_NAME.name == "KEGG_DISEASE_ID"
    assert InputDatabase.KEGG_GENE_ID.name == "KEGG_DRUG_ID"
    assert InputDatabase.KEGG_PATHWAY_ID.name == "KEGG_DRUG_NAME"
    assert InputDatabase.MAIZEGDB_ID.name == "KEGG_GENE_ID"
    assert InputDatabase.MGI_ID.name == "KEGG_PATHWAY_ID"
    assert InputDatabase.MIM_ID.name == "MAIZEGDB_ID"
    assert InputDatabase.MIRBASE_ID.name == "MGI_ID"
    assert InputDatabase.MIRBASE_MATURE_MIRNA_ACC.name == "MIM_ID"
    assert InputDatabase.NCIPID_PATHWAY_NAME.name == "MIRBASE_ID"
    assert InputDatabase.ORGANISM_SCIENTIFIC_NAME.name == "MIRBASE_MATURE_MIRNA_ACC"
    assert InputDatabase.PDB_ID.name == "NCIPID_PATHWAY_NAME"
    assert InputDatabase.PFAM_ID.name == "ORGANISM_SCIENTIFIC_NAME"
    assert InputDatabase.PHARMGKB_ID.name == "PDB_ID"
    assert InputDatabase.PUBCHEM_ID.name == "PFAM_ID"
    assert InputDatabase.REACTOME_PATHWAY_NAME.name == "PHARMGKB_ID"
    assert InputDatabase.REFSEQ_GENOMIC_ACCESSION.name == "PUBCHEM_ID"
    assert InputDatabase.REFSEQ_MRNA_ACCESSION.name == "REACTOME_PATHWAY_NAME"
    assert InputDatabase.REFSEQ_PROTEIN_ACCESSION.name == "REFSEQ_GENOMIC_ACCESSION"
    assert InputDatabase.SGD_ID.name == "REFSEQ_MRNA_ACCESSION"
    assert InputDatabase.TAIR_ID.name == "REFSEQ_PROTEIN_ACCESSION"
    assert InputDatabase.TAXON_ID.name == "SGD_ID"
    assert InputDatabase.UNIGENE_ID.name == "TAIR_ID"
    assert InputDatabase.UNIPROT_ACCESSION.name == "TAXON_ID"
    assert InputDatabase.UNIPROT_ENTRY_NAME.name == "UNIGENE_ID"
    assert InputDatabase.UNIPROT_PROTEIN_NAME.name == "UNIPROT_ACCESSION"
    assert InputDatabase.UNISTS_ID.name == "UNIPROT_ENTRY_NAME"


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
    assert InputDatabase.AFFY_GENECHIP_ARRAY == "Affy GeneChip Array"
    assert InputDatabase.AFFY_ID == "Affy ID"
    assert InputDatabase.AFFY_TRANSCRIPT_CLUSTER_ID == "Affy Transcript Cluster ID"
    assert InputDatabase.AGILENT_ID == "Agilent ID"
    assert InputDatabase.BIOCARTA_PATHWAY_NAME == "Biocarta Pathway Name"
    assert InputDatabase.CODELINK_ID == "CodeLink ID"
    assert InputDatabase.DBSNP_ID == "dbSNP ID"
    assert InputDatabase.DRUGBANK_DRUG_ID == "DrugBank Drug ID"
    assert InputDatabase.DRUGBANK_DRUG_NAME == "DrugBank Drug Name"
    assert InputDatabase.EC_NUMBER == "EC Number"
    assert InputDatabase.ENSEMBL_GENE_ID == "Ensembl Gene ID"
    assert InputDatabase.ENSEMBL_PROTEIN_ID == "Ensembl Protein ID"
    assert InputDatabase.ENSEMBL_TRANSCRIPT_ID == "Ensembl Transcript ID"
    assert InputDatabase.EST_ACCESSION == "EST Accession"
    assert InputDatabase.FLYBASE_GENE_ID == "FlyBase Gene ID"
    assert InputDatabase.GENBANK_NUCLEOTIDE_ACCESSION == "GenBank Nucleotide Accession"
    assert InputDatabase.GENBANK_PROTEIN_ACCESSION == "GenBank Protein Accession"
    assert InputDatabase.GENE_ID == "Gene ID"
    assert InputDatabase.GENE_SYMBOL == "Gene Symbol"
    assert InputDatabase.GENE_SYMBOL_AND_SYNONYMS == "Gene Symbol and Synonyms"
    assert InputDatabase.GENE_SYMBOL_ORDERED_LOCUS == "Gene Symbol Ordered Locus"
    assert InputDatabase.GENE_SYMBOL_ORF == "Gene Symbol ORF"
    assert InputDatabase.GI_NUMBER == "GI Number"
    assert InputDatabase.GO_ID == "GO ID"
    assert InputDatabase.GSEA_STANDARD_NAME == "GSEA Standard Name"
    assert InputDatabase.H_INV_LOCUS_ID == "H-Inv Locus ID"
    assert InputDatabase.H_INV_PROTEIN_ID == "H-Inv Protein ID"
    assert InputDatabase.H_INV_TRANSCRIPT_ID == "H-Inv Transcript ID"
    assert InputDatabase.HGNC_ID == "HGNC ID"
    assert InputDatabase.HMDB_METABOLITE == "HMDB Metabolite"
    assert InputDatabase.HOMOLOGENE_ID == "HomoloGene ID"
    assert InputDatabase.ILLUMINA_ID == "Illumina ID"
    assert InputDatabase.INTERPRO_ID == "InterPro ID"
    assert InputDatabase.IPI_ID == "IPI ID"
    assert InputDatabase.KEGG_COMPOUND_ID == "KEGG Compound ID"
    assert InputDatabase.KEGG_COMPOUND_NAME == "KEGG Compound Name"
    assert InputDatabase.KEGG_DISEASE_ID == "KEGG Disease ID"
    assert InputDatabase.KEGG_DRUG_ID == "KEGG Drug ID"
    assert InputDatabase.KEGG_DRUG_NAME == "KEGG Drug Name"
    assert InputDatabase.KEGG_GENE_ID == "KEGG Gene ID"
    assert InputDatabase.KEGG_PATHWAY_ID == "KEGG Pathway ID"
    assert InputDatabase.MAIZEGDB_ID == "MaizeGDB ID"
    assert InputDatabase.MGI_ID == "MGI ID"
    assert InputDatabase.MIM_ID == "MIM ID"
    assert InputDatabase.MIRBASE_ID == "miRBase ID"
    assert InputDatabase.MIRBASE_MATURE_MIRNA_ACC == "miRBase Mature miRNA Acc"
    assert InputDatabase.NCIPID_PATHWAY_NAME == "NCIPID Pathway Name"
    assert InputDatabase.ORGANISM_SCIENTIFIC_NAME == "Organism Scientific Name"
    assert InputDatabase.PDB_ID == "PDB ID"
    assert InputDatabase.PFAM_ID == "Pfam ID"
    assert InputDatabase.PHARMGKB_ID == "PharmGKB ID"
    assert InputDatabase.PUBCHEM_ID == "PubChem ID"
    assert InputDatabase.REACTOME_PATHWAY_NAME == "Reactome Pathway Name"
    assert InputDatabase.REFSEQ_GENOMIC_ACCESSION == "RefSeq Genomic Accession"
    assert InputDatabase.REFSEQ_MRNA_ACCESSION == "RefSeq mRNA Accession"
    assert InputDatabase.REFSEQ_PROTEIN_ACCESSION == "RefSeq Protein Accession"
    assert InputDatabase.SGD_ID == "SGD ID"
    assert InputDatabase.TAIR_ID == "TAIR ID"
    assert InputDatabase.TAXON_ID == "Taxon ID"
    assert InputDatabase.UNIGENE_ID == "UniGene ID"
    assert InputDatabase.UNIPROT_ACCESSION == "UniProt Accession"
    assert InputDatabase.UNIPROT_ENTRY_NAME == "UniProt Entry Name"
    assert InputDatabase.UNIPROT_PROTEIN_NAME == "UniProt Protein Name"
    assert InputDatabase.UNISTS_ID == "UniSTS ID"
