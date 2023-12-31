from enum import Enum


class InputDatabase(Enum):
    """
    These are valid input database types for the BioDBNet API.
    """
    AFFY_GENECHIP_ARRAY = "Affy GeneChip Array"
    AFFY_ID = "Affy ID"
    AFFY_TRANSCRIPT_CLUSTER_ID = "Affy Transcript Cluster ID"
    AGILENT_ID = "Agilent ID"
    BIOCARTA_PATHWAY_NAME = "Biocarta Pathway Name"
    CODELINK_ID = "CodeLink ID"
    DBSNP_ID = "dbSNP ID"
    DRUGBANK_DRUG_ID = "DrugBank Drug ID"
    DRUGBANK_DRUG_NAME = "DrugBank Drug Name"
    EC_NUMBER = "EC Number"
    ENSEMBL_GENE_ID = "Ensembl Gene ID"
    ENSEMBL_PROTEIN_ID = "Ensembl Protein ID"
    ENSEMBL_TRANSCRIPT_ID = "Ensembl Transcript ID"
    EST_ACCESSION = "EST Accession"
    FLYBASE_GENE_ID = "FlyBase Gene ID"
    GENBANK_NUCLEOTIDE_ACCESSION = "GenBank Nucleotide Accession"
    GENBANK_PROTEIN_ACCESSION = "GenBank Protein Accession"
    GENE_ID = "Gene ID"
    GENE_SYMBOL = "Gene Symbol"
    GENE_SYMBOL_AND_SYNONYMS = "Gene Symbol and Synonyms"
    GENE_SYMBOL_ORDERED_LOCUS = "Gene Symbol Ordered Locus"
    GENE_SYMBOL_ORF = "Gene Symbol ORF"
    GI_NUMBER = "GI Number"
    GO_ID = "GO ID"
    GSEA_STANDARD_NAME = "GSEA Standard Name"
    H_INV_LOCUS_ID = "H-Inv Locus ID"
    H_INV_PROTEIN_ID = "H-Inv Protein ID"
    H_INV_TRANSCRIPT_ID = "H-Inv Transcript ID"
    HGNC_ID = "HGNC ID"
    HMDB_METABOLITE = "HMDB Metabolite"
    HOMOLOGENE_ID = "HomoloGene ID"
    ILLUMINA_ID = "Illumina ID"
    INTERPRO_ID = "InterPro ID"
    IPI_ID = "IPI ID"
    KEGG_COMPOUND_ID = "KEGG Compound ID"
    KEGG_COMPOUND_NAME = "KEGG Compound Name"
    KEGG_DISEASE_ID = "KEGG Disease ID"
    KEGG_DRUG_ID = "KEGG Drug ID"
    KEGG_DRUG_NAME = "KEGG Drug Name"
    KEGG_GENE_ID = "KEGG Gene ID"
    KEGG_PATHWAY_ID = "KEGG Pathway ID"
    MAIZEGDB_ID = "MaizeGDB ID"
    MGI_ID = "MGI ID"
    MIM_ID = "MIM ID"
    MIRBASE_ID = "miRBase ID"
    MIRBASE_MATURE_MIRNA_ACC = "miRBase Mature miRNA Acc"
    NCIPID_PATHWAY_NAME = "NCIPID Pathway Name"
    ORGANISM_SCIENTIFIC_NAME = "Organism Scientific Name"
    PDB_ID = "PDB ID"
    PFAM_ID = "Pfam ID"
    PHARMGKB_ID = "PharmGKB ID"
    PUBCHEM_ID = "PubChem ID"
    REACTOME_PATHWAY_NAME = "Reactome Pathway Name"
    REFSEQ_GENOMIC_ACCESSION = "RefSeq Genomic Accession"
    REFSEQ_MRNA_ACCESSION = "RefSeq mRNA Accession"
    REFSEQ_PROTEIN_ACCESSION = "RefSeq Protein Accession"
    SGD_ID = "SGD ID"
    TAIR_ID = "TAIR ID"
    TAXON_ID = "Taxon ID"
    UNIGENE_ID = "UniGene ID"
    UNIPROT_ACCESSION = "UniProt Accession"
    UNIPROT_ENTRY_NAME = "UniProt Entry Name"
    UNIPROT_PROTEIN_NAME = "UniProt Protein Name"
    UNISTS_ID = "UniSTS ID"
