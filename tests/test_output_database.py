from async_bioservices.output_database import OutputDatabase


def test_name():
    assert OutputDatabase.AFFY_ID.name == "AFFY_ID"
    assert OutputDatabase.AGILENT_ID.name == "AGILENT_ID"
    assert OutputDatabase.ALLERGOME_CODE.name == "ALLERGOME_CODE"
    assert OutputDatabase.APIDB_CRYPTODB_ID.name == "APIDB_CRYPTODB_ID"
    assert OutputDatabase.BIOCARTA_PATHWAY_NAME.name == "BIOCARTA_PATHWAY_NAME"
    assert OutputDatabase.BIOCYC_ID.name == "BIOCYC_ID"
    assert OutputDatabase.CCDS_ID.name == "CCDS_ID"
    assert OutputDatabase.CHROMOSOMAL_LOCATION.name == "CHROMOSOMAL_LOCATION"
    assert OutputDatabase.CLEANEX_ID.name == "CLEANEX_ID"
    assert OutputDatabase.CODELINK_ID.name == "CODELINK_ID"
    assert OutputDatabase.COSMIC_ID.name == "COSMIC_ID"
    assert OutputDatabase.CPDB_PROTEIN_INTERACTOR.name == "CPDB_PROTEIN_INTERACTOR"
    assert OutputDatabase.CTD_DISEASE_INFO.name == "CTD_DISEASE_INFO"
    assert OutputDatabase.CTD_DISEASE_NAME.name == "CTD_DISEASE_NAME"
    assert OutputDatabase.CYGD_ID.name == "CYGD_ID"
    assert OutputDatabase.DBSNP_ID.name == "DBSNP_ID"
    assert OutputDatabase.DICTYBASE_ID.name == "DICTYBASE_ID"
    assert OutputDatabase.DIP_ID.name == "DIP_ID"
    assert OutputDatabase.DISPROT_ID.name == "DISPROT_ID"
    assert OutputDatabase.DRUGBANK_DRUG_ID.name == "DRUGBANK_DRUG_ID"
    assert OutputDatabase.DRUGBANK_DRUG_INFO.name == "DRUGBANK_DRUG_INFO"
    assert OutputDatabase.DRUGBANK_DRUG_NAME.name == "DRUGBANK_DRUG_NAME"
    assert OutputDatabase.EC_NUMBER.name == "EC_NUMBER"
    assert OutputDatabase.ECHOBASE_ID.name == "ECHOBASE_ID"
    assert OutputDatabase.ECOGENE_ID.name == "ECOGENE_ID"
    assert OutputDatabase.ENSEMBL_BIOTYPE.name == "ENSEMBL_BIOTYPE"
    assert OutputDatabase.ENSEMBL_GENE_ID.name == "ENSEMBL_GENE_ID"
    assert OutputDatabase.ENSEMBL_GENE_INFO.name == "ENSEMBL_GENE_INFO"
    assert OutputDatabase.ENSEMBL_PROTEIN_ID.name == "ENSEMBL_PROTEIN_ID"
    assert OutputDatabase.ENSEMBL_TRANSCRIPT_ID.name == "ENSEMBL_TRANSCRIPT_ID"
    assert OutputDatabase.FLYBASE_GENE_ID.name == "FLYBASE_GENE_ID"
    assert OutputDatabase.FLYBASE_PROTEIN_ID.name == "FLYBASE_PROTEIN_ID"
    assert OutputDatabase.FLYBASE_TRANSCRIPT_ID.name == "FLYBASE_TRANSCRIPT_ID"
    assert OutputDatabase.GAD_DISEASE_INFO.name == "GAD_DISEASE_INFO"
    assert OutputDatabase.GAD_DISEASE_NAME.name == "GAD_DISEASE_NAME"
    assert OutputDatabase.GENBANK_NUCLEOTIDE_ACCESSION.name == "GENBANK_NUCLEOTIDE_ACCESSION"
    assert OutputDatabase.GENBANK_NUCLEOTIDE_GI.name == "GENBANK_NUCLEOTIDE_GI"
    assert OutputDatabase.GENBANK_PROTEIN_ACCESSION.name == "GENBANK_PROTEIN_ACCESSION"
    assert OutputDatabase.GENBANK_PROTEIN_GI.name == "GENBANK_PROTEIN_GI"
    assert OutputDatabase.GENE_ID.name == "GENE_ID"
    assert OutputDatabase.GENE_INFO.name == "GENE_INFO"
    assert OutputDatabase.GENE_SYMBOL.name == "GENE_SYMBOL"
    assert OutputDatabase.GENE_SYMBOL_AND_SYNONYMS.name == "GENE_SYMBOL_AND_SYNONYMS"
    assert OutputDatabase.GENE_SYMBOL_ORDERED_LOCUS.name == "GENE_SYMBOL_ORDERED_LOCUS"
    assert OutputDatabase.GENE_SYMBOL_ORF.name == "GENE_SYMBOL_ORF"
    assert OutputDatabase.GENE_SYNONYMS.name == "GENE_SYNONYMS"
    assert OutputDatabase.GENEFARM_ID.name == "GENEFARM_ID"
    assert OutputDatabase.GO_BIOLOGICAL_PROCESS.name == "GO_BIOLOGICAL_PROCESS"
    assert OutputDatabase.GO_CELLULAR_COMPONENT.name == "GO_CELLULAR_COMPONENT"
    assert OutputDatabase.GO_MOLECULAR_FUNCTION.name == "GO_MOLECULAR_FUNCTION"
    assert OutputDatabase.GO_ID.name == "GO_ID"
    assert OutputDatabase.GSEA_STANDARD_NAME.name == "GSEA_STANDARD_NAME"
    assert OutputDatabase.H_INV_LOCUS_ID.name == "H_INV_LOCUS_ID"
    assert OutputDatabase.HAMAP_ID.name == "HAMAP_ID"
    assert OutputDatabase.HGNC_ID.name == "HGNC_ID"
    assert OutputDatabase.HMDB_METABOLITE.name == "HMDB_METABOLITE"
    assert OutputDatabase.HOMOLOG_ALL_ENS_GENE_ID.name == "HOMOLOG_ALL_ENS_GENE_ID"
    assert OutputDatabase.HOMOLOG_ALL_ENS_PROTEIN_ID.name == "HOMOLOG_ALL_ENS_PROTEIN_ID"
    assert OutputDatabase.HOMOLOG_ALL_GENE_ID.name == "HOMOLOG_ALL_GENE_ID"
    assert OutputDatabase.HOMOLOG_HUMAN_ENS_GENE_ID.name == "HOMOLOG_HUMAN_ENS_GENE_ID"
    assert OutputDatabase.HOMOLOG_HUMAN_ENS_PROTEIN_ID.name == "HOMOLOG_HUMAN_ENS_PROTEIN_ID"
    assert OutputDatabase.HOMOLOG_HUMAN_GENE_ID.name == "HOMOLOG_HUMAN_GENE_ID"
    assert OutputDatabase.HOMOLOG_MOUSE_ENS_GENE_ID.name == "HOMOLOG_MOUSE_ENS_GENE_ID"
    assert OutputDatabase.HOMOLOG_MOUSE_ENS_PROTEIN_ID.name == "HOMOLOG_MOUSE_ENS_PROTEIN_ID"
    assert OutputDatabase.HOMOLOG_MOUSE_GENE_ID.name == "HOMOLOG_MOUSE_GENE_ID"
    assert OutputDatabase.HOMOLOG_RAT_ENS_GENE_ID.name == "HOMOLOG_RAT_ENS_GENE_ID"
    assert OutputDatabase.HOMOLOG_RAT_ENS_PROTEIN_ID.name == "HOMOLOG_RAT_ENS_PROTEIN_ID"
    assert OutputDatabase.HOMOLOG_RAT_GENE_ID.name == "HOMOLOG_RAT_GENE_ID"
    assert OutputDatabase.HOMOLOGENE_ID.name == "HOMOLOGENE_ID"
    assert OutputDatabase.HPA_ID.name == "HPA_ID"
    assert OutputDatabase.HPRD_ID.name == "HPRD_ID"
    assert OutputDatabase.HPRD_PROTEIN_COMPLEX.name == "HPRD_PROTEIN_COMPLEX"
    assert OutputDatabase.HPRD_PROTEIN_INTERACTOR.name == "HPRD_PROTEIN_INTERACTOR"
    assert OutputDatabase.ILLUMINA_ID.name == "ILLUMINA_ID"
    assert OutputDatabase.IMGT_GENE_DB_ID.name == "IMGT_GENE_DB_ID"
    assert OutputDatabase.INTERPRO_ID.name == "INTERPRO_ID"
    assert OutputDatabase.IPI_ID.name == "IPI_ID"
    assert OutputDatabase.KEGG_DISEASE_ID.name == "KEGG_DISEASE_ID"
    assert OutputDatabase.KEGG_GENE_ID.name == "KEGG_GENE_ID"
    assert OutputDatabase.KEGG_ORTHOLOGY_ID.name == "KEGG_ORTHOLOGY_ID"
    assert OutputDatabase.KEGG_PATHWAY_ID.name == "KEGG_PATHWAY_ID"
    assert OutputDatabase.KEGG_PATHWAY_INFO.name == "KEGG_PATHWAY_INFO"
    assert OutputDatabase.KEGG_PATHWAY_TITLE.name == "KEGG_PATHWAY_TITLE"
    assert OutputDatabase.LEGIOLIST_ID.name == "LEGIOLIST_ID"
    assert OutputDatabase.LEPROMA_ID.name == "LEPROMA_ID"
    assert OutputDatabase.LOCUS_TAG.name == "LOCUS_TAG"
    assert OutputDatabase.MAIZEGDB_ID.name == "MAIZEGDB_ID"
    assert OutputDatabase.MEROPS_ID.name == "MEROPS_ID"
    assert OutputDatabase.MGC_ZGC_XGC_ID.name == "MGC_ZGC_XGC_ID"
    assert OutputDatabase.MGC_ZGC_XGC_IMAGE_ID.name == "MGC_ZGC_XGC_IMAGE_ID"
    assert OutputDatabase.MGC_ZGC_XGC_INFO.name == "MGC_ZGC_XGC_INFO"
    assert OutputDatabase.MGI_ID.name == "MGI_ID"
    assert OutputDatabase.MIM_ID.name == "MIM_ID"
    assert OutputDatabase.MIM_INFO.name == "MIM_INFO"
    assert OutputDatabase.MIRBASE_ID.name == "MIRBASE_ID"
    assert OutputDatabase.NCIPID_PATHWAY_NAME.name == "NCIPID_PATHWAY_NAME"
    assert OutputDatabase.NCIPID_PROTEIN_COMPLEX.name == "NCIPID_PROTEIN_COMPLEX"
    assert OutputDatabase.NCIPID_PROTEIN_INTERACTOR.name == "NCIPID_PROTEIN_INTERACTOR"
    assert OutputDatabase.NCIPID_PTM.name == "NCIPID_PTM"
    assert OutputDatabase.ORPHANET_ID.name == "ORPHANET_ID"
    assert OutputDatabase.PANTHER_ID.name == "PANTHER_ID"
    assert OutputDatabase.PARALOG_ENS_GENE_ID.name == "PARALOG_ENS_GENE_ID"
    assert OutputDatabase.PBR_ID.name == "PBR_ID"
    assert OutputDatabase.PDB_ID.name == "PDB_ID"
    assert OutputDatabase.PEROXIBASE_ID.name == "PEROXIBASE_ID"
    assert OutputDatabase.PFAM_ID.name == "PFAM_ID"
    assert OutputDatabase.PHARMGKB_DRUG_INFO.name == "PHARMGKB_DRUG_INFO"
    assert OutputDatabase.PHARMGKB_GENE_ID.name == "PHARMGKB_GENE_ID"
    assert OutputDatabase.PIR_ID.name == "PIR_ID"
    assert OutputDatabase.PIRSF_ID.name == "PIRSF_ID"
    assert OutputDatabase.PPTASEDB_ID.name == "PPTASEDB_ID"
    assert OutputDatabase.PRINTS_ID.name == "PRINTS_ID"
    assert OutputDatabase.PRODOM_ID.name == "PRODOM_ID"
    assert OutputDatabase.PROSITE_ID.name == "PROSITE_ID"
    assert OutputDatabase.PSEUDOCAP_ID.name == "PSEUDOCAP_ID"
    assert OutputDatabase.PUBMED_ID.name == "PUBMED_ID"
    assert OutputDatabase.REACTOME_ID.name == "REACTOME_ID"
    assert OutputDatabase.REACTOME_PATHWAY_NAME.name == "REACTOME_PATHWAY_NAME"
    assert OutputDatabase.REBASE_ID.name == "REBASE_ID"
    assert OutputDatabase.REFSEQ_GENOMIC_ACCESSION.name == "REFSEQ_GENOMIC_ACCESSION"
    assert OutputDatabase.REFSEQ_GENOMIC_GI.name == "REFSEQ_GENOMIC_GI"
    assert OutputDatabase.REFSEQ_MRNA_ACCESSION.name == "REFSEQ_MRNA_ACCESSION"
    assert OutputDatabase.REFSEQ_NCRNA_ACCESSION.name == "REFSEQ_NCRNA_ACCESSION"
    assert OutputDatabase.REFSEQ_NUCLEOTIDE_GI.name == "REFSEQ_NUCLEOTIDE_GI"
    assert OutputDatabase.REFSEQ_PROTEIN_ACCESSION.name == "REFSEQ_PROTEIN_ACCESSION"
    assert OutputDatabase.REFSEQ_PROTEIN_GI.name == "REFSEQ_PROTEIN_GI"
    assert OutputDatabase.RFAM_ID.name == "RFAM_ID"
    assert OutputDatabase.RGD_ID.name == "RGD_ID"
    assert OutputDatabase.SGD_ID.name == "SGD_ID"
    assert OutputDatabase.SMART_ID.name == "SMART_ID"
    assert OutputDatabase.STRING_PROTEIN_INTERACTOR.name == "STRING_PROTEIN_INTERACTOR"
    assert OutputDatabase.TAIR_ID.name == "TAIR_ID"
    assert OutputDatabase.TAXON_ID.name == "TAXON_ID"
    assert OutputDatabase.TCDB_ID.name == "TCDB_ID"
    assert OutputDatabase.TIGRFAMS_ID.name == "TIGRFAMS_ID"
    assert OutputDatabase.TUBERCULIST_ID.name == "TUBERCULIST_ID"
    assert OutputDatabase.UCSC_ID.name == "UCSC_ID"
    assert OutputDatabase.UNIGENE_ID.name == "UNIGENE_ID"
    assert OutputDatabase.UNIPROT_ACCESSION.name == "UNIPROT_ACCESSION"
    assert OutputDatabase.UNIPROT_ENTRY_NAME.name == "UNIPROT_ENTRY_NAME"
    assert OutputDatabase.UNIPROT_INFO.name == "UNIPROT_INFO"
    assert OutputDatabase.UNIPROT_PROTEIN_NAME.name == "UNIPROT_PROTEIN_NAME"
    assert OutputDatabase.UNISTS_ID.name == "UNISTS_ID"
    assert OutputDatabase.VECTORBASE_GENE_ID.name == "VECTORBASE_GENE_ID"
    assert OutputDatabase.VEGA_GENE_ID.name == "VEGA_GENE_ID"
    assert OutputDatabase.VEGA_PROTEIN_ID.name == "VEGA_PROTEIN_ID"
    assert OutputDatabase.VEGA_TRANSCRIPT_ID.name == "VEGA_TRANSCRIPT_ID"
    assert OutputDatabase.WORMBASE_GENE_ID.name == "WORMBASE_GENE_ID"
    assert OutputDatabase.WORMPEP_PROTEIN_ID.name == "WORMPEP_PROTEIN_ID"
    assert OutputDatabase.XENBASE_GENE_ID.name == "XENBASE_GENE_ID"
    assert OutputDatabase.ZFIN_ID.name == "ZFIN_ID"


def test_type():
    assert isinstance(OutputDatabase.AFFY_ID, OutputDatabase)
    assert isinstance(OutputDatabase.AGILENT_ID, OutputDatabase)
    assert isinstance(OutputDatabase.ALLERGOME_CODE, OutputDatabase)
    assert isinstance(OutputDatabase.APIDB_CRYPTODB_ID, OutputDatabase)
    assert isinstance(OutputDatabase.BIOCARTA_PATHWAY_NAME, OutputDatabase)
    assert isinstance(OutputDatabase.BIOCYC_ID, OutputDatabase)
    assert isinstance(OutputDatabase.CCDS_ID, OutputDatabase)
    assert isinstance(OutputDatabase.CHROMOSOMAL_LOCATION, OutputDatabase)
    assert isinstance(OutputDatabase.CLEANEX_ID, OutputDatabase)
    assert isinstance(OutputDatabase.CODELINK_ID, OutputDatabase)
    assert isinstance(OutputDatabase.COSMIC_ID, OutputDatabase)
    assert isinstance(OutputDatabase.CPDB_PROTEIN_INTERACTOR, OutputDatabase)
    assert isinstance(OutputDatabase.CTD_DISEASE_INFO, OutputDatabase)
    assert isinstance(OutputDatabase.CTD_DISEASE_NAME, OutputDatabase)
    assert isinstance(OutputDatabase.CYGD_ID, OutputDatabase)
    assert isinstance(OutputDatabase.DBSNP_ID, OutputDatabase)
    assert isinstance(OutputDatabase.DICTYBASE_ID, OutputDatabase)
    assert isinstance(OutputDatabase.DIP_ID, OutputDatabase)
    assert isinstance(OutputDatabase.DISPROT_ID, OutputDatabase)
    assert isinstance(OutputDatabase.DRUGBANK_DRUG_ID, OutputDatabase)
    assert isinstance(OutputDatabase.DRUGBANK_DRUG_INFO, OutputDatabase)
    assert isinstance(OutputDatabase.DRUGBANK_DRUG_NAME, OutputDatabase)
    assert isinstance(OutputDatabase.EC_NUMBER, OutputDatabase)
    assert isinstance(OutputDatabase.ECHOBASE_ID, OutputDatabase)
    assert isinstance(OutputDatabase.ECOGENE_ID, OutputDatabase)
    assert isinstance(OutputDatabase.ENSEMBL_BIOTYPE, OutputDatabase)
    assert isinstance(OutputDatabase.ENSEMBL_GENE_ID, OutputDatabase)
    assert isinstance(OutputDatabase.ENSEMBL_GENE_INFO, OutputDatabase)
    assert isinstance(OutputDatabase.ENSEMBL_PROTEIN_ID, OutputDatabase)
    assert isinstance(OutputDatabase.ENSEMBL_TRANSCRIPT_ID, OutputDatabase)
    assert isinstance(OutputDatabase.FLYBASE_GENE_ID, OutputDatabase)
    assert isinstance(OutputDatabase.FLYBASE_PROTEIN_ID, OutputDatabase)
    assert isinstance(OutputDatabase.FLYBASE_TRANSCRIPT_ID, OutputDatabase)
    assert isinstance(OutputDatabase.GAD_DISEASE_INFO, OutputDatabase)
    assert isinstance(OutputDatabase.GAD_DISEASE_NAME, OutputDatabase)
    assert isinstance(OutputDatabase.GENBANK_NUCLEOTIDE_ACCESSION, OutputDatabase)
    assert isinstance(OutputDatabase.GENBANK_NUCLEOTIDE_GI, OutputDatabase)
    assert isinstance(OutputDatabase.GENBANK_PROTEIN_ACCESSION, OutputDatabase)
    assert isinstance(OutputDatabase.GENBANK_PROTEIN_GI, OutputDatabase)
    assert isinstance(OutputDatabase.GENE_ID, OutputDatabase)
    assert isinstance(OutputDatabase.GENE_INFO, OutputDatabase)
    assert isinstance(OutputDatabase.GENE_SYMBOL, OutputDatabase)
    assert isinstance(OutputDatabase.GENE_SYMBOL_AND_SYNONYMS, OutputDatabase)
    assert isinstance(OutputDatabase.GENE_SYMBOL_ORDERED_LOCUS, OutputDatabase)
    assert isinstance(OutputDatabase.GENE_SYMBOL_ORF, OutputDatabase)
    assert isinstance(OutputDatabase.GENE_SYNONYMS, OutputDatabase)
    assert isinstance(OutputDatabase.GENEFARM_ID, OutputDatabase)
    assert isinstance(OutputDatabase.GO_BIOLOGICAL_PROCESS, OutputDatabase)
    assert isinstance(OutputDatabase.GO_CELLULAR_COMPONENT, OutputDatabase)
    assert isinstance(OutputDatabase.GO_MOLECULAR_FUNCTION, OutputDatabase)
    assert isinstance(OutputDatabase.GO_ID, OutputDatabase)
    assert isinstance(OutputDatabase.GSEA_STANDARD_NAME, OutputDatabase)
    assert isinstance(OutputDatabase.H_INV_LOCUS_ID, OutputDatabase)
    assert isinstance(OutputDatabase.HAMAP_ID, OutputDatabase)
    assert isinstance(OutputDatabase.HGNC_ID, OutputDatabase)
    assert isinstance(OutputDatabase.HMDB_METABOLITE, OutputDatabase)
    assert isinstance(OutputDatabase.HOMOLOG_ALL_ENS_GENE_ID, OutputDatabase)
    assert isinstance(OutputDatabase.HOMOLOG_ALL_ENS_PROTEIN_ID, OutputDatabase)
    assert isinstance(OutputDatabase.HOMOLOG_ALL_GENE_ID, OutputDatabase)
    assert isinstance(OutputDatabase.HOMOLOG_HUMAN_ENS_GENE_ID, OutputDatabase)
    assert isinstance(OutputDatabase.HOMOLOG_HUMAN_ENS_PROTEIN_ID, OutputDatabase)
    assert isinstance(OutputDatabase.HOMOLOG_HUMAN_GENE_ID, OutputDatabase)
    assert isinstance(OutputDatabase.HOMOLOG_MOUSE_ENS_GENE_ID, OutputDatabase)
    assert isinstance(OutputDatabase.HOMOLOG_MOUSE_ENS_PROTEIN_ID, OutputDatabase)
    assert isinstance(OutputDatabase.HOMOLOG_MOUSE_GENE_ID, OutputDatabase)
    assert isinstance(OutputDatabase.HOMOLOG_RAT_ENS_GENE_ID, OutputDatabase)
    assert isinstance(OutputDatabase.HOMOLOG_RAT_ENS_PROTEIN_ID, OutputDatabase)
    assert isinstance(OutputDatabase.HOMOLOG_RAT_GENE_ID, OutputDatabase)
    assert isinstance(OutputDatabase.HOMOLOGENE_ID, OutputDatabase)
    assert isinstance(OutputDatabase.HPA_ID, OutputDatabase)
    assert isinstance(OutputDatabase.HPRD_ID, OutputDatabase)
    assert isinstance(OutputDatabase.HPRD_PROTEIN_COMPLEX, OutputDatabase)
    assert isinstance(OutputDatabase.HPRD_PROTEIN_INTERACTOR, OutputDatabase)
    assert isinstance(OutputDatabase.ILLUMINA_ID, OutputDatabase)
    assert isinstance(OutputDatabase.IMGT_GENE_DB_ID, OutputDatabase)
    assert isinstance(OutputDatabase.INTERPRO_ID, OutputDatabase)
    assert isinstance(OutputDatabase.IPI_ID, OutputDatabase)
    assert isinstance(OutputDatabase.KEGG_DISEASE_ID, OutputDatabase)
    assert isinstance(OutputDatabase.KEGG_GENE_ID, OutputDatabase)
    assert isinstance(OutputDatabase.KEGG_ORTHOLOGY_ID, OutputDatabase)
    assert isinstance(OutputDatabase.KEGG_PATHWAY_ID, OutputDatabase)
    assert isinstance(OutputDatabase.KEGG_PATHWAY_INFO, OutputDatabase)
    assert isinstance(OutputDatabase.KEGG_PATHWAY_TITLE, OutputDatabase)
    assert isinstance(OutputDatabase.LEGIOLIST_ID, OutputDatabase)
    assert isinstance(OutputDatabase.LEPROMA_ID, OutputDatabase)
    assert isinstance(OutputDatabase.LOCUS_TAG, OutputDatabase)
    assert isinstance(OutputDatabase.MAIZEGDB_ID, OutputDatabase)
    assert isinstance(OutputDatabase.MEROPS_ID, OutputDatabase)
    assert isinstance(OutputDatabase.MGC_ZGC_XGC_ID, OutputDatabase)
    assert isinstance(OutputDatabase.MGC_ZGC_XGC_IMAGE_ID, OutputDatabase)
    assert isinstance(OutputDatabase.MGC_ZGC_XGC_INFO, OutputDatabase)
    assert isinstance(OutputDatabase.MGI_ID, OutputDatabase)
    assert isinstance(OutputDatabase.MIM_ID, OutputDatabase)
    assert isinstance(OutputDatabase.MIM_INFO, OutputDatabase)
    assert isinstance(OutputDatabase.MIRBASE_ID, OutputDatabase)
    assert isinstance(OutputDatabase.NCIPID_PATHWAY_NAME, OutputDatabase)
    assert isinstance(OutputDatabase.NCIPID_PROTEIN_COMPLEX, OutputDatabase)
    assert isinstance(OutputDatabase.NCIPID_PROTEIN_INTERACTOR, OutputDatabase)
    assert isinstance(OutputDatabase.NCIPID_PTM, OutputDatabase)
    assert isinstance(OutputDatabase.ORPHANET_ID, OutputDatabase)
    assert isinstance(OutputDatabase.PANTHER_ID, OutputDatabase)
    assert isinstance(OutputDatabase.PARALOG_ENS_GENE_ID, OutputDatabase)
    assert isinstance(OutputDatabase.PBR_ID, OutputDatabase)
    assert isinstance(OutputDatabase.PDB_ID, OutputDatabase)
    assert isinstance(OutputDatabase.PEROXIBASE_ID, OutputDatabase)
    assert isinstance(OutputDatabase.PFAM_ID, OutputDatabase)
    assert isinstance(OutputDatabase.PHARMGKB_DRUG_INFO, OutputDatabase)
    assert isinstance(OutputDatabase.PHARMGKB_GENE_ID, OutputDatabase)
    assert isinstance(OutputDatabase.PIR_ID, OutputDatabase)
    assert isinstance(OutputDatabase.PIRSF_ID, OutputDatabase)
    assert isinstance(OutputDatabase.PPTASEDB_ID, OutputDatabase)
    assert isinstance(OutputDatabase.PRINTS_ID, OutputDatabase)
    assert isinstance(OutputDatabase.PRODOM_ID, OutputDatabase)
    assert isinstance(OutputDatabase.PROSITE_ID, OutputDatabase)
    assert isinstance(OutputDatabase.PSEUDOCAP_ID, OutputDatabase)
    assert isinstance(OutputDatabase.PUBMED_ID, OutputDatabase)
    assert isinstance(OutputDatabase.REACTOME_ID, OutputDatabase)
    assert isinstance(OutputDatabase.REACTOME_PATHWAY_NAME, OutputDatabase)
    assert isinstance(OutputDatabase.REBASE_ID, OutputDatabase)
    assert isinstance(OutputDatabase.REFSEQ_GENOMIC_ACCESSION, OutputDatabase)
    assert isinstance(OutputDatabase.REFSEQ_GENOMIC_GI, OutputDatabase)
    assert isinstance(OutputDatabase.REFSEQ_MRNA_ACCESSION, OutputDatabase)
    assert isinstance(OutputDatabase.REFSEQ_NCRNA_ACCESSION, OutputDatabase)
    assert isinstance(OutputDatabase.REFSEQ_NUCLEOTIDE_GI, OutputDatabase)
    assert isinstance(OutputDatabase.REFSEQ_PROTEIN_ACCESSION, OutputDatabase)
    assert isinstance(OutputDatabase.REFSEQ_PROTEIN_GI, OutputDatabase)
    assert isinstance(OutputDatabase.RFAM_ID, OutputDatabase)
    assert isinstance(OutputDatabase.RGD_ID, OutputDatabase)
    assert isinstance(OutputDatabase.SGD_ID, OutputDatabase)
    assert isinstance(OutputDatabase.SMART_ID, OutputDatabase)
    assert isinstance(OutputDatabase.STRING_PROTEIN_INTERACTOR, OutputDatabase)
    assert isinstance(OutputDatabase.TAIR_ID, OutputDatabase)
    assert isinstance(OutputDatabase.TAXON_ID, OutputDatabase)
    assert isinstance(OutputDatabase.TCDB_ID, OutputDatabase)
    assert isinstance(OutputDatabase.TIGRFAMS_ID, OutputDatabase)
    assert isinstance(OutputDatabase.TUBERCULIST_ID, OutputDatabase)
    assert isinstance(OutputDatabase.UCSC_ID, OutputDatabase)
    assert isinstance(OutputDatabase.UNIGENE_ID, OutputDatabase)
    assert isinstance(OutputDatabase.UNIPROT_ACCESSION, OutputDatabase)
    assert isinstance(OutputDatabase.UNIPROT_ENTRY_NAME, OutputDatabase)
    assert isinstance(OutputDatabase.UNIPROT_INFO, OutputDatabase)
    assert isinstance(OutputDatabase.UNIPROT_PROTEIN_NAME, OutputDatabase)
    assert isinstance(OutputDatabase.UNISTS_ID, OutputDatabase)
    assert isinstance(OutputDatabase.VECTORBASE_GENE_ID, OutputDatabase)
    assert isinstance(OutputDatabase.VEGA_GENE_ID, OutputDatabase)
    assert isinstance(OutputDatabase.VEGA_PROTEIN_ID, OutputDatabase)
    assert isinstance(OutputDatabase.VEGA_TRANSCRIPT_ID, OutputDatabase)
    assert isinstance(OutputDatabase.WORMBASE_GENE_ID, OutputDatabase)
    assert isinstance(OutputDatabase.WORMPEP_PROTEIN_ID, OutputDatabase)
    assert isinstance(OutputDatabase.XENBASE_GENE_ID, OutputDatabase)
    assert isinstance(OutputDatabase.ZFIN_ID, OutputDatabase)


def test_value():
    assert OutputDatabase.AFFY_ID.value == "Affy ID"
    assert OutputDatabase.AGILENT_ID.value == "Agilent ID"
    assert OutputDatabase.ALLERGOME_CODE.value == "Allergome Code"
    assert OutputDatabase.APIDB_CRYPTODB_ID.value == "ApiDB_CryptoDB ID"
    assert OutputDatabase.BIOCARTA_PATHWAY_NAME.value == "Biocarta Pathway Name"
    assert OutputDatabase.BIOCYC_ID.value == "BioCyc ID"
    assert OutputDatabase.CCDS_ID.value == "CCDS ID"
    assert OutputDatabase.CHROMOSOMAL_LOCATION.value == "Chromosomal Location"
    assert OutputDatabase.CLEANEX_ID.value == "CleanEx ID"
    assert OutputDatabase.CODELINK_ID.value == "CodeLink ID"
    assert OutputDatabase.COSMIC_ID.value == "COSMIC ID"
    assert OutputDatabase.CPDB_PROTEIN_INTERACTOR.value == "CPDB Protein Interactor"
    assert OutputDatabase.CTD_DISEASE_INFO.value == "CTD Disease Info"
    assert OutputDatabase.CTD_DISEASE_NAME.value == "CTD Disease Name"
    assert OutputDatabase.CYGD_ID.value == "CYGD ID"
    assert OutputDatabase.DBSNP_ID.value == "dbSNP ID"
    assert OutputDatabase.DICTYBASE_ID.value == "dictyBase ID"
    assert OutputDatabase.DIP_ID.value == "DIP ID"
    assert OutputDatabase.DISPROT_ID.value == "DisProt ID"
    assert OutputDatabase.DRUGBANK_DRUG_ID.value == "DrugBank Drug ID"
    assert OutputDatabase.DRUGBANK_DRUG_INFO.value == "DrugBank Drug Info"
    assert OutputDatabase.DRUGBANK_DRUG_NAME.value == "DrugBank Drug Name"
    assert OutputDatabase.EC_NUMBER.value == "EC Number"
    assert OutputDatabase.ECHOBASE_ID.value == "EchoBASE ID"
    assert OutputDatabase.ECOGENE_ID.value == "EcoGene ID"
    assert OutputDatabase.ENSEMBL_BIOTYPE.value == "Ensembl Biotype"
    assert OutputDatabase.ENSEMBL_GENE_ID.value == "Ensembl Gene ID"
    assert OutputDatabase.ENSEMBL_GENE_INFO.value == "Ensembl Gene Info"
    assert OutputDatabase.ENSEMBL_PROTEIN_ID.value == "Ensembl Protein ID"
    assert OutputDatabase.ENSEMBL_TRANSCRIPT_ID.value == "Ensembl Transcript ID"
    assert OutputDatabase.FLYBASE_GENE_ID.value == "FlyBase Gene ID"
    assert OutputDatabase.FLYBASE_PROTEIN_ID.value == "FlyBase Protein ID"
    assert OutputDatabase.FLYBASE_TRANSCRIPT_ID.value == "FlyBase Transcript ID"
    assert OutputDatabase.GAD_DISEASE_INFO.value == "GAD Disease Info"
    assert OutputDatabase.GAD_DISEASE_NAME.value == "GAD Disease Name"
    assert OutputDatabase.GENBANK_NUCLEOTIDE_ACCESSION.value == "GenBank Nucleotide Accession"
    assert OutputDatabase.GENBANK_NUCLEOTIDE_GI.value == "GenBank Nucleotide GI"
    assert OutputDatabase.GENBANK_PROTEIN_ACCESSION.value == "GenBank Protein Accession"
    assert OutputDatabase.GENBANK_PROTEIN_GI.value == "GenBank Protein GI"
    assert OutputDatabase.GENE_ID.value == "Gene ID"
    assert OutputDatabase.GENE_INFO.value == "Gene Info"
    assert OutputDatabase.GENE_SYMBOL.value == "Gene Symbol"
    assert OutputDatabase.GENE_SYMBOL_AND_SYNONYMS.value == "Gene Symbol and Synonyms"
    assert OutputDatabase.GENE_SYMBOL_ORDERED_LOCUS.value == "Gene Symbol Ordered Locus"
    assert OutputDatabase.GENE_SYMBOL_ORF.value == "Gene Symbol ORF"
    assert OutputDatabase.GENE_SYNONYMS.value == "Gene Synonyms"
    assert OutputDatabase.GENEFARM_ID.value == "GeneFarm ID"
    assert OutputDatabase.GO_BIOLOGICAL_PROCESS.value == "GO - Biological Process"
    assert OutputDatabase.GO_CELLULAR_COMPONENT.value == "GO - Cellular Component"
    assert OutputDatabase.GO_MOLECULAR_FUNCTION.value == "GO - Molecular Function"
    assert OutputDatabase.GO_ID.value == "GO ID"
    assert OutputDatabase.GSEA_STANDARD_NAME.value == "GSEA Standard Name"
    assert OutputDatabase.H_INV_LOCUS_ID.value == "H-Inv Locus ID"
    assert OutputDatabase.HAMAP_ID.value == "HAMAP ID"
    assert OutputDatabase.HGNC_ID.value == "HGNC ID"
    assert OutputDatabase.HMDB_METABOLITE.value == "HMDB Metabolite"
    assert OutputDatabase.HOMOLOG_ALL_ENS_GENE_ID.value == "Homolog - All Ens Gene ID"
    assert OutputDatabase.HOMOLOG_ALL_ENS_PROTEIN_ID.value == "Homolog - All Ens Protein ID"
    assert OutputDatabase.HOMOLOG_ALL_GENE_ID.value == "Homolog - All Gene ID"
    assert OutputDatabase.HOMOLOG_HUMAN_ENS_GENE_ID.value == "Homolog - Human Ens Gene ID"
    assert OutputDatabase.HOMOLOG_HUMAN_ENS_PROTEIN_ID.value == "Homolog - Human Ens Protein ID"
    assert OutputDatabase.HOMOLOG_HUMAN_GENE_ID.value == "Homolog - Human Gene ID"
    assert OutputDatabase.HOMOLOG_MOUSE_ENS_GENE_ID.value == "Homolog - Mouse Ens Gene ID"
    assert OutputDatabase.HOMOLOG_MOUSE_ENS_PROTEIN_ID.value == "Homolog - Mouse Ens Protein ID"
    assert OutputDatabase.HOMOLOG_MOUSE_GENE_ID.value == "Homolog - Mouse Gene ID"
    assert OutputDatabase.HOMOLOG_RAT_ENS_GENE_ID.value == "Homolog - Rat Ens Gene ID"
    assert OutputDatabase.HOMOLOG_RAT_ENS_PROTEIN_ID.value == "Homolog - Rat Ens Protein ID"
    assert OutputDatabase.HOMOLOG_RAT_GENE_ID.value == "Homolog - Rat Gene ID"
    assert OutputDatabase.HOMOLOGENE_ID.value == "HomoloGene ID"
    assert OutputDatabase.HPA_ID.value == "HPA ID"
    assert OutputDatabase.HPRD_ID.value == "HPRD ID"
    assert OutputDatabase.HPRD_PROTEIN_COMPLEX.value == "HPRD Protein Complex"
    assert OutputDatabase.HPRD_PROTEIN_INTERACTOR.value == "HPRD Protein Interactor"
    assert OutputDatabase.ILLUMINA_ID.value == "Illumina ID"
    assert OutputDatabase.IMGT_GENE_DB_ID.value == "IMGT/GENE-DB ID"
    assert OutputDatabase.INTERPRO_ID.value == "InterPro ID"
    assert OutputDatabase.IPI_ID.value == "IPI ID"
    assert OutputDatabase.KEGG_DISEASE_ID.value == "KEGG Disease ID"
    assert OutputDatabase.KEGG_GENE_ID.value == "KEGG Gene ID"
    assert OutputDatabase.KEGG_ORTHOLOGY_ID.value == "KEGG Orthology ID"
    assert OutputDatabase.KEGG_PATHWAY_ID.value == "KEGG Pathway ID"
    assert OutputDatabase.KEGG_PATHWAY_INFO.value == "KEGG Pathway Info"
    assert OutputDatabase.KEGG_PATHWAY_TITLE.value == "KEGG Pathway Title"
    assert OutputDatabase.LEGIOLIST_ID.value == "LegioList ID"
    assert OutputDatabase.LEPROMA_ID.value == "Leproma ID"
    assert OutputDatabase.LOCUS_TAG.value == "Locus Tag"
    assert OutputDatabase.MAIZEGDB_ID.value == "MaizeGDB ID"
    assert OutputDatabase.MEROPS_ID.value == "MEROPS ID"
    assert OutputDatabase.MGC_ZGC_XGC_ID.value == "MGC(ZGC/XGC) ID"
    assert OutputDatabase.MGC_ZGC_XGC_IMAGE_ID.value == "MGC(ZGC/XGC) Image ID"
    assert OutputDatabase.MGC_ZGC_XGC_INFO.value == "MGC(ZGC/XGC) Info"
    assert OutputDatabase.MGI_ID.value == "MGI ID"
    assert OutputDatabase.MIM_ID.value == "MIM ID"
    assert OutputDatabase.MIM_INFO.value == "MIM Info"
    assert OutputDatabase.MIRBASE_ID.value == "miRBase ID"
    assert OutputDatabase.NCIPID_PATHWAY_NAME.value == "NCIPID Pathway Name"
    assert OutputDatabase.NCIPID_PROTEIN_COMPLEX.value == "NCIPID Protein Complex"
    assert OutputDatabase.NCIPID_PROTEIN_INTERACTOR.value == "NCIPID Protein Interactor"
    assert OutputDatabase.NCIPID_PTM.value == "NCIPID PTM"
    assert OutputDatabase.ORPHANET_ID.value == "Orphanet ID"
    assert OutputDatabase.PANTHER_ID.value == "PANTHER ID"
    assert OutputDatabase.PARALOG_ENS_GENE_ID.value == "Paralog - Ens Gene ID"
    assert OutputDatabase.PBR_ID.value == "PBR ID"
    assert OutputDatabase.PDB_ID.value == "PDB ID"
    assert OutputDatabase.PEROXIBASE_ID.value == "PeroxiBase ID"
    assert OutputDatabase.PFAM_ID.value == "Pfam ID"
    assert OutputDatabase.PHARMGKB_DRUG_INFO.value == "PharmGKB Drug Info"
    assert OutputDatabase.PHARMGKB_GENE_ID.value == "PharmGKB Gene ID"
    assert OutputDatabase.PIR_ID.value == "PIR ID"
    assert OutputDatabase.PIRSF_ID.value == "PIRSF ID"
    assert OutputDatabase.PPTASEDB_ID.value == "PptaseDB ID"
    assert OutputDatabase.PRINTS_ID.value == "PRINTS ID"
    assert OutputDatabase.PRODOM_ID.value == "ProDom ID"
    assert OutputDatabase.PROSITE_ID.value == "PROSITE ID"
    assert OutputDatabase.PSEUDOCAP_ID.value == "PseudoCAP ID"
    assert OutputDatabase.PUBMED_ID.value == "PubMed ID"
    assert OutputDatabase.REACTOME_ID.value == "Reactome ID"
    assert OutputDatabase.REACTOME_PATHWAY_NAME.value == "Reactome Pathway Name"
    assert OutputDatabase.REBASE_ID.value == "REBASE ID"
    assert OutputDatabase.REFSEQ_GENOMIC_ACCESSION.value == "RefSeq Genomic Accession"
    assert OutputDatabase.REFSEQ_GENOMIC_GI.value == "RefSeq Genomic GI"
    assert OutputDatabase.REFSEQ_MRNA_ACCESSION.value == "RefSeq mRNA Accession"
    assert OutputDatabase.REFSEQ_NCRNA_ACCESSION.value == "RefSeq ncRNA Accession"
    assert OutputDatabase.REFSEQ_NUCLEOTIDE_GI.value == "RefSeq Nucleotide GI"
    assert OutputDatabase.REFSEQ_PROTEIN_ACCESSION.value == "RefSeq Protein Accession"
    assert OutputDatabase.REFSEQ_PROTEIN_GI.value == "RefSeq Protein GI"
    assert OutputDatabase.RFAM_ID.value == "Rfam ID"
    assert OutputDatabase.RGD_ID.value == "RGD ID"
    assert OutputDatabase.SGD_ID.value == "SGD ID"
    assert OutputDatabase.SMART_ID.value == "SMART ID"
    assert OutputDatabase.STRING_PROTEIN_INTERACTOR.value == "STRING Protein Interactor"
    assert OutputDatabase.TAIR_ID.value == "TAIR ID"
    assert OutputDatabase.TAXON_ID.value == "Taxon ID"
    assert OutputDatabase.TCDB_ID.value == "TCDB ID"
    assert OutputDatabase.TIGRFAMS_ID.value == "TIGRFAMs ID"
    assert OutputDatabase.TUBERCULIST_ID.value == "TubercuList ID"
    assert OutputDatabase.UCSC_ID.value == "UCSC ID"
    assert OutputDatabase.UNIGENE_ID.value == "UniGene ID"
    assert OutputDatabase.UNIPROT_ACCESSION.value == "UniProt Accession"
    assert OutputDatabase.UNIPROT_ENTRY_NAME.value == "UniProt Entry Name"
    assert OutputDatabase.UNIPROT_INFO.value == "UniProt Info"
    assert OutputDatabase.UNIPROT_PROTEIN_NAME.value == "UniProt Protein Name"
    assert OutputDatabase.UNISTS_ID.value == "UniSTS ID"
    assert OutputDatabase.VECTORBASE_GENE_ID.value == "VectorBase Gene ID"
    assert OutputDatabase.VEGA_GENE_ID.value == "VEGA Gene ID"
    assert OutputDatabase.VEGA_PROTEIN_ID.value == "VEGA Protein ID"
    assert OutputDatabase.VEGA_TRANSCRIPT_ID.value == "VEGA Transcript ID"
    assert OutputDatabase.WORMBASE_GENE_ID.value == "WormBase Gene ID"
    assert OutputDatabase.WORMPEP_PROTEIN_ID.value == "WormPep Protein ID"
    assert OutputDatabase.XENBASE_GENE_ID.value == "XenBase Gene ID"
    assert OutputDatabase.ZFIN_ID.value == "ZFIN ID"
