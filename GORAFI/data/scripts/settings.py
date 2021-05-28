"""
Global variables called during the GORAFI data preparation step.

@ May 2021
@ Asloudj Yanis
"""

## Raw data relative path:
raw_path = "./data/raw"
## Rdy2use data relative path:
rdy2use_path = "./data/rdy2use"

## URLs to download data files from and their respective resulting files:
go_obo_url = "http://current.geneontology.org/ontology/go-basic.obo"
go_obo_file = "%s/go-basic.obo" % raw_path

reactome_hierarchy_url = "https://reactome.org/download/current/ReactomePathwaysRelation.txt"
reactome_hierarchy_file = "%s/reactome_hierarchy.txt" % raw_path

reactome_label_url = "https://reactome.org/download/current/ReactomePathways.txt"
reactome_label_file = "%s/reactome_label.txt" % raw_path

reactome_annotation_url = "https://reactome.org/download/current/UniProt2Reactome.txt"
reactome_annotation_file = "%s/reactome_annotation.txt" % raw_path

# reactome obo file is not downloaded but generated.
reactome_obo_file = "%s/reactome.obo" % raw_path

hpo_annotation_url = "http://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt"
hpo_annotation_file = "%s/hpo_annotation.txt" % raw_path

hpo_obo_url = "http://purl.obolibrary.org/obo/hp.obo"
hpo_obo_file = "%s/hp.obo" % rdy2use_path