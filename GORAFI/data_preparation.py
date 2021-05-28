"""
Retrieves the Gene Ontology and Reactome data.
The GO OBO file is filtered according to the terms specific to a species found using the GAF file.
A Reactome OBO file is created using the hierarchy and the label Reactome files.
The Reactome OBO file is filtered according to the terms specific to a species found using the leaves Reactome file.

The files used are indicated on settings.py.

Gene Association File, GAF format (up to version 2.2):
http://geneontology.org/docs/go-annotation-file-gaf-format-2.2/

OBO format (up to version 1.2):
https://owlcollab.github.io/oboformat/doc/GO.format.obo-1_2.html

@ May 2021
@ Asloudj Yanis
"""

from goatools import obo_parser, associations
from Bio.UniProt.GOA import gafiterator
from data.scripts.settings import *

import time as tm
import pandas as pd
import data.scripts.ontology as onto
import data.scripts.common as cmn
import data.scripts.reactome as rc
import data.scripts.hpo as hpo


# the species of interest:
species = "chicken"
# GAF file path:
gaf_file = "%s/goa_%s.gaf" % (raw_path, species)

# - - - - - - - - -- - - - - - -- - - - - - -- - -- - - -- - -- - - -- - -- - -- - - -- - - -- - -  - - - -
# - - - - - - - - -- - - - - - -- - - - - - -- - -- - - -- - -- - - -- - -- - -- - - -- - - -- - -  - - - -
# - - - - - - - - -- - - - - - -- - - - - - -- - -- - - -- - -- - - -- - -- - -- - - -- - - -- - -  - - - -

#_________________________________________ P R E P A R A T I O N __________________________________________

#_________________________________________ D O W N L O A D I N G

begin = tm.time()

associations.dnld_annotation(gaf_file)
cmn.replace_first_line(gaf_file, "!gaf-version: 2.1\n")
cmn.download_url(go_obo_url, go_obo_file)
cmn.download_url(reactome_hierarchy_url, reactome_hierarchy_file)
cmn.download_url(reactome_label_url, reactome_label_file)
cmn.download_url(reactome_annotation_url, reactome_annotation_file)

if species == "human":
    cmn.download_url(hpo_annotation_url, hpo_annotation_file)
    cmn.replace_first_line(hpo_annotation_file, "")
    cmn.download_url(hpo_obo_url, hpo_obo_file)

end_dl = tm.time()

#_________________________________________ L O A D I N G

## GO ANNOTATIONS
gene_go_annotation = {}
gene_symbol_id = set()
with open(gaf_file, 'rt') as gaf:
    for anno in gafiterator(gaf):
        gene_symbol_id.add((anno['DB_Object_Symbol'], anno['DB_Object_ID']))
        try:
            gene_go_annotation[anno['DB_Object_ID']].add(anno['GO_ID'])
        except KeyError:
            gene_go_annotation[anno['DB_Object_ID']] = set([anno['GO_ID']])
gene_symbol_id_dict = dict((x, y) for x, y in gene_symbol_id)

## GO ONTOLOGY
go_onto = obo_parser.GODag(go_obo_file, optional_attrs = "relationship")

## REACTOME ANNOTATIONS
gene_reactome_annotation = rc.load_reactome_annotation(reactome_annotation_file)

## REACTOME ONTOLOGY
reacterm = rc.load_reacterm_dict(reactome_hierarchy_file, reactome_label_file)
onto.save_as_obo(reacterm, reactome_obo_file, "ontology: reactome")
react_onto = obo_parser.GODag(reactome_obo_file)

## Extract common genes between multiple annotations
common_genes = set(gene_go_annotation.keys()).intersection(set(gene_reactome_annotation.keys()))

if species == "human":    
    ## HPO ANNOTATIONS
    gene_hpo_annotation = hpo.load_hpo_annotation(hpo_annotation_file)
    gene_uniprotkb_hpo_annotation = []
    # replace gene symbols by their UniProtKB IDs
    for k in gene_hpo_annotation.keys():
        try:
            uniprotkb_id = gene_symbol_id_dict[k]
            gene_uniprotkb_hpo_annotation.append((uniprotkb_id, gene_hpo_annotation[k]))
        except KeyError:
            pass
    gene_hpo_annotation = dict((x, y) for x, y in gene_uniprotkb_hpo_annotation)
    common_genes = set(gene_hpo_annotation.keys()).intersection(common_genes)

    ## HPO ONTOLOGY
    hpo_onto = obo_parser.GODag(hpo_obo_file)

# terms to keep on each ontology
filtered_go_onto_keys = set()
filtered_react_onto_keys = set()

# get the GO, Reactome (and if species == human, HPO) nodes corresponding to each gene.
rdy2use_data = {}
for gene_id in common_genes:
    ancestry_goid = onto.get_ancestry_id(gene_id, gene_go_annotation, go_onto, relationship = True)
    filtered_go_onto_keys.update(ancestry_goid)

    ancestry_reactid = onto.get_ancestry_id(gene_id, gene_reactome_annotation, react_onto)
    filtered_react_onto_keys.update(ancestry_reactid)

    # associate each gene with its nodes.
    rdy2use_data[gene_id] = {'GO': ancestry_goid, 'Reactome': ancestry_reactid}

    if species == "human":
        ancestry_hpoid = onto.get_ancestry_id(gene_id, gene_hpo_annotation, hpo_onto)
        rdy2use_data[gene_id]['HPO'] = ancestry_hpoid

end_load = tm.time()

#_________________________________________ E X P O R T

## Filter the GO terms and the REACTerms according to the species
delete_go_keys = set(go_onto.keys()).difference(filtered_go_onto_keys)
for k in delete_go_keys:
    go_onto.pop(k)
delete_react_keys = set(react_onto.keys()).difference(filtered_react_onto_keys)
for k in delete_react_keys:
    react_onto.pop(k)

## Export the data as rdy2use files.
cmn.save_as_json(rdy2use_data, "%s/%s_gene_annotation.json" % (rdy2use_path, species))
onto.save_as_obo(go_onto, "%s/%s_go-basic.obo" % (rdy2use_path, species), "ontology: go")
onto.save_as_obo(react_onto, "%s/%s_reactome.obo" % (rdy2use_path, species), "ontology: reactome")


end_exp = tm.time()

print("\nPreparation completed in %s seconds:\n \
    .downloading: %ss\n \
    .loading: %ss\n \
    .exporting: %ss\n"  % (
        (end_exp - begin), 
        (end_dl - begin), 
        (end_load - end_dl),
        (end_exp - end_load)))
