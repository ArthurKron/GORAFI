"""
For a given set of genes, identify the frequent itemsets involving multiple ontologies terms.
"""

from goatools import obo_parser
from settings import *
from mlxtend.frequent_patterns import apriori, fpgrowth, association_rules
from mlxtend.preprocessing.transactionencoder import TransactionEncoder
# from itemset_mining.two_phase_huim import TwoPhase

import csv
import json
import time as tm
import pandas as pd
import scripts.common as cmn
import scripts.ontology as onto

# - - - - - - - - -- - - - - - -- - - - - - -- - -- - - -- - -- - - -- - -- - -- - - -- - - -- - -  - - - -
# - - - - - - - - -- - - - - - -- - - - - - -- - -- - - -- - -- - - -- - -- - -- - - -- - - -- - -  - - - -
# - - - - - - - - -- - - - - - -- - - - - - -- - -- - - -- - -- - - -- - -- - -- - - -- - - -- - -  - - - -

#_________________________________________ I T E M S E T S __________________________________________

#_________________________________________ L O A D I N G

begin = tm.time()

# Convert gene ids into symbols
if symbol:
    reader = csv.reader(open("%s/%s_gene_symbol.csv" % (rdy2use_path, species), 'rt'))
    gene_symbol_id_dict = {}
    for row in reader:
        k, v = row
        gene_symbol_id_dict[k] = v
    genes = cmn.get_values_from_keys(genes, gene_symbol_id_dict)

ontologies = {
    'GO': obo_parser.GODag("%s/%s_go-basic.obo" % (rdy2use_path, species)), 
    'R-': obo_parser.GODag("%s/%s_reactome.obo" % (rdy2use_path, species))}
if species == "human":
    ontologies['HP'] = obo_parser.GODag(hpo_obo_file)


## extract the annotations corresponding to the selected genes:
with open("%s/%s_gene_annotation.json" % (rdy2use_path, species), 'rt') as jn:
    annotation = json.load(jn)
print(annotation[gene_symbol_id_dict["A8MVA2"]]['GO'])

common_keys = set(genes).intersection(set(annotation.keys()))
transactions = cmn.get_values_from_keys(common_keys, annotation)
te = TransactionEncoder()
trans_bool = te.fit(transactions).transform(transactions)
df = pd.DataFrame(trans_bool, columns=te.columns_)

for col in df.columns:
    term = ontologies[col[:2]][col]
    if col[:2] != "R-":
        if term.depth == 0:
            df = df.drop(columns=col)
    else:
        if term.depth * term.level <= term.level + term.depth:
            df = df.drop(columns=col)

end_load = tm.time()

#_________________________________________ M I N I N G

def has_same_ontology(itemset):
    for e in itemset:
        try:
            if prev_e[:2] != e[:2]:
                return False
        except NameError:
            pass
        prev_e = e
    return True

def str_from_iterable(iterable):
    iterable = list(iterable)
    string = iterable[0]
    for i in range(1, len(iterable)):
        string += ",%s" % iterable[i]
    return string

itemsets = fpgrowth(df, min_support = 0.25, max_len = 2, use_colnames = True)
rules = association_rules(itemsets)
rules["antecedents"] = rules["antecedents"].apply(lambda x: str_from_iterable(x))
rules["consequents"] = rules["consequents"].apply(lambda x: str_from_iterable(x))
rules.to_csv("rules.csv", index=False)

# interest_isets = list()
# for iset in itemsets:
#     if len(iset) > 1:
#         x, y = iset
#         x_term = ontologies[x[:2]][x]
#         y_term = ontologies[y[:2]][y]
#         if x[:2] != y[:2]:
#             interest_isets.append(iset)
# print(len(interest_isets))

end_min = tm.time()

print("\nitemset mining completed in %s seconds:\n \
    .loading: %ss\n \
    .mining: %ss\n" % (
        (end_min - begin), 
        (end_load - begin),
        (end_min - end_load)))