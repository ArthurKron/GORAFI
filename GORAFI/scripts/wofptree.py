"""
Weighted Ordered FP-Tree Algorithm implementation.

WOFP:
Yuanyuan Li, Shaohong Yin (2020). 
Mining algorithm for weighted FP-tree frequent item sets based on two-dimensional table.
Journal of Physics: Conference Series.
DOI:10.1088/1742-6596/1453/1/012002

@ June 2021
@ Asloudj Yanis
"""

class Node:
    """Node of the FP tree."""

    def __init__(self, node_id, weight, parent, children):
        """A node has 4 domains: its node_id, its summed weight, its parent node and its children nodes."""
        self.node_id = node_id
        self.weight = weight
        self.parent = parent
        self.children = children

    def __str__(self):
        """Displays the node_id and the weight of the node"""
        return "%s {%s}" % (self.node_id, self.weight)

    def __repr__(self):
        """Display all attributes."""
        ret = [
            "node id: %s\n" % self.node_id,
            "weight: %s\n" % self.weight,
            "parent: %s\n" % self.parent,
            "children: %s" % self.children
        ]
        return "".join(ret)

# class FPTree(dict):
#     """FP tree dict-like object containing all the nodes."""



test = init_root
print()

def init_root():
    """
    # Description
    Returns the root of the FP Tree.

    # Usage
    >>> print(init_root())
    ... root {None}
    """
    return Node("root", None, None, None)

def get_all_nodes(node_id, fptree):
    """
    # Description
    Returns a list of nodes given a node id.

    #
    """