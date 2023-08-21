import numpy as np
import pandas as pd
import seaborn as sns
from ete3 import Tree
from ete3 import TreeStyle, NodeStyle, TreeNode, TextFace
import re
import random
import os

## this is a collection of ETE3 functions for custom use ##

## General functions for midpoint rooting, node collapse based on bootstrap values and polytomy resolving ##
def midpoint(input):    
    tree = Tree(input, format = 1)
    
    ## get midpoint root of tree ##
    midpoint = tree.get_midpoint_outgroup()

    ## set midpoint root as outgroup ##
    tree.set_outgroup(midpoint)
    tree.write(format=1, outfile=input+".midpoint_rooted")
    return


def bootstrap_collapse(tree, threshold):
        t=Tree(tree)
        for node in t.traverse():
                if node.support < threshold:
                        return node.delete()
                else:
                        return node

                   
def resolve_polytomies(input):   
    tree = Tree(input, format = 1)   
    tree.resolve_polytomy(recursive=True) ## resolve polytomies in tree ##
    tree.write(format=1, outfile=input+".resolved_polytomies")
    return



## Tree visualization functions ##

def color_node(node): ## function to color node and all descendants ##
    node.img_style["fgcolor"] = node_color
    for child in node.children:
        color_node(child)


## main visualization function for leaf coloring and bootstrap support ##
def visualize_tree(tree, reference, layout = "c"):
        t=Tree(tree)
        ts = TreeStyle()
        ts.show_leaf_name = False
        ts.mode = layout

        species_colors = {}

        ## Assign a unique color to each species ##
        for leaf in t.iter_leaves():
                if not re.search (reference, leaf.name ):
                        species = re.sub("_.*","", leaf.name)

                if species not in species_colors:
                        color = "#{:06x}".format(random.randint(0, 0xFFFFFF))
                        species_colors[species] = color

        thresholds = {
                50: "grey",
                75: "darkgrey",
                100: "black"
        }

        for node in t.traverse():
                nstyle=NodeStyle()
                nstyle["size"] = 0
                node.set_style(nstyle)

                ## Species associated with different colors ##
                if node.is_leaf():
                        if re.search(reference,node.name):
                                reference_face = TextFace(node.name,fgcolor="black", fsize=10,ftype="Arial")
                                node.add_face(reference_face, column=0, position='branch-right')

                        else:
                                species_n=re.sub("_.*","",node.name)
                                color_n = species_colors[species_n]
                                species_face = TextFace(node.name,fgcolor=color_n, fsize=10,ftype="Arial")
                                node.add_face(species_face, column=1, position='branch-right')

                ## Circles representing bootstrap value ##
                else:
                        for threshold, col in thresholds.items():
                                if node.support <= 50:
                                        color_b = "lightgrey"
                                elif node.support >= threshold:
                                        color_b = col

                        if color_b:
                                node_style=NodeStyle()
                                node_style["fgcolor"] = color_b
                                node_style["size"] = 5
                                node.set_style(node_style)


        # save tree to svg and show
        t.render(tree+".svg", w=500, units="mm", tree_style=ts)
        t.show(tree_style=ts)


def visualize_subtrees ( file ):
    if file.endswith(".nwk"):
        visualize_tree ( file, layout= "c" )

## Replace geneids in newick with names from a gene ID list ##

def fasta_names ( fasta ):
    gene_names = {}
    fasta=fasta.strip ("\n")
    id=re.sub(">|_.*","",fasta)
    name=re.sub(">|/","",fasta)
    gene_names[id]=name
    return gene_names


def sub_names_nwk(names, newick):
    gene_names = {}
    for file in open(names, "r").readlines():
        fasta = file.strip()
        gene_names.update(fasta_names(fasta))

    with open(newick+".names", "w") as file:

        new_nwk=""
        for file2 in open ( newick, "r" ).readlines():
            nwk = file2.strip()

            for ids in gene_names.keys():
                nwk = re.sub(ids,gene_names[ids],nwk)
                new_nwk = nwk

        return file.write(new_nwk)
