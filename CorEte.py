from ete3 import Tree, TreeStyle, NodeStyle, TreeNode, TextFace
import re
import random

## this is a collection of ETE3 functions for custom use ##

def visualize_tree(tree, reference, layout = "c"):
        t=Tree(tree)
        ts = TreeStyle()
        ts.show_leaf_name = False
        ts.mode = layout

        species_colors = {}

        ## Assign a unique color to each species ##
        for leaf in t.iter_leaves():
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
                                node_style["size"] = 5  # Hide default node shape
                                node.set_style(node_style)

        t.render(tree+".svg", w=500, units="mm", tree_style=ts)
        t.show(tree_style=ts)
 

def midpoint(input):    
    tree = Tree(input, format = 1)
    
    ## get midpoint root of tree ##
    midpoint = tree.get_midpoint_outgroup()

    ## set midpoint root as outgroup ##
    tree.set_outgroup(midpoint)
    tree.write(format=1, outfile=input+".midpoint_rooted")
