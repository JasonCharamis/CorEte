import seaborn as sns
from ete3 import Tree, TreeStyle, NodeStyle, TextFace
import re, random, os
import argparse

## Collection of functions to manipulate and visualize phylogenetic trees using the ETE3 toolkit ##

def parseCmdArguments():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Library for efficiently manipulating gff3 files.')
    parser.add_argument('-t','--tree', type=str, help='Nwk tree input file.')
    parser.add_argument('-m','--midpoint', action="store_true",help='Midpoint root tree.')
    parser.add_argument('--collapse', action="store_true",help='Collapse nodes based on bootstrap support.')
    parser.add_argument('--cutoff', type=str,help='Bootstrap support cutoff for collapsing nodes.')
    parser.add_argument('-r','--resolve', action="store_true",help='Resolve polytomies.')
    parser.add_argument('-v','--visualize', action="store_true",help='Visualize tree.')
    parser.add_argument('-l','--layout', type=str,help='Layout to visualize tree. Default: circular')
    parser.add_argument('-c','--count', action="store_true",help='Count leaves')
    parser.add_argument('-n','--names', action="store_true",help='Option to replace names in nwk')
    parser.add_argument('-nf','--names_file', type = str ,help='File with names to replace. Default is second column in tab-separated format.')
    args = parser.parse_args()

    if not any(vars(args).values()):
        parser.print_help()
        print("Error: No arguments provided.")

    return parser.parse_args()

## Implementation ##

def main():
    
    parser = argparse.ArgumentParser(description='Library for efficiently manipulating newick files and trees.')
    args = parseCmdArguments()
    
    if args.tree:

        inp = re.sub (".nwk$","",args.tree)
        
        if args.midpoint:
            midpoint(args.tree)

        elif args.collapse:
            if args.cutoff:
                bootstrap_collapse(args.tree,args.cutoff)
            else:
                print ( "Please provide a bootstrap cutoff.")
            
        elif args.visualize:
            visualize_tree(args.tree, layout=args.layout if args.layout else 'c' )

        elif args.resolve:
            resolve_polytomies(args.tree)

        elif args.names:
            if args.names_file:
                print ( sub_names_nwk(args.tree, args.names_file ) )
            else:
                print ( "Please provide a list with gene names to replace.")

        elif args.count:
            print ( count_leaves(args.tree) )
            
    else:
        print ("Please provide a nwk file as input.")


## Midpoint rooting ##
def midpoint(input):    
    tree = Tree(input, format = 1)
    
    ## get midpoint root of tree ##
    midpoint = tree.get_midpoint_outgroup()

    ## set midpoint root as outgroup ##
    tree.set_outgroup(midpoint)
    tree.write(format=1, outfile=input+".midpointed")
    return


def bootstrap_collapse(tree, threshold=50):
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


def count_leaves ( tree ):
    nleaves = []

    t = Tree(tree)
    
    for leaf in t.iter_leaves():
        nleaves.append(leaf)

    counts = len(nleaves)
    
    return counts
        

def color_node(node): ## function to color node and all descendants ##
    node.img_style["fgcolor"] = node_color
    for child in node.children:
        color_node(child)

## Main tree visualization function ##
def visualize_tree(tree, layout = "c"):
    t=Tree(tree)

    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.mode = layout

    species_colors = {}

    for leaf in t.iter_leaves(): ## Assign a unique color to each species ##
        leaf.name=re.sub ( "^g","Llongipalpis_g", leaf.name )
        leaf.name = re.sub ( "\..*", "", leaf.name )

        if re.search ( "_", leaf.name ): ## If _ is present, get the prefix as species name
            species = re.sub("_.*","", leaf.name)

        else: ## If no _ is found, get the species names from the first four letters of the gene ID
            species = leaf.name[0:3]
            
        if species not in species_colors:
            color = "#{:06x}".format(random.randint(0, 0xFFFFFF))
            if color != "#000000": ## keep black color for special identifiers ##
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

        if node.is_leaf(): ## if nodes are leaves, get color name based on species ##   
            if re.search ( "_", node.name ): ## If _ is present, get the prefix as species name
                species_n = re.sub("_.*","", node.name)
                
            else: ## If no _ is found, get the species names from the first four letters of the gene ID
                species_n = node.name[0:3]
                
            color_n = species_colors[species_n]               
            species_face = TextFace(node.name,fgcolor=color_n, fsize=500,ftype="Arial")
            node.add_face(species_face, column=1, position='branch-right')

        if node.support < 1: ## convert all bootstrap values to hundreds
            node.support = node.support*100

        for threshold, col in thresholds.items(): ## bootstrap values in internal nodes represented as circles ##
            if node.support <= 50:
                color_b = "lightgrey"
            elif node.support >= threshold:
                color_b = col
                
        if color_b:
            node_style=NodeStyle()
            node_style["fgcolor"] = color_b
            node_style["size"] = 500
            node.set_style(node_style)


    t.render(tree+".svg", w=5000, units="mm", tree_style=ts)
#    t.show(tree_style=ts)

    return


## Leaf counting functions ##
def count_descendant_leaves ( tree, node ):
    t=Tree(tree)
    descendant_leaves = []

    for node in t.traverse ("preorder"):
        descendant_leaves.append ( node.get_leaf_names() )
        
    print ("You have", len(descendant_leaves), "leaves" )
    return


def count_leaves_by_taxon ( tree, taxon_ID ):
    t=Tree(tree)
    descendant_leaves = []

    for node in t.traverse ("preorder"):
        if node.is_leaf():
            if re.search ( taxon_ID, node.name ):
                descendant_leaves.append ( node.get_leaf_names() )
                
    print ("You have", len(descendant_leaves), "leaves for", taxon_ID )
    return


def count_descendant_leaves_by_taxon ( tree, node, taxon_ID ):
    t=Tree(tree)
    descendant_leaves = []

    for node in t.traverse ("preorder"):
        if node.is_leaf():
            if re.search ( taxon_ID, node.name ):
                descendant_leaves.append ( node.get_leaf_names() )

    print ("You have", len(descendant_leaves), "descendant leaves for", taxon_ID )
    return


## Replace geneids in newick with names from a gene ID list ##

def fasta_names ( fasta ):
    gene_names = {}
    fasta=fasta.strip ("\n")
    id=re.sub(">|_.*","",fasta)
    name=re.sub(">|/","",fasta)
    gene_names[id]=name
    return gene_names


def get_newick(node):   
    if node.is_leaf():
        return f"{node.name}:{node.dist}"
    
    else:
        children_newick = ",".join([get_newick(child) for child in node.children])

        if hasattr(node, "support"):
            return f"({children_newick}){node.name}:{node.dist}[&support={node.support}]"
        else:
            return f"({children_newick}){node.name}:{node.dist}"

    return f"({children_newick}){node.name}:{node.dist}"


if __name__ == "__main__":
    main()

