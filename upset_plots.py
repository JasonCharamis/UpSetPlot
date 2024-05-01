import os
import re
import argparse
import pandas as pd
import numpy as np
from upsetplot import UpSet
from matplotlib import pyplot as plt
from matplotlib import cm

# Parse command-line arguments

def parse_arguments():

    parser = argparse.ArgumentParser(description='Create Upset plots from a collection of lists.')
    parser.add_argument('-de', '--de_files', type=str, required=False, help='List of de files to extract subsets of up-regulated and down-regulated genes.')
    parser.add_argument('-f', '--field', type=str, required=False, help='Field in de files to find Fold Change (FC) values and define up-regulated and down-regulated genes. Default: second field (1).')
    parser.add_argument('-ls', '--lists_of_strings', type=str,required=False, help='Lists of strings, in the case of not provided de files.')
    parser.add_argument('-out', '--outfile', type=str, required=False, help='Output file to save results supporting the UpSet plot. Default: "UpSet_genes.tsv" ')
    parser.add_argument('-plt', '--plot', type=str, required=False, help='Output png file to print UpSet plot. Default: "UpSet_plot.svg" ')
    parser.add_argument('-img', '--image_format', required=False, type=str, help='Format to save image. Default: "svg" ')
    parser.add_argument('-st', '--stacked', required=False, type=str, help='List of regular expressions for identifying categories for stacked barplots. The second column corresponds to the classifying group.')
    parser.add_argument('-cv', '--categorical_variable', required=False, type=str, help='Categorical column to classify groups. Default: annotation')
        
    args = parser.parse_args()

    if not any(vars(args).values()):
        parser.print_help()

    return args


def main():

    args = parse_arguments()

    if args.de_files and not args.lists_of_strings:        
        if not os.path.exists(args.de_files):
            print(f"Error: File '{args.de_files}' not found.")
            return

        else:
            if args.stacked:
                upset_plots(
                    input_files = args.de_files.strip('\n'),
                    plot = args.plot if args.plot else 'UpSet_plot.svg',
                    image_format = args.image_format if args.image_format else 'svg',
                    de = True,
                    stacked = args.stacked if args.stacked else None,
                    categorical_variable = args.categorical_variable if args.categorical_variable else 'Annotation'
                )

            else:
                if args.categorical_variable:
                    print (f"Warning: Categorical variable '{args.categorical_variable}' is ignored when not using stacked barplots.")
                    pass
                else:
                    upset_plots(
                        input_files = args.de_files.strip('\n'),
                        plot = args.plot if args.plot else 'UpSet_plot.svg',
                        image_format = args.image_format if args.image_format else 'svg',
                        de = True
                    )
                
    elif not args.de_files and args.lists_of_strings:
        upset_plots(
            input_files = args.de_files,
            plot = args.plot if args.plot else 'UpSet_plot.svg',
            image_format = args.image_format if args.image_format else 'svg',
            field = args.field if args.field else '0',
            de = False,
            stacked = args.stacked if args.stacked else None,
            categorical_variable = args.categorical_variable if args.categorical_variable else 'Annotation'
        )

    elif args.de_files and args.lists_of_strings:
        print ( "Please select either de files or lists of strings.")

    else:
        print ("Please provide a list of files as input.")
        


def classify_subsets(categorical_variable, stacked_categories):
    with open(stacked_categories, "r") as file:
        categories = [line.strip().split('\t') for line in file]

    classified_subsets = []
    for cv in categorical_variable:
        category_found = False
        
        for category in categories:
            if re.search(category[0], cv):
                classified_subsets.append(category[1])
                category_found = True
                break
        if not category_found:
            classified_subsets.append('Other')

    return classified_subsets


def generate_lists_of_de_genes(filenames, stacked=None, categorical_variable = 'Annotation'):
    subsets = []
    dfs = {}

    with open(filenames, "r") as fls:
        fl = fls.readlines()

        for filename in fl:
            filename = filename.strip('\n')
            filename_regex = re.compile(r"\\..\\/|\\\..\*")
            comparison = re.sub(filename_regex, "", filename)  # Get name of the comparison from filename                                                                                                                                                             

            # Load pandas df with the contents of the de file and remove lines with Up- and Down-regulate                                                                                                                                                             
            de_file = pd.read_csv(filename, sep='\t', header=0)
            de_file['Geneid'] = de_file['Geneid'].astype(str)
            filtered_de_file = de_file[~de_file['Geneid'].str.contains('Upregulated|Downregulated')]

            # Identify rows corresponding to up- and down-regulated genes                                                                                                                                                                                             
            filtered_de_file['Fold_Change'] = filtered_de_file['Fold_Change'].astype(float)

            upregulated = filtered_de_file[filtered_de_file['Fold_Change'] > 0]
            downregulated = filtered_de_file[filtered_de_file['Fold_Change'] < 0]

            if stacked:
                up_annot = classify_subsets(upregulated[categorical_variable].tolist(), stacked)
                down_annot = classify_subsets(downregulated[categorical_variable].tolist(), stacked)
                dfs[comparison + "_UP"] = pd.DataFrame({
                    'Geneid': upregulated['Geneid'], categorical_variable: up_annot
                }).apply( sorted, axis=0).reset_index(drop=True)
                
                dfs[comparison + "_DOWN"] = pd.DataFrame({
                    'Geneid': downregulated['Geneid'], categorical_variable: down_annot
                }).apply( sorted, axis=0).reset_index(drop=True)
                
            else:
                dfs[comparison + "_UP"] = upregulated[['Geneid']].copy().apply(sorted, axis=0).reset_index(drop=True)
                dfs[comparison + "_DOWN"] = downregulated[['Geneid']].copy().apply(sorted, axis=0).reset_index(drop=True)

    return dfs

    
def upset_plots(
        input_files,
        outfile = "UpSet_genes.tsv",
        plot = "UpSet_plot",
        image_format = "svg",
        field = 1,
        de = True,
        stacked = False,
        categorical_variable = 'Annotation'
):

    # Unpack strings to be compared
    sets = []

    if de:
        if stacked:
            dfs = generate_lists_of_de_genes(filenames = input_files, stacked = stacked)
        else:
            dfs = generate_lists_of_de_genes(filenames = input_files)

        if dfs:
            
            # Unpack the contents of the dictionary of Pandas dataframes into list
            sets = {subsets: series['Geneid'].tolist() for subsets, series in dfs.items()}
            names = sets.keys()

        else:
            print("Problem with input dictionary of Pandas dataframes.")

    else:
        names = [re.sub(r"\\..*", "", os.path.splitext(os.path.basename(file_path))[0]) for file_path in input_files]
        sets = pd.read_csv(input_files, skip_blank_lines = True)

    # Get the unique elements in all subsets
    unique_elems = list(set().union(*sets.values()))

    # Identify presence/absence of each element in each subset
    df = pd.DataFrame([[
        e in st
        for comp, st in sets.items()
        for name in names
        if name == comp
    ]
        for e in unique_elems
    ],
        columns = names,  
        index = unique_elems
    )

    # Clean column names for readability
    garbage = re.compile('_geneids|.annotated.sorted|.tsv')
    df.columns = df.columns.str.replace(garbage, '', regex=True)
    df.columns = df.columns.str.replace('_', ' ', regex=True)

    # Save results of per gene presence in tab-separated file
    df.to_csv(outfile, sep = '\t', index = unique_elems, index_label = "GeneID")
    df.reset_index(inplace = True) # Move geneids from index to column
    df.rename(columns = {'index': 'Geneid'}, inplace = True)

    # Associate geneids in classified subsets with categorical variable from original parsing
    if stacked:
        concatenated_dfs = pd.concat(dfs.values()).drop_duplicates(subset='Geneid') 
        dfm = pd.merge(df, concatenated_dfs[['Geneid', str(categorical_variable)]], on='Geneid', how='left')
    else:
        dfm = df


    # Use subset combinations as index: expected input by UpSetPlot library
    subsets = [ col for col in df.columns if col != 'Geneid']
    dfm.set_index(subsets, inplace = True) 

    # Create UpSet plot with stacked bars
    plt.rc('font', family = 'Arial')
    upset = UpSet(dfm, orientation = 'horizontal', show_counts = True, intersection_plot_elements=0)
    upset.add_stacked_bars(by="Annotation", colors=cm.Pastel1, title="Intersections", elements=10)
    upset.plot()

    # Plot post-processing (MatPlotLib)
    ax = plt.gca()

    # Disable grid in all axes
    all_axes = plt.gcf().get_axes()
    for axis in all_axes:
        axis.grid(False)

    # Add bold ylabel
    ax.set_ylabel('Intersections', fontfamily = 'Arial', fontsize = 14, fontweight='bold')

    # Add horizontal line in x axis
    plt.axhline(y=0, color='black', linestyle='-', linewidth = 2)

    # Adjust plot margins
#    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.2, top=0.8, wspace=0.4, hspace=0.4)
    
    plt.savefig(plot, format = image_format, dpi = 600)
    plt.show()

    
if __name__  ==  "__main__":
    main()
