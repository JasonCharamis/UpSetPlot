 
import pandas as pd

# Parse command-line arguments
def parse_arguments() -> None:
    
    import argparse
    
    parser = argparse.ArgumentParser(description='Create Upset plots from a collection of lists.')
    
    parser.add_argument(
        '-de',
        '--de_files',
        type=str,
        required=False,
        help="List of de files to extract subsets of up-regulated and down-regulated genes."
    )
    parser.add_argument(
        '-ls',
        '--lists_of_strings',
        type=str,
        required=False,
        help="Lists of strings, in the case of not provided de files."
    )
    parser.add_argument(
        '-out',
        '--outfile',
        type=str,
        required=False,
        help="Output file to save results supporting the UpSet plot. Default: 'UpSet_genes.tsv' "
    )
    parser.add_argument(
        '-plt',
        '--plot',
        type=str,
        required=False,
        help="Output filename to print UpSet plot. Default: 'UpSet_plot.svg' "
    )
    parser.add_argument(
        '-img',
        '--image_format',
        required=False,
        type=str,
        help="Format to save image. Default: 'svg' "
    )
    parser.add_argument(
        '-vlist',
        '--variable_list',
        required=False,
        type=str,
        help="List of variables and categorical columns for stacking bars of UpSetPlot."
    )
    parser.add_argument(
        '-p_sum',
        '--pie_sum',
        required=False,
        type=str,
        help="Column to use for pie charts. Default: 'Orthology' "
    )
    parser.add_argument(
        '-fl',
        '--font_list',
        required=False,
        type=str,
        help="List of priority-ordered font names. Default: ['Aileron', 'Arial']."
    )
    
    args = parser.parse_args()

    if not any(vars(args).values()):
        parser.print_help()

    return args



# If a member is found in >1 categories, structure the 'stacked_categories' file as such,
# with the highest priority categories on top and the lowest priority categories at the end,
# to ensure that you keep the highest-priority category when dropping duplicate values.
# For example, I have an input file with 'Annotation' (categorical_variable) categories.
# To ensure that Uncharacterized will be the last assigned category (lowest priority),
# I structured it as following:
# ABC	ABC transporters
# cath	Cathepsins
# UDP	UGTs
# ase	Metabolic enzymes
# No hits	Uncharacterized
# hypothetical	Uncharacterized
# uncharacterized Uncharacterized

# The general idea is to go from as much specific to more broad/general categories (only the 1st column matters for the pattern matching and order identification)


def _import_stacking_variables(variable_list) -> dict:
    stacks = {}
    variables = {}

    with open(variable_list, 'r') as vlist:
        for lines in vlist:
            categories = lines.strip().split('\t')[0]
            stacks[categories] = {}

            files = lines.strip().split('\t')[1]

            with open(files, 'r') as category:
                for line in category.readlines():
                    columns = line.strip().split('\t')
                    string = columns[0]
                    group = columns[1]
                    stacks[categories][string] = group
    return stacks


# Classify subsets to stack bars
def classify_subsets(
        stacks: dict,
        column: pd.Series
) -> list:

    import re
    
    column_name = column.name
    column_values = column

    classified_subsets = []

    for stack, string in stacks.items():
        if re.search(column_name, stack, re.I):  # Match the column name with the stack name
            classified_stack = []
            for value in column_values:
                matched = False
                for string_, group in string.items():
                    if re.search(str(string_), str(value)):
                        classified_stack.append(group)
                        matched = True
                        break
                if not matched:
                    classified_stack.append('Other')

    return classified_stack



def generate_lists_of_genes(
        filenames,
        variable_list=None,
        de=True,
        pie_sum='Orthology'
) -> dict:

    import re
                            

    variables = {}
    
    if variable_list:
        variables = _import_stacking_variables(variable_list)
        
    dfs = {}

    with open(filenames, "r") as fls:
        fl = fls.readlines()
        
        for filename in fl:           
            filename = filename.strip('\n')
            filename_regex = re.compile(r"\\..\\/|\\\..\*|\..*")
            subset = re.sub(filename_regex, "", filename)

            if de:
                # Load pandas df with the contents of the de file and remove lines with Up- and Down-regulated flags
                de_file = pd.read_csv(filename, sep='\t', header=0, na_values=['NA', 'na'], keep_default_na=True)
                de_file['Geneid'] = de_file['Geneid'].astype(str)
                filtered_de_file = de_file[~de_file['Geneid'].str.contains('Upregulated|Downregulated')]

                # Isolate rows corresponding to up- and down-regulated genes
                filtered_de_file['Fold_Change'] = filtered_de_file['Fold_Change'].astype(float)             
                upregulated = pd.DataFrame(filtered_de_file[filtered_de_file['Fold_Change'] > 0])
                downregulated = pd.DataFrame(filtered_de_file[filtered_de_file['Fold_Change'] < 0])

                if variable_list:
                    up_dfs = []
                    down_dfs = []
                    dfs[subset + "_UP"] = {}
                    dfs[subset + "_DOWN"] = {}

                    for cat_var in variables.keys():
                        # Create DataFrame for upregulated genes with classification
                        up_df = pd.DataFrame({
                            'Geneid': upregulated['Geneid'].tolist(),
                            cat_var: classify_subsets(variables, upregulated[cat_var])
                        })
                        
                        # Create DataFrame for downregulated genes with classification
                        down_df = pd.DataFrame({
                            'Geneid': downregulated['Geneid'].tolist(),
                            cat_var: classify_subsets(variables, downregulated[cat_var])
                        })

                        # Append the DataFrames to the lists
                        up_dfs.append(up_df)
                        down_dfs.append(down_df)

                        # Merge all DataFrames in the lists on 'Geneid'
                    if up_dfs:
                        combined_up_df = upregulated[['Geneid']].copy()
                        
                        for up_df in up_dfs:
                            combined_up_df = pd.merge(
                                combined_up_df,
                                up_df,
                                on='Geneid',
                                how='left'
                            )

                        dfs[subset + "_UP"] = combined_up_df

                    if down_dfs:
                        combined_down_df = downregulated[['Geneid']].copy()
                    
                        for down_df in down_dfs:
                            combined_down_df = pd.merge(
                                combined_down_df,
                                down_df,
                                on='Geneid',
                                how='left'
                            )
                            
                        dfs[subset + "_DOWN"] = combined_down_df

                        
            elif pie_sum:
                dfs[subset + "_UP"] = pd.DataFrame({
                    'Geneid': upregulated['Geneid'],
                    pie_sum: upregulated[pie_sum].tolist()
                })

                dfs[subset + "_DOWN"] = pd.DataFrame({
                    'Geneid': downregulated['Geneid'],
                    pie_sum: downregulated[pie_sum].tolist()
                })
                
            else:
                dfs[subset + "_UP"] = pd.DataFrame({ 'Geneid': upregulated['Geneid'].copy() })
                dfs[subset + "_DOWN"] = pd.DataFrame({ 'Geneid': downregulated['Geneid'].copy() })

       
    return dfs



def _identify_fonts(font_list):

    """
    Identify and add fonts from a given list of font names to Matplotlib's font manager.
    
    Parameters:
    font_list (list): List of font names to identify and add.
    """

    from matplotlib import font_manager

    if isinstance(font_list, str):
        return font_list

    elif isinstance(font_list, list):
    
        # Find all system font files
        font_files = font_manager.findSystemFonts(fontpaths=None, fontext='otf')
    
        # Add specified fonts to the font manager
        for font_file in font_files:    
            for font_name in font_list:
                if font_name in font_file:
                    font_manager.fontManager.addfont(font_file)

            
def upset_plots(
        input_files,
        variable_list,
        de = True,
        outfile = "UpSet_genes.tsv",
        plot = "UpSet_plot",
        image_format = "svg",
        font_list = ['Aileron', 'Arial'],
        fontsize = 12
) -> None:

    import re
    from upsetplot import UpSet
    from matplotlib import pyplot as plt
    from matplotlib import cm
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib import font_manager
    from matplotlib.font_manager import FontProperties

    # Open input files
    if de and variable_list:
        dfs = generate_lists_of_genes(filenames = input_files, variable_list = variable_list, de = True)
        variables = _import_stacking_variables(variable_list)        

    elif de and not variable_list:
        dfs = generate_lists_of_genes(filenames=input_files, de=True)

    elif not de and variable_list:
        dfs = generate_lists_of_genes(filenames=input_files, de=False, stacked=first_categories)
        variables = _import_stacking_variables(variable_list)

    elif not de and not variable_list:
        dfs = generate_lists_of_genes(filenames=input_files, de=False)

        
    sets = []

    if dfs:
        sets = {subset: df[df.columns[0]].tolist() for subset, df in dfs.items()} # Unpack the contents of the dictionary of Pandas dataframes into list
        names = sets.keys()

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
        garbage = re.compile('_geneids|.annotated.sorted|P1.*.tsv')
        df.columns = df.columns.str.replace(garbage, '', regex=True)
        df.columns = df.columns.str.replace('_', ' ', regex=True)

        # Save results of per gene presence in tab-separated file
        df.to_csv(outfile, sep = '\t', index = unique_elems, index_label = 'GeneID' if de else 'Index')
        df.reset_index(inplace = True) # Move geneids from index to column
        df.rename(columns = {'index': 'Geneid' if de else 'Index'}, inplace = True)

        # Create UpSet plot depending on provided options
        # Associate geneids in classified subsets with categorical variable from original parsing
       
        if de:
            if variables:
                concatenated_dfs = pd.concat(dfs.values()).drop_duplicates(subset='Geneid')
                
                dfm = pd.merge(
                    df,
                    concatenated_dfs,
                    on = 'Geneid',
                    how = 'left'
                )

                subsets = [ col for col in df.columns if re.search ('UP|DOWN', col) ]
                dfm.set_index(subsets, inplace = True)
                upset = UpSet(dfm, orientation = 'horizontal', show_counts = True, intersection_plot_elements = 0)
                
                for stv in variables.keys():
                    upset.add_stacked_bars(by = stv, colors = cm.Accent, title = "Intersections", elements = 8)
                            
            else:
                concatenated_dfs = pd.concat(dfs.values()).drop_duplicates(subset='Geneid') 
                          
                dfm = pd.merge(
                    df,
                    concatenated_dfs[['Geneid']],
                    on = 'Geneid',
                    how = 'left'
                )
                
                subsets = [ col for col in df.columns if not re.search ('Geneid', col) ]
                dfm.set_index(subsets, inplace = True)
                upset = UpSet(dfm, orientation = 'horizontal', show_counts = True)

        else:
            dfm = df
            subsets = [ col for col in df.columns if not re.search ('Index', col) ]
            dfm.set_index(subsets, inplace = True)

            if stacked:
                upset = UpSet(dfm, orientation = 'horizontal', show_counts = True, intersection_plot_elements = 0)
                upset.add_stacked_bars(by = first_variable, colors = colormap, title = "Intersections", elements = 10)

            else:
                upset = UpSet(dfm, orientation = 'horizontal', show_counts = True)

        # Create the UpSet plot        
        _identify_fonts(font_list)
        plt.rcParams['font.family'] = font_list

        upset.plot()
        fig = plt.gcf()

        # Disable grid in all axes
        all_axes = fig.get_axes()
        for axis in all_axes:
            axis.grid(False)

        # Add bold ylabel
        ax = fig.gca()
        ax.set_ylabel('Intersections', fontfamily = font_list, fontsize = fontsize, fontweight='bold')
            
        # Add borders in UpSet plot bars
        for patch in ax.patches:
            patch.set_edgecolor('black')
            patch.set_linewidth(0.3)

        # Add horizontal line in x axis
        plt.axhline(y=-0.1, color='black', linestyle='-', linewidth = 1)
                
        plt.savefig(plot, format = image_format, dpi = 600)
        plt.show()

    else:
        print("Problem with input dictionary of Pandas dataframes.")



def pie_charts(
        input_files,
        column,
        plot,
        title,
        de = True,
        num_cols = 2,
        image_format = 'svg',
        font_list = [ 'Aileron', 'Arial']
) -> None:

    import re
    import numpy as np
    from matplotlib import pyplot as plt
    from matplotlib import cm
    from matplotlib.font_manager import FontProperties
    from adjustText import adjust_text
    
    if de:
        dfs = generate_lists_of_genes(filenames=input_files, pie_sum=column, de=True)
    else:
        dfs = generate_lists_of_genes(filenames=input_files, pie_sum=column, de=False)

    if dfs:
        num_subsets = len(dfs)
        num_rows = (num_subsets + num_cols - 1) // num_cols  # Calculate the number of rows needed

        # Create a single figure and set of axes
        fig, axes = plt.subplots(num_rows, num_cols, figsize = (20, 14))
        axes = axes.flatten()  # Flatten the axes array to simplify the loop

        # Consolidate all categories to define a common color map
        all_categories = sorted(set(cat for df in dfs.values() for cat in df[column].unique()))
        colors = plt.cm.tab20c(range(len(all_categories)))
        color_map = dict(zip(all_categories, colors))

        # Initialize a list to hold wedge handles for the legend
        wedge_handles = []

        for idx, (subset, df) in enumerate(dfs.items()):
            if idx >= len(axes):  # In case there are more subsets than axes available
                break

            # Count the occurrences of each 'Orthology' type
            counts = df[column].value_counts()

            # Adjust pie chart size based on the total number of elements 
            total_sum = counts.sum()
            radius = np.sqrt(total_sum) / 15
            
            # Plot the pie charts
            ax = axes[idx]
            
            wedges, texts, autotexts = ax.pie(
                counts,
                radius=radius,
                autopct='%1.1f%%',
                pctdistance=0.85,
                explode=[0.06] * len(counts),
                shadow=False,
                colors=[color_map[cat] for cat in counts.index],
                startangle=30,
                wedgeprops={'linewidth': 1, 'edgecolor': 'black'},
                labeldistance=1.1,  # Position labels slightly outside the pie
                textprops={'horizontalalignment': 'center', 'verticalalignment': 'center'}
            )

            # Adjust the labels
            for autotext, wedge in zip(autotexts, wedges):
                # Calculate the angle at the middle of the wedge
                angle = (wedge.theta2 + wedge.theta1) / 2
                
                # Rotate the label to align with the wedge
                autotext.set_rotation(angle * 180 / np.pi)
                
                # Set the font properties
                autotext.set_fontsize(8)
                autotext.set_fontfamily(font_list)

                x, y = np.cos(angle), np.sin(angle)
                    
                horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
                connectionstyle = "arc,rad=0"  # Set to straight line

                autotext.set_position((x * 1.05, y * 1.05))  # Adjust the radial distance as needed
                
                annot_texts = []
                
                annot = ax.annotate(
                    value,
                    xy=(1.02*x, 1.02*y),  # Coordinates at the edge of the wedge
                    xytext=(1.08* x, 1.08* y),  # Adjust the offset as needed
                    horizontalalignment=horizontalalignment,
                    bbox=dict(
                        boxstyle="round,pad=0.1",
                        edgecolor="white",
                        facecolor="white",
                        alpha=0
                    ),
                    
                arrowprops=dict(
                    arrowstyle="-",
                    connectionstyle=connectionstyle,
                ),
                    fontproperties=FontProperties(family=font_list, style='normal', size=8)
                )

                annot_texts.append(annot)
                
            adjust_text(annot_texts, ax = ax)

            cleaned_title = subset.split('.tsvf')[0].replace('_', ' ')
            ax.set_title(cleaned_title, fontdict = {'family': font_list, 'size': 14}, pad = 2)

        # Hide unused axes
        for i in range(idx + 1, len(axes)):
            axes[i].axis('off')

        # Create a common legend from the handles collected
        sorted_wedge_handles = sorted(wedge_handles, key=lambda x: (len(x[1].split()), x[1]))
        
        fig.legend(
            [handle for handle, label in sorted_wedge_handles], 
            [label for handle, label in sorted_wedge_handles],
            loc = 'center left',
            prop = FontProperties(family = font_list ,style = 'italic', size = 10),
            title = title,
            title_fontproperties = FontProperties(family = font_list, weight = 'bold', style = 'normal', size = 13),
            frameon = False
        )


        plt.savefig(plot, format = image_format, dpi = 600)
        plt.show()

    else:
        print("Problem with input dictionary of Pandas dataframes.")
        
        
def main() -> None:

    import os
    
    args = parse_arguments()

    if args.de_files and not args.lists_of_strings:        
        if not os.path.exists(args.de_files):
            print(f"Error: File '{args.de_files}' not found.")
            return

        if args.pie_sum:
            pie_charts(
                input_files = args.de_files.strip('\n'),
                plot = args.plot if args.plot else 'Pie_charts.svg',
                column = args.pie_sum if args.pie_sum else 'Orthology',
                image_format = args.image_format if args.image_format else 'svg',
                font_list = font_list if args.font_list else ['Aileron', 'Arial']
            )

        else:
            upset_plots(
                input_files = args.de_files.strip('\n'),
                variable_list = args.variable_list,
                plot = args.plot if args.plot else 'UpSet_plot.svg',
                image_format = args.image_format if args.image_format else 'svg',
                de = True,
                font_list = font_list if args.font_list else ['Aileron', 'Arial']
            )

            
    elif not args.de_files and args.lists_of_strings:
        if not os.path.exists(args.lists_of_strings):
            print(f"Error: File '{args.de_files}' not found.")
            return

        if pie_sum:
            pie_charts(
                input_files = args.lists_of_strings.strip('\n'),
                plot = args.plot if args.plot else 'Pie_charts.svg',
                column = args.pie_sum if args.pie_sum else 'Orthology',
                image_format = args.image_format if args.image_format else 'svg',
                font_list = font_list if args.font_list else ['Aileron', 'Arial']
            )
        
        else:
            upset_plots(
                input_files = args.lists_of_strings.strip('\n'),
                variable_list = args.variable_list,
                plot = args.plot if args.plot else 'UpSet_plot.svg',
                image_format = args.image_format if args.image_format else 'svg',
                de = False,
                font_list = font_list if args.font_list else ['Aileron', 'Arial']
            )        

    elif args.de_files and args.lists_of_strings:        
        print ( "Please select either de files or lists of strings.")

    else:
        print ("Please provide a list of files as input.")
        
    
if __name__  ==  "__main__":
    main()
