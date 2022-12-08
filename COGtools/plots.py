import matplotlib.pyplot as plt
import pandas as pd
import pkg_resources
from os import listdir
from os.path import isfile, join
import re


def categories_barplot(path_to_data, names=None, draft=False, cog_palette=True, include_unknown=True):
    """
    Visualizes the relative abundance of cog categories in the given genomes using barplots.
    :param path_to_data: the path to directory with genomes to be plotted
    :param names: names of the genomes, if not given, the file names in the folder will be used
    :type draft: bool
    :type cog_palette: bool
    :param draft: plotting draft genomes
    :param cog_palette: use palette from COG database
    :param include_unknown: include COGs unknown
    :return: relative abundance of cog categories in the given genomes in csv file
    :return: barplots of relative abundance of cog categories in the given genomes
    """
    header = "bacterium;J;A;K;L;B;D;Y;V;T;M;N;Z;W;U;O;X;C;G;E;F;H;I;P;Q;R;S;unknown\n"
    organisms = [f for f in listdir(path_to_data) if isfile(join(path_to_data, f)) and f.endswith(".txt")]
    if names is None:
        names = organisms

    my_string = header
    for index in range(len(organisms)):
        cat_dic = {"J": 0, "A": 0, "K": 0, "L": 0, "B": 0, "D": 0, "Y": 0, "V": 0, "T": 0, "M": 0, "N": 0,
                   "Z": 0, "W": 0, "U": 0, "O": 0, "X": 0, "C": 0, "G": 0, "E": 0, "F": 0, "H": 0, "I": 0, "P": 0,
                   "Q": 0, "R": 0, "S": 0, "-": 0}
        organism = organisms[index]
        organism_data = path_to_data + "/" + organism

        if draft:
            CogTools_data = pd.read_csv(organism_data, sep="\t", header=None, comment="#",
                                        names=("protein_id", "source", "cog", "cat"))

            for row in CogTools_data.index:
                cat = CogTools_data["cat"][row]
                cat_dic[cat] = cat_dic[cat] + 1

        else:
            CogTools_data = pd.read_csv(organism_data, sep="\t", header=None, comment="#",
                                        names=("seqname", "source", "type", "start",
                                               "end", "score", "strand", "frame",
                                               "attribute"), dtype=str)

            for row in CogTools_data.index:
                # get only useful information about each CDS: feature_id, name, COG, COG category
                type = CogTools_data["type"][row]
                if type == "CDS" or type == "pseudogene":
                    attribute = re.split(r"[=,;]", CogTools_data["attribute"][row])
                    cat_dic[attribute[5]] = cat_dic[attribute[5]] + 1

        if include_unknown is False:
            cat_dic["-"] = 0

        my_string = my_string + names[index].replace('\n','')
        number = sum(cat_dic.values())
        for category in cat_dic.keys():
            cat_dic[category] = cat_dic[category] * 100 / number
            my_string = my_string + ";" + str(cat_dic[category])

        my_string = my_string + "\n"

    with open(path_to_data + '/categories_barplot.csv', "w") as file_to_save:
        file_to_save.write(my_string)
        file_to_save.close()

    # plotting
    index = 1 if cog_palette else 3
    colors = []
    file = pkg_resources.resource_filename(__name__, 'COGtools-data/fun-20.tab.txt')
    features = (open(file, "r")).read()
    features = features.split('\n')

    for feature in features:
        feature = feature.split('\t')
        if "-" not in feature[0] and "RNA" not in feature[0]:
            colors.append('#' + feature[index])
        elif "-" in feature[0]:
            colors.append("#" + feature[1])

    data = pd.read_csv(path_to_data + '/categories_barplot.csv', sep=";")

    my_plot = data.plot(kind='bar', stacked=True, color=colors, legend=False)
    my_plot.set_xlabel("")
    my_plot.set_ylabel("Relative abundance", fontsize=20)

    plt.legend(title="Category", title_fontsize=15, fontsize=12, bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
    plt.yticks([0, 20, 40, 60, 80, 100], ["0%", "20%", "40%", "60%", "80%", "100%"], fontsize=20)
    plt.xticks([i for i in range(len(names))], names, fontstyle='italic', fontsize=18, rotation=0)
    plt.show()
