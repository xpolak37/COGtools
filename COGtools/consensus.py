import pandas as pd
from re import split, search
import os
from Bio import SeqIO
from random import randint
import pkg_resources


def read_file(file):
    """
    read processed file and get only location of feature and assigned COG
    """
    data = pd.read_csv(file, sep='\t', comment='#')
    new_data = data[["start", "end"]]
    new_data.insert(2, "COG", [search('COG=(.*);CAT', og).group(1) for og in data["attribute"]], True)
    return [data, new_data]


def consensus(organism_name, em_file=None, om_file=None, batch_file=None, fasta_file=None, get_pseudo=False,
              get_ncrna=False, gff_file=None, cat_choice=1, output_dir=os.getcwd()):
    """
    Improves the functional annotation of the bacterial genome using a consensus of three programs:
    eggNOG-mapper, Operon-mapper and Batch CD-Search. Function saves all predicted features and COG assignments
    and prepares the outputs file for visualization with DNAPlotter
    :type organism_name: str
    :type cat_choice: int
    :param organism_name: organism name
    :param em_file: the path to Eggnog-mapper processed file
    :param om_file: the path to Operon-mapper processed file
    :param batch_file: the path to Batch CD-Search processed file
    :param fasta_file: the path to genomic sequence
    :param output_dir: the output file
    :param cat_choice: select the option to assign a category (1-4)
    :type get_pseudo: bool
    :type get_ncrna: bool
    :param gff_file: the path to gff file where all features are stored
    :return:  file with functional annotation of the bacterial genome
    """
    # how many files are given
    nones = [em_file, om_file, batch_file]
    if nones.count(None) == 2:
        for tool in nones:
            if tool is not None:
                df = pd.read_csv(tool, sep='\t', comment='#')
                df = categories_choice(df, cat_choice=cat_choice)

    elif nones.count(None) == 0:
        df = create_consensus(em_file, om_file, batch_file)
        df = categories_choice(df, cat_choice=cat_choice)

    else:
        print("Three files are needed to create consensus.")

    # add pseudogenes and/or ncRNA
    if get_pseudo or get_ncrna:
        df = get_features(gff_file, df, get_pseudo, get_ncrna)

    # save the created dataframe into new file and add genomic sequence
    f = open(output_dir + '/' + organism_name + '_file_to_plot.txt', 'w')
    f.write('# created with COGtools 1.0.0\n'
            '# AC number: ' + df["seqname"][0] + '\n'
                                                 '# COG annotation\n')
    f.close()
    df.to_csv(output_dir + '/' + organism_name + '_file_to_plot.txt', sep='\t', index=False, header=False, mode='a')

    if fasta_file is not None:
        fasta_data = (open(fasta_file, "r")).read()
        with open(output_dir + '/' + organism_name + '_file_to_plot.txt', 'a') as my_file:
            my_file.write('\n' + fasta_data)
        my_file.close()


def get_features(gff_file, df, get_pseudo, get_ncrna):
    """
    change the feature type to a pseudogene according to information in gff_file
    and add ncRNA feature to the dataframe
    """
    gff_file = pd.read_csv(gff_file, sep="\t", header=None, comment="#", names=("seqname", "source", "type", "start",
                                                                                "end", "score", "strand", "frame",
                                                                                "attribute"))
    if get_pseudo:
        pseudogenes = gff_file.loc[gff_file['type'] == 'pseudogene']
        for row in pseudogenes.index:
            start = pseudogenes.start[row]
            df.loc[df['start'] == start, 'type'] = 'pseudogene'
    if get_ncrna:
        ncrnas = gff_file.loc[gff_file['type'] == 'ncRNA']
        df = pd.concat([df, ncrnas], ignore_index=True)

    return df


def consensus_draft(organism_name, proteins=None, em_file=None, om_file=None, batch_file=None, cat_choice=1,
                    output_dir=os.getcwd()):
    """
        Improves the functional annotation of the draft bacterial genome using a consensus of three programs:
        eggNOG-mapper, Operon-mapper and Batch CD-Search.
        :type organism_name: str
        :type cat_choice: int
        :param organism_name: organism name
        :param cat_choice: select the option to assign a category (1-4)
        :param proteins: the path to downloaded proteins of the draft genome (fasta)
        :param em_file: the path to Eggnog-mapper processed file
        :param om_file: the path to Operon-mapper processed file
        :param batch_file: the path to Batch CD-Search processed file
        :param output_dir: output file
        :return:  file with COG assignments
        """
    nones = [em_file, om_file, batch_file]
    if nones.count(None) == 2:
        for tool in nones:
            if tool is not None:
                df = pd.read_csv(tool, sep='\t', comment='#')
                df = categories_choice_draft(df, cat_choice=cat_choice)

    elif nones.count(None) == 0:
        df = create_consensus_draft(proteins, em_file, om_file, batch_file)
        df = categories_choice_draft(df, cat_choice=cat_choice)

    else:
        print("Three files are needed to create consensus.")

    # save the created dataframe into new file and add genomic sequence
    f = open(output_dir + '/consensus_' + organism_name + '.txt', 'w')
    f.write('# created with COGtools 1.0.0\n# AC number: unknown\n# COG annotation\n')
    f.close()
    return df.to_csv(output_dir + '/consensus_' + organism_name + '.txt', sep='\t', index=False, mode='a')


def categories_choice(df, cat_choice=1):
    # categories choice
    if cat_choice == 1:
        df['attribute'] = df['attribute'].apply(
            lambda x: x.replace(search('CAT=(.*);', x).group(0), "CAT=" + search('CAT=(.*);', x).group(0)[4] + ";")
            if len(search('CAT=(.*);', x).group(1)) > 1 else x)
    elif cat_choice == 2:
        df['attribute'] = df['attribute'].apply(
            lambda x: x.replace(search('CAT=(.*);', x).group(0), "CAT=" + search('CAT=(.*);', x).group(0)
                                [randint(0, len(search('CAT=(.*);', x).group(1)) - 1)] + ";")
            if len(search('CAT=(.*);', x).group(1)) > 1 else x)
    else:
        cat_dic = {"-": 0, "J": 0, "A": 0, "K": 0, "L": 0, "B": 0, "D": 0, "Y": 0, "V": 0, "T": 0, "M": 0, "N": 0,
                   "Z": 0, "W": 0, "U": 0, "O": 0, "X": 0, "C": 0, "G": 0, "E": 0, "F": 0, "H": 0, "I": 0, "P": 0,
                   "Q": 0, "R": 0, "S": 0}

        for row in df.index:
            # get only useful information about each CDS: feature_id, name, COG, COG category
            type = df["type"][row]
            if type == "CDS":
                attribute = split(r"[=,;]", df["attribute"][row])
                if len(attribute[5]) == 1:
                    cat_dic[attribute[5]] = cat_dic[attribute[5]] + 1
                else:
                    for i in range(len(attribute[5])):
                        cat_dic[attribute[5][i]] = cat_dic[attribute[5][i]] + 1

        if cat_choice == 3:
            df['attribute'] = df['attribute'].apply(
                lambda x: x.replace(
                    search('CAT=(.*);', x).group(0), "CAT=" + max(
                        {search('CAT=(.*);', x).group(1)[i]: cat_dic[search('CAT=(.*);', x).group(1)[i]]
                         for i in range(len(search('CAT=(.*);', x).group(1)))},
                        key={search('CAT=(.*);', x).group(1)[i]: cat_dic[search('CAT=(.*);', x).group(1)[i]]
                             for i in range(len(search('CAT=(.*);', x).group(1)))}.get) + ";")
                if len(search('CAT=(.*);', x).group(1)) > 1 else x)

        elif cat_choice == 4:
            df['attribute'] = df['attribute'].apply(
                lambda x: x.replace(
                    search('CAT=(.*);', x).group(0), "CAT=" + min(
                        {search('CAT=(.*);', x).group(1)[i]: cat_dic[search('CAT=(.*);', x).group(1)[i]]
                         for i in range(len(search('CAT=(.*);', x).group(1)))},
                        key={search('CAT=(.*);', x).group(1)[i]: cat_dic[search('CAT=(.*);', x).group(1)[i]]
                             for i in range(len(search('CAT=(.*);', x).group(1)))}.get) + ";")
                if len(search('CAT=(.*);', x).group(1)) > 1 else x)
    return df


def categories_choice_draft(df, cat_choice=1):
    # categories choice
    if cat_choice == 1:
        df['cat'] = df['cat'].apply(lambda x: x[0])  # the first category
    elif cat_choice == 2:
        df['cat'] = df['cat'].apply(lambda x: x[randint(0, len(x) - 1)] if len(x) > 1 else x)  # random category
    else:
        # category count
        cat_dic = {"-": 0, "J": 0, "A": 0, "K": 0, "L": 0, "B": 0, "D": 0, "Y": 0, "V": 0, "T": 0, "M": 0, "N": 0,
                   "Z": 0, "W": 0, "U": 0, "O": 0, "X": 0, "C": 0, "G": 0, "E": 0, "F": 0, "H": 0, "I": 0, "P": 0,
                   "Q": 0, "R": 0, "S": 0}

        for row in df.index:
            cat = df["cat"][row]
            if len(cat) == 1:
                cat_dic[cat] = cat_dic[cat] + 1
            else:
                for i in range(len(cat)):
                    cat_dic[cat[i]] = cat_dic[cat[i]] + 1

        if cat_choice == 3:  # the most numerous category
            df['cat'] = df['cat'].apply(lambda x: max({x[i]: cat_dic[x[i]] for i in range(len(x))},
                                                      key={x[i]: cat_dic[x[i]] for i in range(len(x))}.get)
                        if len(x) > 1 else x)

        elif cat_choice == 4:  # the least numerous category
            df['cat'] = df['cat'].apply(lambda x: min({x[i]: cat_dic[x[i]] for i in range(len(x))},
                                                      key={x[i]: cat_dic[x[i]] for i in range(len(x))}.get)
                        if len(x) > 1 else x)

    return df


def create_consensus(em_file, om_file, batch_file):
    # read eggnog_mapper file
    [em_data, new_em_data] = read_file(em_file)
    # read operon_mapper file
    [om_data, new_om_data] = read_file(om_file)
    # read batch cd search file
    [batch_data, new_batch_data] = read_file(batch_file)

    # merge the three dataframes to save all predicted features into one dataframe
    # cog_x = eggnog-mapper #cog_y = operon-operon-mapper #cog = batch cd-search
    new_df = pd.merge(new_em_data, new_om_data, on="start", how="outer")
    new_df = (pd.merge(new_df, new_batch_data, on="start", how="outer")).fillna("-")

    # create new DataFrame to store the necessary information
    df = pd.DataFrame(columns=["seqname", "source", "type", "start", "end", "score", "strand", "frame", "attribute"])

    # data from COG database
    cogs_file = pkg_resources.resource_filename(__name__, 'COGtools-data/cogs.txt')
    cogs_data = (open(cogs_file, "r")).readlines()

    # iterate through features
    for row in new_df.index:
        # get assigned COGs, number of NaN and start of the feature
        cogs = [new_df["COG_x"][row], new_df["COG_y"][row], new_df["COG"][row]]
        nan = cogs.count("-")
        start = new_df["start"][row]
        # find out which tools match in the COG assignment
        idxs = [[cogs[:idx].index(item), idx] for idx, item in enumerate(cogs) if
                item in cogs[:idx]]

        # if all three tools have assigned the COG
        if nan == 0:
            # if all three tools match in assignment, add Batch CD-Search
            if len(idxs) != 1:
                df = pd.concat([df, batch_data.loc[batch_data.start == start, :]], ignore_index=True)

            # two tools match in assignment
            else:
                #  if batch and eggnog or operon -> add batch, if operon and eggnog -> add eggnog
                if idxs[0] == [1, 2] or idxs[0] == [0, 2]:
                    df = pd.concat([df, batch_data.loc[batch_data.start == start, :]], ignore_index=True)
                else:
                    df = pd.concat([df, em_data.loc[em_data.start == start, :]], ignore_index=True)
                    cog = search('COG=(.*);CAT', df["attribute"][row]).group(1)
                    # if cog is from COG database
                    try:
                        index = [i for i, s in enumerate(cogs_data) if cog in s][0]
                        cat = cogs_data[index][8:-1]
                        df["attribute"][row] = df["attribute"][row].replace(
                            search('CAT=(.*);', df["attribute"][row]).group(0), "CAT=" + cat + ";")

                    # else - it is from eggNOG
                    except IndexError:
                        continue

        # if only one tool has assigned the COG
        elif nan == 2:
            # find out which one was it and add it
            which = [number for number in [0, 1, 2] if number not in idxs[0]][0]
            programs = ['em_data', 'om_data', 'batch_data']
            df = pd.concat([df, (locals()[programs[which]]).loc[(locals()[programs[which]]).start == start, :]],
                           ignore_index=True)
            if which == 0 or which == 1:
                cog = search('COG=(.*);CAT', df["attribute"][row]).group(1)
                # if cog is from COG database
                try:
                    index = [i for i, s in enumerate(cogs_data) if cog in s][0]
                    cat = cogs_data[index][8:-1]
                    df["attribute"][row] = df["attribute"][row].replace(
                        search('CAT=(.*);', df["attribute"][row]).group(0), "CAT=" + cat + ";")

                # else - it is from eggNOG or ROG
                except IndexError:
                    continue

        # if two tools have assigned the COG
        elif nan == 1:
            nan_position = cogs.index("-")
            # add batch
            if nan_position == 0 or nan_position == 1:
                df = pd.concat([df, batch_data.loc[batch_data.start == start, :]], ignore_index=True)
            # add eggnog
            else:
                df = pd.concat([df, em_data.loc[em_data.start == start, :]], ignore_index=True)
                cog = search('COG=(.*);CAT', df["attribute"][row]).group(1)
                # if cog is from COG database
                try:
                    index = [i for i, s in enumerate(cogs_data) if cog in s][0]
                    cat = cogs_data[index][8:-1]
                    df["attribute"][row] = df["attribute"][row].replace(
                        search('CAT=(.*);', df["attribute"][row]).group(0), "CAT=" + cat + ";")

                # else - it is from eggNOG or ROG
                except IndexError:
                    continue

        # no tool has assigned the COG, either it is not assigned or it is not a CDS -> add operon-mapper
        else:
            df = pd.concat([df, om_data.loc[om_data.start == start, :]], ignore_index=True)

    return df


def create_consensus_draft(proteins, em_file, om_file, batch_file):
    # get id of downloaded proteins
    proteins = list(SeqIO.parse(proteins, "fasta"))
    id_downloaded = [i.id for i in proteins]

    # processed files
    em_data = pd.read_csv(em_file, sep='\t', comment='#')
    om_data = pd.read_csv(om_file, sep='\t', comment='#')
    batch_data = pd.read_csv(batch_file, sep='\t', comment='#')

    # get only tables with protein_id and cog from processed files and merged them
    new_em_data = em_data[["protein_id", "cog"]]
    new_om_data = om_data[["protein_id", "cog"]]
    new_batch_data = batch_data[["protein_id", "cog"]]
    new_downloaded_data = pd.DataFrame(id_downloaded, columns=['protein_id'])

    new_df = pd.merge(new_em_data, new_om_data, on="protein_id", how="outer")
    new_df = (pd.merge(new_df, new_batch_data, on="protein_id", how="outer")).fillna("-")
    new_df = (pd.merge(new_df, new_downloaded_data, on="protein_id", how="outer")).fillna("-")

    # table for final data
    df = pd.DataFrame(columns=["protein_id", "source", "cog", "cat"])

    # data from COG database
    cogs_file = pkg_resources.resource_filename(__name__, 'COGtools-data/cogs.txt')
    cogs_data = (open(cogs_file, "r")).readlines()

    for row in new_df.index:
        # get assigned COGs, number of NaN and start of the feature
        cogs = [new_df["cog_x"][row], new_df["cog_y"][row], new_df["cog"][row]]
        nan = cogs.count("-")
        protein_id = new_df["protein_id"][row]
        # find out which tools match in the COG assignment
        idxs = [[cogs[:idx].index(item), idx] for idx, item in enumerate(cogs) if
                item in cogs[:idx]]

        # if all three tools have assigned the COG
        if nan == 0:
            # if all three tools match in assignment, add eggnog-mapper
            if len(idxs) != 1:
                df = pd.concat([df, batch_data.loc[batch_data.protein_id == protein_id, :]], ignore_index=True)
            # two tools match in assignment
            else:
                #  if batch and eggnog or operon -> add batch, if operon and eggnog -> add eggnog
                if idxs[0] == [1, 2] or idxs[0] == [0, 2]:
                    df = pd.concat([df, batch_data.loc[batch_data.protein_id == protein_id, :]], ignore_index=True)
                else:
                    df = pd.concat([df, em_data.loc[em_data.protein_id == protein_id, :]], ignore_index=True)
                    cog = df["cog"][row]
                    # if cog is from COG database
                    try:
                        index = [i for i, s in enumerate(cogs_data) if cog in s][0]
                        df["cat"][row] = cogs_data[index][8:-1]

                    # else - it is from eggNOG
                    except IndexError:
                        continue

        # if only one tool has assigned the COG
        elif nan == 2:
            # find out which one was it and add it
            which = [number for number in [0, 1, 2] if number not in idxs[0]][0]
            programs = ['em_data', 'om_data', 'batch_data']
            df = pd.concat([df, (locals()[programs[which]]).loc[(locals()[programs[which]]).protein_id == protein_id, :]
                            ], ignore_index=True)

            if which == 0 or which == 1:
                cog = df["cog"][row]
                # if cog is from COG database
                try:
                    index = [i for i, s in enumerate(cogs_data) if cog in s][0]
                    df["cat"][row] = cogs_data[index][8:-1]

                # else - it is from eggNOG or ROG
                except IndexError:
                    continue

        # if two tools have assigned the COG
        elif nan == 1:
            nan_position = cogs.index("-")
            # add batch cd search of eggnog_mapper
            if nan_position == 0 or nan_position == 1:
                df = pd.concat([df, batch_data.loc[batch_data.protein_id == protein_id, :]], ignore_index=True)
            else:
                df = pd.concat([df, em_data.loc[em_data.protein_id == protein_id, :]], ignore_index=True)
                cog = df["cog"][row]
                # if cog is from COG database
                try:
                    index = [i for i, s in enumerate(cogs_data) if cog in s][0]
                    df["cat"][row] = cogs_data[index][8:-1]

                # else - it is from eggNOG or ROG
                except IndexError:
                    continue

        # no tool has assigned the COG, either it is not assigned or it is not a CDS -> add operon-mapper
        else:
            new_row = pd.DataFrame({"protein_id": [protein_id], "source": ["-"], "cog": ["-"], "cat": ["-"]})
            df = pd.concat([df, new_row], ignore_index=True)

    return df
