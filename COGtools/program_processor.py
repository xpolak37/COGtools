from re import search, split
import pandas as pd
from Bio import SeqIO
import os
import pkg_resources
from Bio import Align
import warnings


def em_processor(organism_name, em_file, cds_file, cogs_only=False, output_dir=os.getcwd()):
    """
    Process the output file (decorated.gff) from eggNOG-mapper tool into more structured COGtools-data.
    The outputs of this function is file in gff format that contains a suitable header with information about CDSs with
    assigned COG by eggNOG-mapper
    :type organism_name: str
    :type cogs_only: bool
    :param em_file: the path to eggNOG-mapper output file
    :param cogs_only: include only COGs, orthologous groups from eggNOG will be considered as unknown
    :param cds_file: the path to eggNOG-mapper input file
    :param output_dir: the output directory
    :return: processed file
    """
    em_data = pd.read_csv(em_file, sep="\t", header=None, comment="#", names=("seqname", "source", "type", "start",
                                                                              "end", "score", "strand", "frame",
                                                                              "attribute"))
    cds_data = SeqIO.parse(cds_file, "fasta")
    # get only header of each CDS
    desc = [record.description for record in cds_data]
    # iterate through CDSs header

    for record in desc:
        seq_id = record[0:record.index(" [")]
        # include all possible locations forward/reverse strand and join
        dic = {'pattern': "complement\((.*?)\)", 'strand': "-"} if "complement" in record else \
            {'pattern': "location=(.*?)]", 'strand': "+"}
        location = search(dic['pattern'], record).group(1)

        if "<" in location:
            location = location[1:]

        if "join" in location:
            location = location[:-1] if location[-1] == ")" else location
            location = (location[5:].split(".."))
            location = [location[0], location[len(location) - 1]]
        else:
            location = location.split("..")

        if ">" in location[1]:
            location[1] = location[1][1:]

        # add the information about location to the corresponding row of eggnog-mappers outputs
        em_data.loc[em_data.seqname == seq_id, ["seqname", "strand", "start", "end"]] = \
            [seq_id[seq_id.index("|") + 1:seq_id.index("_prot")], dic['strand'], location[0], location[1]]

    for row in em_data.index:
        # get only useful information about each CDS: feature_id, name, COG, COG category
        attribute = em_data["attribute"][row]
        try:
            where = search(r'Bacteria', attribute).span(0)[0]
            dic = {'COG': attribute[where - 10: where - 3]}
            if "COG" in dic['COG']:
                attribute = split(r"[=,;]", attribute)
                dic['seq_id'] = attribute[attribute.index("ID") + 1][4:]
                dic['cat'] = attribute[attribute.index("em_COG_cat") + 1]
                dic['desc']= attribute[attribute.index("em_desc") + 1]
            else:
                if cogs_only:
                    em_data = em_data.drop(row)
                    continue
                else:
                    dic['COG'] = attribute[where-8:where-3]
                    attribute = split(r"[=,;]", attribute)
                    dic['seq_id'] = attribute[attribute.index("ID") + 1][4:]
                    dic['cat'] = attribute[attribute.index("em_COG_cat") + 1]
                    dic['desc'] = attribute[attribute.index("em_desc") + 1]
        except:
            # no bacteria group in annotation
            em_data = em_data.drop(row)
            continue

        dic['cat'] = 'S' if dic['cat'] == 'None' else dic['cat']

        em_data.loc[row, "attribute"] = "".join(["ID=", dic['seq_id'], ";COG=", dic['COG'], ";CAT=", dic['cat'],
                                                 ";desc=", dic['desc']])

    f = open(output_dir + '/em_' + organism_name + '.gff', 'w')
    f.write('# created with COGtools 1.0.0\n# AC number: ' + em_data["seqname"][0] + "\n# Processed data from eggNOG-mapper\n")
    f.close()
    return em_data.to_csv(output_dir + '/em_' + organism_name + '.gff', sep='\t', index=False, mode = "a")


def em_processor_draft(organism_name, em_file, cogs_only=False, output_dir=os.getcwd()):
    """
      Process the output file (decorated.gff) from eggNOG-mapper tool into more structured COGtools-data.
      The output of this function is file in txt format that contains a suitable header with assigned cogs
      and their categories.
      :type organism_name: str
      :type cogs_only: bool
      :param em_file: the path to eggNOG-mapper output file
      :param cogs_only: neglect other orhologous groups than COGs
      :param output_dir: the output directory
      :return: processed file
      """
    # data needed: annotated file by eggNOG-mapper
    em_data = pd.read_csv(em_file, sep="\t", header=None, comment="#", names=("seqname", "source", "type", "start",
                                                                              "end", "score", "strand", "frame",
                                                                              "attribute"))
    # table for processed data
    em_table = pd.DataFrame(columns=["protein_id","source","cog","cat"])

    # cog categories in COG database
    cogs_file = pkg_resources.resource_filename(__name__, 'COGtools-data/cogs.txt')
    cat_data = (open(cogs_file, "r")).readlines()

    for row in em_data.index:
        # get only useful information about each CDS: feature_id, COG, COG category
        attribute = em_data["attribute"][row]

        try:
            where = search(r'Bacteria', attribute).span(0)[0]
            COG = attribute[where-10: where-3]
            if "COG" in COG:
                attribute = split(r"[=,;]", attribute)
                seq_id = attribute[attribute.index("ID") + 1]
                CAT = attribute[attribute.index("em_COG_cat") + 1]
            else:
                if cogs_only:
                    continue
                else:
                    COG = attribute[where-8:where-3]
                    attribute = split(r"[=,;]", attribute)
                    seq_id = attribute[attribute.index("ID") + 1]
                    CAT = attribute[attribute.index("em_COG_cat") + 1]

        except AttributeError:
            continue

        CAT = 'S' if CAT == 'None' else CAT

        # save this info into the table
        new_row = pd.DataFrame(
            {"protein_id": [seq_id], "source": ["eggnog_mapper"], "cog": [COG], "cat": [CAT]})
        em_table = pd.concat([em_table, new_row], ignore_index=True)

    f = open(output_dir + '/em_' + organism_name + '.txt', 'w')
    f.write('# created with COGtools 1.0.0\n# AC number: unknown\n# Processed data from eggNOG-mapper\n')
    f.close()
    return em_table.to_csv(output_dir + '/em_' + organism_name + '.txt', sep='\t', index=False, mode='a')


def om_processor(organism_name, orf_file, cog_file, output_dir=os.getcwd()):
    """
    Process the outputs files (ORF_coordinates.txt and predicted_COGs.txt) from Operon-mapper into more structured COGtools-data.
    The outputs of this function is file in gff format that contains a suitable header with information about all
    predicted features
    :type organism_name: str
    :param orf_file: the path to Operon-mapper outputs file ORFs_coordinates.txt
    :param cog_file: the path to Operon-mapper outputs file predicted_COGs.txt
    :param output_dir: the output directory
    :return: processed file
    """
    orf_data = pd.read_csv(orf_file, sep="\t", header=None, comment="#", names=("seqname", "source", "type", "start",
                                                                                "end", "score", "strand", "frame",
                                                                                "attribute"))
    cog_data = pd.read_csv(cog_file, sep="\t", header=None, comment="#", names=("ID", "COG", "category"))

    for row in orf_data.index:
        # iterate through all features in ORF file and save the relevant information from the COG file
        attribute = split(r"[=,;]", orf_data["attribute"][row])
        feature_id = attribute[attribute.index("ID") + 1]
        try:
            cog = cog_data.loc[cog_data["ID"] == feature_id, "COG"].values[0]
            CATS = cog_data.loc[cog_data["ID"] == feature_id, "category"].values[0]
            CAT = search('\w+', CATS).group(0)
            dic = {'cat': CAT, 'desc': cog_data.loc[cog_data["ID"] == feature_id, "category"].values[0].split("] ")[1]}

            orf_data.loc[row, "attribute"] = "".join(
                ["ID=", feature_id, ";COG=", cog, ";CAT=", dic['cat'], ";desc=", dic['desc']])

        except:
            orf_data.loc[row, "attribute"] = "".join(["ID=", feature_id, ";COG=", "-", ";CAT=", "-", ";desc=", "-"])

    f = open(output_dir + '/om_' + organism_name + '.gff', 'w')
    f.write('# created with COGtools 1.0.0\n# AC number: ' + orf_data["seqname"][0] +
            "\n# Processed data from Operon-mapper\n")
    f.close()
    return orf_data.to_csv(output_dir + '/om_' + organism_name + '.gff', sep='\t', index=False, mode = "a")
    # return orf_data.to_csv(output_dir + '/om_' + organism_name + '.gff', sep='\t', index=False)


def batch_splitter(organism_name, gene_file, output_dir=os.getcwd()):
    """
    :type organism_name: str
    :param gene_file: the path to file to be split
    :param output_dir: the output directory
    :return: two split files
    """
    records = list(SeqIO.parse(gene_file, "fasta"))
    if len(records) > 1000:
        num_of_parts = int(len(records) / 1000) + (len(records) % 1000 > 0)
        start = 0
        end = 999
        for i in range(num_of_parts):
            SeqIO.write(records[start:end], output_dir + "/" + organism_name + "_" + str(i) + ".fasta", "fasta")
            start = end
            end = end + 999


def batch_merger(organism_name, files, output_dir=os.getcwd()):
    """
    :type organism_name: str
    :param files: paths to annotated files in list
    :param output_dir: the output directory
    :return: a merged file
    """
    data1 = open(files[0]).read()
    data1 = data1[data1.index("Q#"):len(data1)]
    for i in range(1, len(files)):
        data = open(files[i]).read()
        data = data[data.index("Q#"):len(data)]
        data1 += data

    with open(output_dir + "/" + organism_name + "_merged_hitdata.txt", "w") as file:
        file.write(data1)


def batch_processor(organism_name, batch_file, output_dir=os.getcwd()):
    """
    Process the outputs file (hitdata.txt) from Batch CD-Search tool into more structured COGtools-data.
    The outputs of this function is file in gff format that contains a suitable header with information about CDSs with
    assigned COG by Batch CD-Search
    :type organism_name: str
    :param batch_file: the path to Batch CD-Search outputs file hitdata.txt
    :param output_dir: the output file
    :return: processed file
    """

    # create new DataFrame to store the necessary information
    batch_gff = pd.DataFrame(
        columns=["seqname", "source", "type", "start", "end", "score", "strand", "frame", "attribute"])
    batch_data = (open(batch_file).read())
    batch_data = (batch_data[batch_data.index("Q#"):len(batch_data) - 1]).split('\n')
    query = ''
    cogs_file = pkg_resources.resource_filename(__name__, 'COGtools-data/cogs.txt')
    cogs_data = (open(cogs_file, "r")).readlines()

    # iterate through queries
    for row in batch_data:
        new_query = search('Q#\d+', row).group(0)

        # if COG is not assigned or new_line is duplicate, continue to next iteration
        if not ('\tspecific\t' in row) or query == new_query:
            continue

        query = new_query
        seq_id = row[row.index("|") + 1:row.index("_prot")]
        id = row[row.index("|") + 1:row.index(" [")]
        id = "".join(["ID=", id])

        # include all possible locations forward/reverse strand and join
        dic = {'pattern': "complement\((.*?)\)", 'strand': '-'} if 'complement' in row else \
            {'pattern': "location=(.*?)]", 'strand': '+'}
        location = search(dic['pattern'], row).group(1)

        if "<" in location:
            location = location[1:]

        if "join" in location:
            location = location[:-1] if location[-1] == ")" else location
            location = (location[5:].split(".."))
            location = [location[0], location[len(location) - 1]]

        else:
            location = location.split("..")
            location = [location[0], location[1]]

        if ">" in location[1]:
            location[1] = location[1][1:]

        try:
            COG = "".join(["COG=", search("(COG\d+)", row).group(1)])

        except AttributeError:
            continue

        # update in COG 2021
        if COG == "COG=COG3512":
            COG = "COG=COG1343"

        # adding categories
        index = [i for i, s in enumerate(cogs_data) if COG[4:] in s][0]
        CAT = search('\t\w+', cogs_data[index]).group(0)[1:]
        CAT = "".join(["CAT=", CAT])

        attribute = id + ";" + COG + ";" + CAT + ";"
        new_row = pd.DataFrame(
            {"seqname": [seq_id], "source": ["batch-cd-search"], "type": ["CDS"], "start": [location[0]], "end": [location[1]],
             "score": ["."], "strand": [dic["strand"]], "frame": ["0"], "attribute": [attribute]})

        batch_gff = pd.concat([batch_gff, new_row], ignore_index=True)

    f = open(output_dir + '/batch_' + organism_name + '.gff', 'w')
    f.write('# created with COGtools 1.0.0\n# AC number: ' + seq_id +
            "\n# Processed data from Batch CD-Search\n")
    f.close()
    return batch_gff.to_csv(output_dir + '/batch_' + organism_name + '.gff', sep='\t', index=False, mode = "a")


def om_processor_draft(organism_name, proteins, operon_proteins, operon_cogs, gff_included=True, output_dir=os.getcwd()):
    """
    Process the outputs files (predicted_protein_sequences.txt and predicted_COGs.txt) from Operon-mapper into more
    structured COGtools-data. The output of this function is a file in txt format that contains a suitable header with
    predicted cogs and categories
    :type organism_name: str
    :type gff_included: bool
    :param proteins: the path to downloaded proteins
    :param operon_proteins: the path to Operon-mapper output file predicted_protein_sequences.txt
    :param operon_cogs: the path to Operon-mapper output file predicted_COGs.txt
    :param gff_included: a gff file was used in the Operon-mapper
    :param output_dir: the output directory
    :return: processed file
    """
    # data needed: downloaded proteins, proteins predicted by Operon-mapper, COGs prediction by Operon-mapper
    proteins = list(SeqIO.parse(proteins, "fasta"))
    operon_proteins = list(SeqIO.parse(operon_proteins, "fasta"))
    cog_data = pd.read_csv(operon_cogs, sep="\t", header=None, comment="#", names=("ID", "COG", "category"))
    sequences = [i.seq for i in proteins]

    # table for saving processed data
    operon_table = pd.DataFrame(columns=["protein_id", "source", "cog", "cat"])

    if gff_included is False:
        #
        warnings.warn('The analyzed loci do not match and will be evaluated based on sequence alignment')
        aligner = Align.PairwiseAligner(mode='local')

    last_index = -1
    # iterate through all downloaded proteins
    for i in range(len(sequences)):
        found = False
        my_protein = sequences[i]

        # search for similar protein predicted by Operon-mapper
        for j in range(1, 10):
            if last_index+j < len(operon_proteins):
                protein_operon = operon_proteins[last_index+j].seq[:-1]

                if gff_included:
                    if my_protein == protein_operon:
                        found = True

                else:
                    # calculate alignment score
                    score = aligner.score(my_protein, protein_operon)

                    lengthDownloaded = len(my_protein)
                    lengthOperon = len(protein_operon)

                    # if proteins are almost identical
                    if score/min(lengthDownloaded, lengthOperon) > 0.90:
                        found = True

                if found:
                    # get necessary info - id, cog, category
                    id_operon = operon_proteins[last_index+j].id
                    id_downloaded = proteins[i].id

                    # if cog was assigned
                    try:
                        cog = cog_data.loc[cog_data["ID"] == id_operon, "COG"].values[0]
                        CATS = cog_data.loc[cog_data["ID"] == id_operon, "category"].values[0]
                        CAT = search('\w+', CATS).group(0)
                    # else assign "-"
                    except:
                        cog = "-"
                        CAT = "-"

                    # save this info into the table
                    new_row = pd.DataFrame(
                        {"protein_id": [id_downloaded], "source": ["operon_mapper"], "cog": [cog], "cat": [CAT]})
                    operon_table = pd.concat([operon_table, new_row], ignore_index=True)
                    break

        # if similar protein was found, save his position to continue were we left off
        if found:
            last_index = last_index + j

    f = open(output_dir + '/om_' + organism_name + '.txt', 'w')
    f.write('# created with COGtools 1.0.0\n# AC number: unknown\n# Processed data from Operon-mapper\n')
    f.close()
    return operon_table.to_csv(output_dir + '/om_' + organism_name + '.txt', sep='\t', index=False, mode = 'a')


def batch_processor_draft(organism_name, batch_file, output_dir=os.getcwd()):
    """
    Processes the outputs file (hitdata.txt) from Batch CD-Search tool into more structured COGtools-data.
    The output of this function is file in txt format that contains a suitable header with assigned COGs and their
    categories.
    :type organism_name: str
    :param batch_file: the path to Batch CD-Search outputs file hitdata.txt
    :param output_dir: the output file
    :return: processed file
    """
    # table for processed data
    batch_table = pd.DataFrame(columns=["protein_id", "source", "cog", "cat"])

    # data needed: annotated file by Batch CD-Search
    batch_data = (open(batch_file).read())
    batch_data = (batch_data[batch_data.index("Q#"):len(batch_data) - 1]).split('\n')
    query = ''

    # cog categories in COG database
    cogs_file = pkg_resources.resource_filename(__name__, 'COGtools-data/cogs.txt')
    cogs_data = (open(cogs_file, "r")).readlines()

    # iterate through queries
    for row in batch_data:
        new_query = search('Q#\d+', row).group(0)

        # if COG is not assigned or new_line is duplicate, continue to next iteration
        if not ('\tspecific\t' in row) or query == new_query:
            continue

        query = new_query
        seq_id = row[row.index(">") + 1:row.index(".")+2]

        try:
            COG = search("(COG\d+)", row).group(1)

        except AttributeError:
            COG = "-"

        # update in COG 2021
        if COG == "COG3512":
            COG = "COG1343"

        index = [i for i, s in enumerate(cogs_data) if COG[4:] in s][0]
        CAT = search('\t\w+', cogs_data[index]).group(0)[1:]

        # save this info into the table
        new_row = pd.DataFrame(
            {"protein_id": [seq_id], "source": ["batch cd-search"], "cog": [COG], "cat": [CAT]})
        batch_table = pd.concat([batch_table, new_row], ignore_index=True)

    f = open(output_dir + '/batch_' + organism_name + '.txt', 'w')
    f.write('# created with COGtools 1.0.0\n# AC number: unknown\n# Processed data from Batch CD-Search\n')
    f.close()
    return batch_table.to_csv(output_dir + '/batch_' + organism_name + '.txt', sep='\t', index=False, mode='a')
