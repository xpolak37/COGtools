import sys
import program_processor
import track_manager
import consensus
import os
import argparse


def cogtools(organism_name, input_dir, output_dir, cogs_only, cat_choice, manager, cog_colors, draft, gff_included):
    try:
        # complete genome
        em_file = None
        om_file = None
        batch_file = None
        if draft is False:
            files = [f for f in os.listdir(input_dir)]
            fasta_file = input_dir + "/" + organism_name + ".fasta" if organism_name + ".fasta" in files else None
            gff_file = input_dir + "/" + organism_name + ".gff3" if organism_name + ".gff3" in files else None

            # Program processor
            if organism_name + "_eggnog.gff" in files:
                program_processor.em_processor(organism_name, input_dir + "/" + organism_name + "_eggnog.gff",
                                               input_dir + "/" +
                                               organism_name + "_cds.txt", cogs_only, output_dir)
                em_file = output_dir + "/em_" + organism_name + ".gff"

            if organism_name + "_orf_operon.txt" in files:
                program_processor.om_processor(organism_name, input_dir + "/" + organism_name + "_orf_operon.txt",
                                               input_dir +
                                               "/" + organism_name + "_cogs_operon.txt", output_dir)
                om_file = output_dir + "/om_" + organism_name + ".gff"

            if organism_name + "_batch.txt" in files:
                program_processor.batch_processor(organism_name, input_dir + "/" + organism_name + "_batch.txt", output_dir)
                batch_file = output_dir + "/batch_" + organism_name + ".gff"

            # Consensus
            if gff_file is None:
                get_pseudo=False
                get_ncrna=False
            else:
                get_pseudo = True
                get_ncrna = True

            consensus.consensus(organism_name, em_file=em_file,
                                om_file=om_file,
                                batch_file=batch_file, fasta_file=fasta_file,
                                get_pseudo=get_pseudo, get_ncrna=get_ncrna,
                                gff_file=gff_file, output_dir=output_dir, cat_choice=cat_choice)

            # Track manager
            if manager:
                track_manager.get_track_template(cog_palette=cog_colors, output_dir=output_dir)
                track_manager.get_legend(output_dir=output_dir, cog_palette=cog_colors)

        # draft genome
        else:
            files = [f for f in os.listdir(input_dir)]
            proteins_file = input_dir + "/" + organism_name + "_proteins.fsa_aa" \
                if organism_name + "_proteins.fsa_aa" in files else None
            # Program processor
            if organism_name + "_eggnog.gff" in files:
                program_processor.em_processor_draft(organism_name, input_dir + "/" + organism_name + "_eggnog.gff",
                                                     output_dir=output_dir, cogs_only=cogs_only)
                em_file = output_dir + "/em_" + organism_name + ".txt"

            if organism_name + "_proteins_operon.txt" in files:
                program_processor.om_processor_draft(organism_name, input_dir + "/" + organism_name + "_proteins.fsa_aa",
                                                     input_dir + "/" + organism_name + "_proteins_operon.txt",
                                                     input_dir + "/" + organism_name + "_cogs_operon.txt",
                                                     output_dir=output_dir,gff_included=gff_included)
                om_file = output_dir + "/om_" + organism_name + ".txt"

            if organism_name + "_batch.txt" in files:
                program_processor.batch_processor_draft(organism_name, input_dir + "/" + organism_name + "_batch.txt",
                                                        output_dir=output_dir)
                batch_file = output_dir + "/batch_" + organism_name + ".txt"

            # Consensus
            consensus.consensus_draft(organism_name,
                                      proteins=proteins_file,
                                      em_file=em_file,om_file=om_file, batch_file=batch_file,
                                      cat_choice=cat_choice,
                                      output_dir=output_dir)

    except Exception as e:
        print(e)
        sys.exit(2)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='COGtools 1.0.0')
    parser.add_argument("-n", "--name", action="store", dest="organism_name", default=None)
    parser.add_argument("-i", "--input", action="store", dest="input_dir", default=os.getcwd())
    parser.add_argument("-o", "--output", action="store", dest="output_dir", default=os.getcwd())
    parser.add_argument("-c", action="store_true", dest="cogs_only")
    parser.add_argument("-ch", "--choice", action="store", dest="cat_choice", default=1, type=int)
    parser.add_argument("-t", "--track", action="store_true", dest="track_manager")
    parser.add_argument("-p", action="store_true", dest="cogs_palette")
    parser.add_argument("-d", "--draft", action="store_true", dest="draft")
    parser.add_argument("-g", "--gff", action="store_true", dest="gff_included")
    arguments = parser.parse_args()

    cogtools(arguments.organism_name, arguments.input_dir, arguments.output_dir,arguments.cogs_only,
             arguments.cat_choice, arguments.track_manager,arguments.cogs_palette, arguments.draft,
             arguments.gff_included)