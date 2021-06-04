#### This pipeline can help trimming all genes alignment and remove the putative paralogs.
#### Any questions, please contact Wenbin Zhou. wzhou10@ncsu.edu

import os
from Bio import SeqIO
import numpy as np
import argparse
from Bio import AlignIO
from collections import Counter
import sys
from glob import glob
from itertools import tee
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC, Gapped

####

iupac = ['m', 'r', 'w', 's', 'y', 'k', 'M', 'R', 'W', 'S', 'Y', 'K', "v", "h", "d", "b", "V", "H", "D", "B"]
### define iupac code for Ploidy over 2
missing_data = ['-', '?', "n", "N"]
##### read each fasta sequences



### return species list as an input for step one.
def input_namelist(namelist_file):
    try:
        text1 = []
        f = open(namelist_file, 'rU')
        line = f.readline().strip()
        while line:
            # print(line)
            text1.append(line)
            line = f.readline().strip()
        f.close()
        print("Namelist file is ready!")
        # print (text1)
        return text1

    except:
        print("ERROR: Wrong path of namelist. Please check and reload it again!")



### s1 changed supercontig or exon files sorted by species to sorted by genes.
def sort_to_genes(seq_directory, output_directory, sp_namelist, seq_type = "supercontig"):

    ### use the output from output_directory()
    output_for_sorted_gene_s1 = output_directory + seq_type + "/s1_Gene/"
    if os.path.isdir(output_for_sorted_gene_s1) == False:
        os.makedirs(output_for_sorted_gene_s1)
        # print(output_for_sorted_gene_s1)
    ### use the sp_namelist from input_namelist()
    # print(sp_namelist)

    try:

        ### Please name you files from GATK and picard with .generated.fasta ending.
        for file in sp_namelist:
            # print(file)
            for file2 in os.listdir(seq_directory):
                # print(file2)
                if file in file2:
                    input_name = seq_directory + str(file2)
                    print("species " + file +" is ready")

                    for test1 in SeqIO.parse(input_name, 'fasta'):
                        print(test1)
                        sequence_name = test1.description.split(" ")[1].split("-")[0]
                        file_name = "Gene" + test1.description.split(" ")[1].split("-")[1].split(":")[0]

                        test1.id = sequence_name
                        test1.name = sequence_name
                        test1.description = sequence_name

                        out_file = output_for_sorted_gene_s1 + file_name + ".fasta"

                        #print out_file
                        f1 = open(out_file, 'a')
                        SeqIO.write(test1, f1, "fasta")
                        f1.close()
    except:
        print("Please check the input files.")
            # else:
        #     print("please check the name in the namelist are uniquely matching to a part of degenerated fasta file.")
    return output_for_sorted_gene_s1





### s2 seq_directory is from sort_to_genes()
def rm_Hs_sequence(seq_directory, namelist_path, Hsite_value =0.05):

    # ### check each file
    output_directory = seq_directory.replace("s1_Gene/", "")

    ### mkdir the output directory for s2
    out_file0 = output_directory + "s2_filtered_Gene/Gene/"

    ### status file for s2
    H_site_summary = output_directory + "s2_filtered_Gene/Hsite_status.txt"

    if os.path.isdir(out_file0) == False:
        os.makedirs(out_file0)


    # ### loop all fasta files
    sum = ""
    for file in os.listdir(seq_directory):
         ### This is for mac user.
        if file != ".DS_Store":
            # print(file)

            ### subset is the input fasta file from s1.
            subset = seq_directory + file
            # print(subset)

            # ### read fasta file
            for sequence in SeqIO.parse(subset, 'fasta'):

                if sequence.id in namelist_path:
                    # print(sequence)

                    sequence_total_len = len(sequence)
                    # print(sequence_total_len)

                    ### change all the characters in to an array.
                    seq_array = np.array([list(rec) for rec in sequence.seq])
                    # print(seq_array)

                    ### if find the bases are from iupac code, put it in the list.
                    iupac_in_alignment = [ele for ele in seq_array if (ele in iupac)]

                    ### count the length of the list, i.e., the number of polymorphic sites.
                    H_site_number = len(iupac_in_alignment)
                    # print (H_site_number)

                    ### put the missing data in the list and count them.
                    missingdata_number = len([ele for ele in seq_array if (ele in missing_data)])
                    # print (missingdata_number)

                    ### calculate the heterozygous site percentage.
                    H_site_percent = float(H_site_number) / float(sequence_total_len - missingdata_number)
                    # print (H_site_percent)

                    ### note down a table for the heterozyous status.
                    sum = sum + file.split(".")[0] + "\t" + sequence.id + "\t" + str(H_site_number) + "\t" + str(
                        H_site_percent) + "\n"
                    # print(sum)

                    ### if the Hsite percentage is less than a value, keep the sequence in the loci.
                    ### write the result to a new file.
                    if H_site_percent <= Hsite_value:
                        out_file = out_file0 + file
                        # print(out_file)

                        f1 = open(out_file, 'a')
                        SeqIO.write(sequence, f1, "fasta")
                        f1.close()

    ### generate final table with an header.
    H_site_info = "Gene" + "\t" + "Accession" + "\t" + "Ha_max_over50_number" + "\t" + "Hs_max" + "\n" + sum
    # print(H_site_info)

    f1 = open(H_site_summary, 'a')
    f1.write(H_site_info)
    f1.close()




### s3 add reference to each gene file based on s2
def add_ref_to_each_gene(seq_directory, reference_name):
    ### /gpfs_common/share02/qyxiang/wzhou10/Cornus//result/supercontig/s1_Gene/

    ### the input files are from s2
    genes_result_s2 = seq_directory.replace("s1_Gene/","s2_filtered_Gene/Gene/*")
    # print(genes_result_s2)

    ### mkdir directory for s3.
    genes_result_s3 = seq_directory.replace("s1_Gene/", "s3_add_ref_Gene/")
    if os.path.isdir(genes_result_s3) == False:
        os.makedirs(genes_result_s3)
        print(genes_result_s3)
    # print genes_result_s3

    ### if the s2 directory is not empty, run mafft.
    try:
        ### copy the genes from s2 to s3 folder for adding the ref to the current genes
        cmd = "cp -R " + genes_result_s2 + " " + genes_result_s3
        print(cmd)
        os.system(cmd)

        ### reference name should be like "SPECIESNAME-GENENAME". Please Capitalize them.
        ### add the ref genes to the matching genes, with same gene names, which is from Hybpiper.
        for test1 in SeqIO.parse(reference_name, 'fasta'):

            sequence_name = test1.description.split("-")[0]
            file_name = "Gene" + test1.description.split("-")[-1]

            test1.id = sequence_name
            test1.name = sequence_name
            test1.description = sequence_name

            out_file = genes_result_s3 + file_name + ".fasta"

            if os.path.isfile(out_file) == True:
                #        print test1
                #        print out_file
                f1 = open(out_file, 'a')
                SeqIO.write(test1, f1, "fasta")
                f1.close()

    except:
        print("s2 is not finished yet. Please check the previous step.")



### s4 align with reference to improve the quality of alignment
### make sure you install the mafft. Recommend to use mminiconda for installation.
def mafft_alignment(seq_directory):

    ### input files are from s3
    genes_result_s3 = seq_directory.replace("s1_Gene/", "s3_add_ref_Gene/")

    ### mkdir output file directory
    genes_result_s4 = seq_directory.replace("s1_Gene/", "s4_alignments/")
    if os.path.isdir(genes_result_s4) == False:
        os.makedirs(genes_result_s4)
        print(genes_result_s4)

    ### if the s3 directory is not empty, run mafft.
    try:
        for file in os.listdir(genes_result_s3):
            fname = genes_result_s3 + file
            output_name = genes_result_s4 + file
            if file != ".DS_Store":
                cmd = "mafft --adjustdirection --maxiterate 1000 --globalpair " + fname + " > " + output_name
                print(cmd)
                os.system(cmd)

    except:
        print("s3 is not finshed yet. Please check the previous step.")




### s5 remove reference sequence
def rm_reference_from_align(seq_directory, reference_name):

    ### input files are from s5
    genes_result_s4 = seq_directory.replace("s1_Gene/", "s4_alignments/")

    ### mkdir output directory for s6
    genes_result_s5 = seq_directory.replace("s1_Gene/", "s5_rm_ref_seq/")
    if os.path.isdir(genes_result_s5) == False:
        os.makedirs(genes_result_s5)
        # print(genes_result_s6)

    ### keep the sequences from namelist.
    try:
        for file1 in os.listdir(genes_result_s4):
            if file1 != ".DS_Store":
                fasta_name1 = genes_result_s4 + file1
                for gene_name in reference_name:
                    print(gene_name)
                    for test1 in SeqIO.parse(fasta_name1, 'fasta'):
                        if gene_name == test1.name:
                            print(test1)
                            output_file = genes_result_s5 + file1
                            f1 = open(output_file, 'a')
                            SeqIO.write(test1, f1, "fasta")
                            f1.close()
        ### replace all n to -, otherwise trimal can not run successfully.

        cmd_replace_N = "FILES=" + genes_result_s5 + " && for f in $FILES*; do sed -i '' -e 's/n/-/g' $f; done"
        os.system(cmd_replace_N)
    except:
        print("s4 is not finshed yet. Please check the previous step.")



### s6 trimal the alignments again.
def trimal_align_2(seq_directory, gap_fraction = 0.51):

    ### input files are from s6
    genes_result_s5 = seq_directory.replace("s1_Gene/", "s5_rm_ref_seq/")

    ### mkdir output directory for s7
    genes_result_s6 = seq_directory.replace("s1_Gene/", "s6_trimal/")
    if os.path.isdir(genes_result_s6) == False:
        os.makedirs(genes_result_s6)
        print(genes_result_s6)

    ### run trimal again.
    try:
        for file in os.listdir(genes_result_s5):
            fname = genes_result_s5 + file
            output_name = genes_result_s6 + file
            if file != genes_result_s5 + ".DS_Store":
                cmd = "trimal -in " + fname + " -out " + output_name + " -gt " + str(gap_fraction)
                print(cmd)
                os.system(cmd)

    except:
        print("s5 is not finshed yet. Please check the previous step.")



### s7 many functions for a better trimmed result

### sliding window function.
def window(iterable, size):
    iters = tee(iterable, size)
    for i in range(1, size):
        for each in iters[i:]:
            next(each, None)
    return zip(*iters)


### read outgroups species and return the list of outgroups.
def input_outgroup(outgroup_file):

    try:
        text1 = []
        f = open(outgroup_file)
        line = f.readline().strip()
        while line:
            # print(line)
            text1.append(line)
            line = f.readline().strip()
        f.close()
        # self.text1 = text1
        # return text1
        if "\r" in text1:
            print("Please reinput all species by manually")
        else:
            print("Namelist file is ready!")
        return text1
    except:
        print("ERROR: Wrong path of namelist. Please check and reload it again!")


### find the likely wrong sites with high polymorphic sites percentage in a sliding window.
### then use trimal to remove them.
def rm_wrong_polymorphism_sites(seq_directory, outgroup_path, window_size = 20, Max_p_sites = 4):
    ### define iupac
    iupac_bases = ['m', 'r', 'w', 's', 'y', 'k', 'M', 'R', 'W', 'S', 'Y', 'K', "v", "h", "d", "b", "V", "H", "D", "B"]

    ### input files are from s6
    genes_result_s6 = seq_directory.replace("s1_Gene/", "s6_trimal/")

    ### mkdir output directory for s7
    genes_result_s7 = seq_directory.replace("s1_Gene/", "s7_well_trimal/")

    ### return outgroup list
    outgroups = input_outgroup(outgroup_path)

    output_directory = genes_result_s7 + "/s1_rm_polymorphism_sites/"
    if os.path.isdir(output_directory) == False:
        os.makedirs(output_directory)

    ### iterate each gene
    for file in os.listdir(genes_result_s6):
        if file != ".DS_Store":
            output_directory_file = output_directory + file
            fasta_name = genes_result_s6 + file

            sequences = glob(fasta_name)
            ### read each alignment sequences
            for sequence in sequences:
                print("sequence: " +sequence)

                alignment = AlignIO.read(sequence, 'fasta')
                # print(alignment)

                ### generate a new alignment sequences without outgroups.
                align = MultipleSeqAlignment([], Gapped(IUPAC.unambiguous_dna, "-"))

                for record in alignment:
                    if record.id not in outgroups:
                        # print(record.id)
                        # print(record.seq)
                        align.add_sequence(str(record.id), str(record.seq))


                print(align)
                # print(align.get_alignment_length())


                total_wrong_poly_sites = []
                ### change alignment to an array.
                align_array = np.array([list(rec) for rec in align])
                ### , np.character
                # print(align_array)

                ### calculate the whole length of the alignment
                total_length = align.get_alignment_length()



                ### using 20bp-long sliding windows.
                for each in window(range(total_length), window_size):
                    # print(list(each))
                    poly_site_no_iupac = 0
                    poly_site_number = 0

                    column_position = []

                    ### for each block calculate the polymorphism sites number.
                    for column in each:
                        ### calculate each site (each column).
                        counter = Counter(align_array[:, column])

                        ### sorted by frequency
                        sorted_bases = counter.most_common()

                        # print(counter)
                        # print(sorted_bases)
                        # print(len(counter))

                        ### count the sites with different situations.
                        gap_yes = 0

                        if len(counter) ==1:
                            poly_site_number = poly_site_number + 0
                            poly_site_no_iupac = poly_site_no_iupac + 0


                        elif len(counter) == 2:
                            for i in sorted_bases:
                                if i[0] == "-":
                                    gap_yes = 1
                                else:
                                    gap_yes = 0
                            # print("gap is 1 or 0:" + str(gap_yes))

                            if gap_yes == 1:
                                # print counter
                                poly_site_number = poly_site_number + 0
                                poly_site_no_iupac = poly_site_no_iupac + 0

                            else:
                                iupac_in_alignment = [ele for ele in sorted_bases if (ele[0] in iupac_bases)]
                                # print(iupac_in_alignment)
                                if len(iupac_in_alignment) == 1:
                                    poly_site_number = poly_site_number + 1
                                    poly_site_no_iupac = poly_site_no_iupac + 0

                                if len(iupac_in_alignment) == 0:
                                    poly_site_number = poly_site_number + 1
                                    poly_site_no_iupac = poly_site_no_iupac + 1
                                    # print(column)
                                    column_position.append(column)

                        elif len(counter) == 3:
                            for i in sorted_bases:
                                if i[0] == "-":
                                    gap_yes = 1
                                else:
                                    gap_yes = 0
                            # print("gap is 1 or 0:" + str(gap_yes))

                            if gap_yes == 1:
                                # print counter

                                iupac_in_alignment = [ele for ele in sorted_bases if (ele[0] in iupac_bases)]
                                # print(iupac_in_alignment)
                                if len(iupac_in_alignment) == 1:
                                    # poly_site_no_iupac = poly_site_no_iupac + 1
                                    poly_site_number = poly_site_number + 1
                                    poly_site_no_iupac = poly_site_no_iupac + 0

                                else:
                                    poly_site_number = poly_site_number + 1
                                    poly_site_no_iupac = poly_site_no_iupac + 1
                                    # print(column)
                                    column_position.append(column)

                            else:
                                poly_site_number = poly_site_number + 1
                                poly_site_no_iupac = poly_site_no_iupac + 1
                                # print(column)
                                column_position.append(column)


                        else:
                            poly_site_number = poly_site_number + 1
                            poly_site_no_iupac = poly_site_no_iupac + 1
                            # print(column)
                            column_position.append(column)


                    # print("column_position: " + str(column_position))
                    # print(len(column_position))

                    ### if there are more than 4 polymorphic sites in 20 base pairs, select those sites positions.
                    if len(column_position) > float(Max_p_sites):
                        print(column_position)
                        total_wrong_poly_sites = total_wrong_poly_sites + column_position

                #print(total_wrong_poly_sites)

                ### generate the unique positions

                total_wrong_poly_sites = total_wrong_poly_sites + list(range(10))
                total_wrong_poly_sites = total_wrong_poly_sites + list(range(total_length-10, total_length))
                ### extract the polymorphic sites from alignment data, might be useful for delete the first 2 species.
                unique_wrong_sites = list(np.unique(total_wrong_poly_sites))
                print(len(unique_wrong_sites))
                # sum2 = alignment[:, total_length:total_length + 1]
                # for i in unique_wrong_sites:
                #     sum2 = sum2 + alignment[:, i:i+1]
                # print(sum2)
                # SeqIO.write(sum2, "/Users/zhouwenbin/Downloads/result/M40_total.phy", "phylip")


                ### operating: if any window has more than 3 polymorphic sites, use trimal to remove those sites.
                ### otherwise, copy the gene to the new folder.
                if len(unique_wrong_sites) > 0:

                    print(str(unique_wrong_sites).replace(" ", "").replace("[", "\{ ").replace("]", " \}"))

                    cmd_selected_col = str(unique_wrong_sites).replace(" ", "").replace("[", "\{ ").replace("]", " \}")

                    cmd = "trimal -in " + fasta_name + " -out " + output_directory_file + " -selectcols " + cmd_selected_col

                    print(cmd)
                    os.system(cmd)

                else:
                    cmd_2 = "cp " + fasta_name + " " + output_directory_file
                    print(cmd_2)
                    os.system(cmd_2)



def replace_outgroup_with_gap(seq_directory, outgroup_path, window_size = 20, Max_p_sites_o = 8):
    ### define iupac
    iupac_bases = ['m', 'r', 'w', 's', 'y', 'k', 'M', 'R', 'W', 'S', 'Y', 'K', "v", "h", "d", "b", "V", "H",
                   "D", "B"]

    ### input directory from s7
    genes_result_s7 = seq_directory.replace("s1_Gene/", "s7_well_trimal/")

    ### return outgroup list
    outgroups = input_outgroup(outgroup_path)

    output_directory_1 = genes_result_s7 + "/s1_rm_polymorphism_sites/"
    output_directory_2 = output_directory_1.replace("/s1_rm_polymorphism_sites/","/s2_rm_polymorphism_in_outgroups/")

    if os.path.isdir(output_directory_2) == False:
        os.makedirs(output_directory_2)

    ### iterate each gene
    for file in os.listdir(output_directory_1):
        if file != ".DS_Store":
            output_directory_file = output_directory_2 + file
            fasta_name = output_directory_1 + file

            sequences = glob(fasta_name)
            ### read each alignment sequences
            for sequence in sequences:
                print("sequence: " + sequence)

                alignment = AlignIO.read(sequence, 'fasta')

                ### calculate the polymorphism in outgroup
                ### change alignment to an array.
                total_wrong_poly_sites_outgroup = []

                align_array_outgroup = np.array([list(rec) for rec in alignment])
                ### , np.character
                # print(align_array)

                ### calculate the whole length of the alignment
                total_length = alignment.get_alignment_length()
                # alignment = AlignIO.read(sequence, 'fasta')
                for each in window(range(total_length), window_size):
                    # print(list(each))
                    poly_site_no_iupac = 0
                    poly_site_number = 0

                    column_position_outgroup = []

                    ### for each block calculate the polymorphism sites number.
                    for column in each:
                        ### calculate each site (each column).
                        counter = Counter(align_array_outgroup[:, column])

                        ### sorted by frequency
                        sorted_bases = counter.most_common()

                        # print(counter)
                        # print(sorted_bases)
                        # print(len(counter))

                        ### count the sites with different situations.
                        gap_yes = 0

                        if len(counter) ==1:
                            poly_site_number = poly_site_number + 0
                            poly_site_no_iupac = poly_site_no_iupac + 0


                        elif len(counter) == 2:
                            for i in sorted_bases:
                                if i[0] == "-":
                                    gap_yes = 1
                                else:
                                    gap_yes = 0
                            # print("gap is 1 or 0:" + str(gap_yes))

                            if gap_yes == 1:
                                # print counter
                                poly_site_number = poly_site_number + 0
                                poly_site_no_iupac = poly_site_no_iupac + 0

                            else:
                                iupac_in_alignment = [ele for ele in sorted_bases if (ele[0] in iupac_bases)]
                                # print(iupac_in_alignment)
                                if len(iupac_in_alignment) == 1:
                                    poly_site_number = poly_site_number + 1
                                    poly_site_no_iupac = poly_site_no_iupac + 0

                                if len(iupac_in_alignment) == 0:
                                    poly_site_number = poly_site_number + 1
                                    poly_site_no_iupac = poly_site_no_iupac + 1
                                    # print(column)
                                    column_position_outgroup.append(column)

                        elif len(counter) == 3:
                            for i in sorted_bases:
                                if i[0] == "-":
                                    gap_yes = 1
                                else:
                                    gap_yes = 0
                            # print("gap is 1 or 0:" + str(gap_yes))

                            if gap_yes == 1:
                                # print counter

                                iupac_in_alignment = [ele for ele in sorted_bases if (ele[0] in iupac_bases)]
                                # print(iupac_in_alignment)
                                if len(iupac_in_alignment) == 1:
                                    # poly_site_no_iupac = poly_site_no_iupac + 1
                                    poly_site_number = poly_site_number + 1
                                    poly_site_no_iupac = poly_site_no_iupac + 0

                                else:
                                    poly_site_number = poly_site_number + 1
                                    poly_site_no_iupac = poly_site_no_iupac + 1
                                    # print(column)
                                    column_position_outgroup.append(column)

                            else:
                                poly_site_number = poly_site_number + 1
                                poly_site_no_iupac = poly_site_no_iupac + 1
                                # print(column)
                                column_position_outgroup.append(column)


                        else:
                            poly_site_number = poly_site_number + 1
                            poly_site_no_iupac = poly_site_no_iupac + 1
                            # print(column)
                            column_position_outgroup.append(column)


                    # print("column_position: " + str(column_position))
                    # print(len(column_position))

                    ### if there are more than 8 polymorphic sites in 20 base pairs, select those sites positions.
                    if len(column_position_outgroup) > float(Max_p_sites_o):
                        print(column_position_outgroup)
                        total_wrong_poly_sites_outgroup = total_wrong_poly_sites_outgroup + column_position_outgroup


                unique_wrong_sites_ougroup = list(np.unique(total_wrong_poly_sites_outgroup))
                print(unique_wrong_sites_ougroup)
                print("outgroup")


                align_2 = MultipleSeqAlignment([], Gapped(IUPAC.unambiguous_dna, "-"))
                for record in alignment:
                    new_seq = ""

                    if record.id in outgroups:
                        print(record.seq)
                        for i in range(total_length):
                            if i in unique_wrong_sites_ougroup:
                                new_seq = new_seq + "-"
                            else:
                                new_seq = new_seq + str(record.seq[i])

                        align_2.add_sequence(str(record.id), str(new_seq))

                    else:
                        align_2.add_sequence(str(record.id), str(record.seq))

                print(align_2)

                AlignIO.write(align_2, output_directory_file, "fasta")





### s8 romove putative paralogs.
def paralog_test(seq_directory, Hs_threshold = 0.5, nh = 1):

    ### input files are from s7
    genes_result_s7 = seq_directory.replace("s1_Gene/", "s7_well_trimal/s2_rm_polymorphism_in_outgroups/")

    ### mkdir output directory for s8
    genes_result_s8 = seq_directory.replace("s1_Gene/", "s8_rm_paralogs/")

    ### mkidir for orthologs, paralogs and short sequences.
    out_file0 = genes_result_s8 + "Final_kept_genes/"
    out_file1 = genes_result_s8 + "paralogs_warning/"
    rest_files = genes_result_s8 + "discard_file/"
    paralog_summary = genes_result_s8 + "paralogs_status.txt"
    # print(paralog_summary)

    if os.path.isdir(out_file0) == False:
        os.makedirs(out_file0)
    if os.path.isdir(out_file1) == False:
        os.makedirs(out_file1)
    if os.path.isdir(rest_files) == False:
        os.makedirs(rest_files)

    try:
        info_paralogs = ""

        ### loop all alignment files from s7
        for file in os.listdir(genes_result_s7):
            ### This is for mac user.
            if file != ".DS_Store":
                #print (file)

                ### subset is the input fasta file
                subset = genes_result_s7 + file

                Hs_max = 0
                H_number = 0

                ### read fasta file
                alignment = AlignIO.read(subset, 'fasta')
                # print(alignment)

                ### change all the characters in to an array.
                align_array = np.array([list(rec) for rec in alignment])
                # print(align_array)

                ### calculate number of species.
                number_of_species = len(align_array[:, 0])
                # print("number of species: " + str(number_of_species))

                ### get total length of alignment.
                total_length = alignment.get_alignment_length()
                # print total_length

                ### loop each site to find heterozygosity.
                for column in range(alignment.get_alignment_length()):
                    #print(align_array[:,column])

                    ### 1. calculate the heterzygous site number
                    iupac_in_alignment = [ele for ele in align_array[:, column] if (ele in iupac)]
                    #print(iupac_in_alignment)

                    if len(iupac_in_alignment) > 0:
                        H_number = H_number + 1
                        # print(H_number)

                    ### 2. calculate the heterzygous site number with over than 50% individuals in an alignment
                    counter = Counter(align_array[:, column])
                    # print(counter)
                    sorted_bases = counter.most_common()

                    ### calculate most common ambiguity percentage at each site.
                    for i in sorted_bases:
                        # print(i)
                        if i[0] in iupac:
                            # print(i[1])
                            # calculate the shared heterozygous percentage.
                            base_percentage = float(i[1])/float(number_of_species)

                            ### if the shared heterozygous percentage is over threshold, the loci whill be removed.
                            ### default setting is 50%.
                            if base_percentage >= Hs_threshold:

                                ### count the total number of shared heterozygous sites.
                                Hs_max = Hs_max + 1
                                # print(i[0] + " " + file)
                            ### it can be changed to a number of shared heterozygous sites.
                            break

                ### calculate the percentage of Maximum heterozygous site.
                H_percentage = float(H_number)/float(total_length)
                H_percentage_3 = round(H_percentage,3)
                #print(Hs_max)

                ### output for status file.
                sum = file.split(".")[0] + "\t" + str(H_percentage_3) + "\t" + str(Hs_max) + "\n"
                # print(sum)

                ### categorize genes to new folders.
                ### aligned length over 90 bases, and at least 4 sequences.
                ### kept genes.
                if total_length > 90 and len(align_array[:, 1]) > 3:
                    ### orthologs
                    if Hs_max <= int(nh)-1:
                        cmd = "cp " + subset + " " + out_file0
                        print(cmd)
                        os.system(cmd)
                    ### paralogs with 50% shared heterozygous sites.
                    if Hs_max > int(nh)-1:
                        cmd = "cp " + subset + " " + out_file1
                        print(cmd)
                        os.system(cmd)

                ### discarded genes.
                else:
                    cmd = "cp " + subset + " " + rest_files
                    print(cmd)
                    os.system(cmd)
                    # print(sum)

                info_paralogs = info_paralogs + sum
                # print(info_paralogs)
        final_info_paralogs = "Gene" + "\t" + "Ha_max_over50_number" + "\t" + "Hs_max" + "\n" + info_paralogs
        # print(final_info_paralogs)

        f1 = open(paralog_summary, 'a')
        f1.write(final_info_paralogs)
        f1.close()

    except:
        print("s7 is not finished yet. Please check the previous step.")






def main():
    parser = argparse.ArgumentParser(description="To find putative paralogs in 353 enrichment data. It requires three input components and one output component. The '-he' and '-hs' are not necessary. Dependencies: mafft and trimal installed by conda. The '-ifa', '-ina', '-iref', '-io', '-o' are required arguments.",
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-ifa", "--inputfastadir", dest='input_dir', type=str,
                        help='the directory where degenerated fasta files deposits, which are generated by Step2')

    parser.add_argument("-ina", "--inputspna", dest='input_sp_name', type=str,
                        help='the full path of namelist.txt with every accession in one line')

    parser.add_argument("-iref", "--inputref", dest='input_ref_seq', type=str,
                        help='the full path of reference seq, which is as same as the one in Hybpiper. Fasta format. The reference name should be like "SPECIESNAME-GENENAME". Please Capitalize them.')
    #### modify
    parser.add_argument("-io", "--inputoutgroup", dest='input_outgroup', type=str,
                        help='the full path of outgroups.txt. Every line has one accession')
    # parser.add_argument("-s", "--step", dest='input_steps', choices=range(1,9), type=int,
    #                     help='input steps you want to run with. Range from step 1 to 8')

    parser.add_argument("-o", "--output", dest='output_dir', type=str,
                        help='output path. You can run supercontig and exon with same path, but you should change your path when you are running another group of taxa.')

    parser.add_argument("-t", "--seqtype", dest='seq_type', type=str, default="supercontig",
                        help='It can be supercontig or exon data type. Default is supercontig')

    parser.add_argument("-he", "--Hsites", dest='Hsite_value', type=float, default=0.05,
                        help='Maximum percentage of heterozygous site in a sequence. The default value is 0.05 (means 5%%).')
    
    parser.add_argument("-gt", "--gapthreshold", dest='gap_threshold', type=float, default=0.51,
                        help='1 - (fraction of sequences with a gap allowed). See details in trimmAl -gt')

    parser.add_argument("-hs", "--Hs", dest= 'Hs_max_value',type=float, default=0.5,
                        help='Maximum shared heterozygous percentage at a site. The default value is 0.5 (means 50%%).')

    parser.add_argument("-nh", "--numberH", dest= 'Hs_site_number',type=int, default=1,
                        help='The number of sites with maximum shared heterozygous percentage. The default value is 1, which means if one or more sites were detected, the locus will be considered as paralogs.')

    parser.add_argument("-w", "--windowsize", dest='window_size', type=int, default=20,
                        help='setting the sliding window size. The default is 20 bp.')

    parser.add_argument("-mi", "--maxpsite", dest='max_p_site', type=int, default=4,
                        help='The allowed maximum number of polymorphism site in each sliding window without outgroup. Default is 4.')
    
    parser.add_argument("-mo", "--maxpsiteo", dest='max_p_site_o', type=int, default=8,
                        help='The allowed maximum number of polymorphism site of outgroup in each sliding window. Default is 8.')

    args = parser.parse_args()
    # print(args)
    f = open("PPD.log", "w")
    sys.stdout = f



    if args.input_dir and args.input_sp_name and args.input_ref_seq and args.output_dir and args.input_outgroup:
        # print("inputs are ready!")


        input_directory = os.path.realpath(args.input_dir) + "/"
        print("input directory is ready: " + input_directory)

        #if args.output_dir[-1] != "/":
        output_directory = os.path.realpath(args.output_dir) + "/result/"
        print("output directory is ready: " + output_directory)
        # sys.stdout.write("output directory is ready: " + output_directory)

        sp_name = input_namelist(os.path.realpath(args.input_sp_name))
        print("species names are: " + str(sp_name))

        seq_result_path = sort_to_genes(input_directory, output_directory, sp_name, seq_type = args.seq_type)
        print("s1 results are here: " + seq_result_path)


        ### all the following steps are based on s1.
        if len(os.listdir(seq_result_path)) <= 1:
            print("s1 is not finished!!! Please check your inputs again.")

        else:
            print("s1 finished and ready for s2")

            ### s2 filter the individual gene sequence with high heterozygous sites, threshold is 0.05.
            rm_Hs_sequence(seq_result_path, sp_name, Hsite_value = args.Hsite_value)
            print("s2 finished and ready for s3")

            ### s3 input path and reference name
            add_ref_to_each_gene(seq_result_path, os.path.realpath(args.input_ref_seq))
            print("s3 finished and ready for s4")

            ### s4
            mafft_alignment(seq_result_path)
            print("s4 finished and ready for s5")

            ### s5 is based on s1 and s2 results.
            rm_reference_from_align(seq_result_path, sp_name)
            print("s6 finished and ready for s7")

            ### s6 trimal gap
            trimal_align_2(seq_result_path, gap_fraction = args.gap_threshold)
            print("s7 finished and ready for s8")


            ### s7 well trimmed
            rm_wrong_polymorphism_sites(seq_result_path, os.path.realpath(args.input_outgroup), args.window_size,
                                        args.max_p_site)

            replace_outgroup_with_gap(seq_result_path, os.path.realpath(args.input_outgroup), args.window_size, args.max_p_site_o)

            ### s8
            paralog_test(seq_result_path, Hs_threshold = args.Hs_max_value, nh = args.Hs_site_number)
            print("s8 finished and ready for s9")


    else:
        print("Please check the parameter settings with -h argument!!!")

    f.close()

if __name__ == '__main__':
    main()
