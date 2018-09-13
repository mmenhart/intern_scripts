'''
This script is to filter out duplicate lines within the vcf file prior to submission for imputation.
It keeps in SNPs that are of the same position but measuring different minor alleles.
'''


from sys import argv
import os.path as op

vcf_file_dir = argv[1]   # ~/science/Cairo/checkVCF-20140116/from_scratch_vcf/switched2/
dup_file_dir = argv[2] # ~/science/Cairo/checkVCF-20140116/from_scratch_vcf/switched2/check_files/

CHR = argv[3]
chr = 'chr%s' % CHR

# This specifies the files in the directories
vcf_file = op.join(vcf_file_dir, '%s_switched4.vcf' % chr)
dup_file = op.join(dup_file_dir, '%s.check.dup' % chr)
out_file = 'Cairo/checkVCF-20140116/from_scratch_vcf1/switched/duplicates_filtered/%s_duplicate_check.vcf' % chr
uhoh_file = 'Cairo/checkVCF-20140116/from_scratch_vcf1/switched/duplicates/%s_inspect_duplicates.vcf' % chr
nomatch_file = 'Cairo/checkVCF-20140116/from_scratch_vcf1/switched/discordant_allele_calls/%s_mismatching_allele_calls.txt' % chr

uhohs = []
throwaway = []
garbage = []
nomatch = []
true_dup = []

# Open the files
with open(dup_file, 'r') as f:
    duplicates = [l.split(':')[1].strip('\n') for l in f]
    duplicates = list(set(duplicates))


with open(out_file, 'w') as fo:
    vcf = open(vcf_file).readlines()

    for i in range(0, len(vcf)):
        # Write the header to the vcf file
        if vcf[i].startswith('#'):
            continue

        # Specify the position in the vcf file
        split_line = vcf[i].split()
        position = split_line[1]
        # Check if the allele is present in the dup file
        for allele in duplicates:
            if position == allele:
                # Put somewhere else:
                uhohs.append(vcf[i]) # append the whole line, containing all information


    temp = []
    # Now check that duplicates in the uhoh file have the same allele calls
    for x in range(len(uhohs)-1):
        row1 = uhohs[x].split()
        row2 = uhohs[x+1].split()


        # First check that row1 and 2 have the same position and minor allele call
        if row1[1] == row2[1]: # If the positions are the same

            if len(temp)==0:
                temp.append(row1)
                temp.append(row2)
            else:
                temp.append(row2)


        else:

            for y in range(len(temp)):
                # You want to identify rows with common minor alleles, regardless of whether they're next to each other
                current = temp[y]

                for p in range(len(temp)-1):

                    if current[4] == temp[p][4] and current[2] != temp[p][2]: # If the minor alleles are the same and the refIDs are not
                        # Check the calls
                        for j in range(9, 261): # 9-260 are the genotype calls
                            # If it gets to the end of the alleles and hasn't been appended to throwaway, then all the calls are the same.

                            if j==260:
                                # Making it to the end means all the calls were the same and it's safe to throw one row out.
                                throwaway.append(current[2])
                                garbage.append(current[1])
                                temp.remove(temp[p])
                                temp.append(current) # to keep the temp length the same


                            elif current[j] != temp[p][j]: # append the whole line with all the information, get rid of both
                                throwaway.append(current[2])
                                throwaway.append(temp[p][2])
                                garbage.append(temp[p][1])
                                garbage.append(current[1])
                                nomatch.append(current[2])
                                nomatch.append(temp[p][2])

                                break
            temp=[]

    # Now go through the vcf file and write any line that isn't identified in throwawayself.
    for i in range(0, len(vcf)):
        # Write the header to the vcf file
        if vcf[i].startswith('#'):
            fo.write(vcf[i])


        else:
            ids = vcf[i].split()

            for p in range(len(throwaway)+1):
                # If you make it to the end of the throwaway list then there is no match and it should be written in the file
                if p == len(throwaway):
                    fo.write(vcf[i])
                # If it's in the throwaway list it should not be written. Skip that row/SNP.
                elif ids[2] == throwaway[p]:
                    break

with open(uhoh_file, 'w') as fu:
    for i in range(0, len(vcf)):
        # Write the header to the vcf file
        if vcf[i].startswith('#'):
            fu.write(vcf[i])

    for i in range(len(uhohs)):
        fu.write(uhohs[i])


with open(nomatch_file, 'w') as fn:
    for i in range(len(vcf)):
        if vcf[i].startswith('##'):
            continue
        elif vcf[i].startswith('#C'):
            fn.write(vcf[i])

        split = vcf[i].split()
        ref_id = split[2]

        for j in range(len(nomatch)):
            if ref_id == nomatch[j]:
                fn.write(vcf[i])


print "SNPs removed because of discordant allele calls are printed in %s_mismatching_allele_calls.vcf" % chr
print "Duplicate positions are printed in %s_inspect_duplicates.vcf" % chr
print "The final vcf file is present as %s_duplicate_check.vcf" % chr
print "Please run a final check on duplicate_check file using VCFtools."
