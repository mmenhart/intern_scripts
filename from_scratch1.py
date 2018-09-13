'''
This script will take the original genomestudio file calls per sample as well
as the snp table provided by genomestudio and produce a vcf file.
'''

from sys import argv
import os.path as op


CHR = argv[1]
chr = 'chr%s' % CHR

snp_calls_dir = argv[2]   # ~/shared/meijer/CAIRO_SNP-data/illumina_snp_calls
snp_table = 'Cairo/SNP_Table.csv'

snp_file = op.join(snp_calls_dir,'%s_illumina_snpcalls_sorted.csv' % chr)
out_file = 'Cairo/checkVCF-20140116/from_scratch_vcf1/%s_calls.vcf' % chr
vcf_header = 'Cairo/vcf_header.txt'

# Create a snp dictionary
with open(snp_table,'r') as f:
	dict = {str(l.split(',')[0]) : str(l.split(',')[3]) for l in f}
	# Lookup by SNP NAME not position (since some are duplicated). Name:[A/G]

with open(out_file,'w') as fo:
    # write header of vcf file
	with open(vcf_header, 'r') as vh:
		for line in vh:
			fo.write(line)
	fo.write('\n')

	with open(snp_file,'r') as fi:
		for l in fi:
			ls1 = l.strip()
			ls = ls1.split(',')
			# Skip the header of the snp file
			if ls[0] == 'Name':
				continue

			chrom = str(ls[1])
			pos = str(ls[2])
			snp_name = str(ls[0].strip())

			# Change the BB/AA/AB notation into 0/0 0/1 1/1 calls
			snp_alleles = dict[snp_name].strip('[]').split('/')
			A_allele = snp_alleles[0]
			B_allele = snp_alleles[1]

			new_line = chrom+'\t'+pos+'\t'+snp_name+'\t'

			# Now write the proper vcf line calls
			new_line = new_line + A_allele + "\t" + B_allele + '\t'
			# Add necessary extras
			new_line = new_line+'.'+'\t'+'.'+'\t'+'PR'+'\t'+'GT'
			# Change the calls so AA = '0/0', BB = '1/1', etc. We will change the major/minor if needed after running checkVCF.py
			for i in range(2, 254):

				if ls[i] == 'AA':
					new_line = new_line + '\t' + '0/0'
				elif ls[i] == 'AB':
					new_line = new_line + '\t' + '0/1'
				elif ls[i] == 'BB':
					new_line = new_line + '\t' + '1/1'
				elif ls[i] == 'NC':
					new_line = new_line + '\t' + './.'

			fo.write(new_line+'\n')
