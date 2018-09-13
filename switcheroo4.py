# check line format:
# MismatchRefBase	1:565508:G-A/.

# VCF line format:
# Comments/header start with #
# REF and ALT are field number 4 and 5

from sys import argv
import os.path as op

vcf_file_dir = argv[1]   # ~/science/Cairo/checkVCF-20140116/from_scratch_vcf/ ()
check_file_dir = argv[2] # ~/science/Cairo/checkVCF-20140116/from_scratch_vcf/check_files/

CHR = argv[3]
chr = 'chr%s' % CHR
vcf_file = op.join(vcf_file_dir,'%s_calls.vcf' % chr)
check_file = op.join(check_file_dir,'%s.check.ref' % chr)
out_file = 'Cairo/checkVCF-20140116/from_scratch_vcf1/switched/%s_switched4.vcf' % chr

with open(check_file,'r') as f:
	ref_dict = {l.split(':')[1] : l.split(':')[2][0] for l in f}

with open(vcf_file,'r') as fi:
	with open(out_file,'w') as fo:
		for line in fi:
			# skip the header
			if line.startswith('#'):
				fo.write(line)
				continue
			ls1 = line.strip('\n')
			ls = ls1.split('\t')
			pos = ls[1]

      		#  Must flip the calls for particular instances
			if pos in ref_dict:

				if ref_dict[pos]=='C':
					# If the alt is the same, also keep the same ref but switch places.
					if ls[4] == 'C':
						ls[4] = ls[3]

					elif ls[4]=='G':
						if ls[3]=='A':
							ls[4]='T'

						elif ls[3]=='T':
							ls[4]='A'

					# Switch allele calls
					for x in range(9, len(ls)):
						if ls[x] == '0/0':
							ls[x] = '1/1'
						elif ls[x] == '1/1':
							ls[x] = '0/0'

				if ref_dict[pos]=='G':
					# If the alt is the same, also keep the same ref but switch places.
					if ls[4] == 'G':
						ls[4] = ls[3]

					elif ls[4]=='C':
						if ls[3]=='A':
							ls[4]='T'

						elif ls[3]=='T':
							ls[4]='A'

					# Switch allele calls
					for x in range(9, len(ls)):
						if ls[x] == '0/0':
							ls[x] = '1/1'
						elif ls[x] == '1/1':
							ls[x] = '0/0'

				if ref_dict[pos]=='T':
					if ls[3]=='A':
						if ls[4]=='G':
							ls[4]='C'
						elif ls[4]=='C':
							ls[4]='G'
						elif ls[4]=='T':
							ls[4]=ls[3]
							# Switch allele calls
							for x in range(9, len(ls)):
								if ls[x] == '0/0':
									ls[x] = '1/1'
								elif ls[x] == '1/1':
									ls[x] = '0/0'

				if ref_dict[pos]=='A':
					if ls[3]=='T':
						if ls[4]=='G':
							ls[4]='C'
						elif ls[4]=='C':
							ls[4]='G'
						elif ls[4]=='A':
							ls[4]=ls[3]
							# Switch allele calls
							for x in range(9, len(ls)):
								if ls[x] == '0/0':
									ls[x] = '1/1'
								elif ls[x] == '1/1':
									ls[x] = '0/0'

				# Make the ref match the ref for all instances.
				ls[3]= ref_dict[pos]


				fo.write('\t'.join(ls)+'\n')
				continue
			else:
				fo.write(line)
