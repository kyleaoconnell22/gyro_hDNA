import os
import sys
import subprocess as sp
import shutil
import argparse
import random
'''
Script will subsample each replicate type from min to max size and then calculate genetic diversity using STACKS module for the number of 
reps specified by bootstrap parameter

need replicate_list.txt which needs to list the replicates like this:
supernatant
pellet
formalin-fixed
RAD

And these need to correspond to existing files called:
supernatant_keep.txt
pellet_keep.txt 
etc. Which need to list all the samples for this replicate-type, and match the names in the vcffile

Kyle O'Connell
kyleaoconnell22@gmail.com
'''

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--in_dir", required=True, help="REQUIRED: The full path to a directory which contains the input alignment files.")
    parser.add_argument("-b", "--boots", required=True, help="REQUIRED: Number of bootstraps.")
    parser.add_argument("-min", "--min_ind", required=True, help="REQUIRED: Min number of individuals to randomly sample = smallest dataset")
    parser.add_argument("-max", "--max_ind", required=True, help="REQUIRED: Max number of individuals to randomly sample = largest dataset")
    parser.add_argument("-v", "--in_vcf", required=True, help="REQUIRED: VCF with full dataset")
    parser.add_argument("-R", "--replicates", required=True, help="REQUIRED: list of replicate types coresponding to existing keep files, see documentation")
    return parser.parse_args()
	
def subsample_loop(in_dir,boots,min_ind,max_ind,in_vcf,replicates):
	os.chdir(in_dir)
	
	rep_list = [i.strip() for i in open(replicates)]

	for rep in rep_list:
		#downsample vcf
		vcf	= in_vcf
		keep_file = rep+'_keep.txt'			
		outfile = rep+'_sumstats_summary_combined.txt'
		fh_out = open(outfile,'w')
		fh_out.write('# Pop ID'+'\t'+'Private'+'\t'+'Sites'+'\t'+'Variant_Sites'+'\t'+"Polymorphic_Sites"+'\t'+"%Polymorphic_Loci"\
		+'\t'+"Num_Indv"+'\t'+"Var_N"+'\t'+"StdErr_N"+'\t'+"P"+'\t'+"Var_P"+'\t'+"StdErr_P"+'\t'+"Obs_Het"+'\t'+"Var_OHe"+'\t'+"StdErr_OHe"\
		+'\t'+"Obs_Hom"+'\t'+"Var_OHo"+'\t'+"StdErr_OHo"+'\t'+"Exp_Het"+'\t'+"Var_EHe"+'\t'+"StdErr_EHe"+'\t'+"Exp_Hom"+'\t'+"Var_EHo"+'\t'+"StdErr_EHo"\
		+'\t'+"Pi"+'\t'+"Var_Pi"+'\t'+"StdErr_Pi"+'\t'+"Fis"+'\t'+"Var_Fis"+'\t'+"StdErr_Fis"+"\n")

		for i in range(int(boots)): #create #boots loops
			j= random.randrange(int(min_ind),int(max_ind))
			#print i,j
			out = vcf.split('.')[0]+'_temp.'+str(i)+'_'+str(j)
			#call vcftools to subsample randomly, max-indv randomly thins
			vcf_sub = "vcftools --vcf {0} --keep {1} --max-indv {2} --recode --out {3} ".format(vcf,keep_file,j,out)
			proc_sub = sp.call(vcf_sub,shell=True)
				
			pops_in = out+'.recode.vcf'
			populations_call = "populations -V {0} -O .".format(pops_in)
			proc_populations = sp.call(populations_call,shell=True)
					
					
			sumstats = pops_in.split('.')[0] + '.'+ pops_in.split('.')[1]+'.recode' + '.p.sumstats_summary.tsv'
			last_line = open(sumstats, "r").readlines()[-1]
					
			#write last line to file
			fh_out.write(last_line)
																
			#remove log files
			for filetype in os.listdir('.'):
					if filetype.endswith('.log') or filetype.endswith('.tsv') or 'temp' in filetype:
						os.remove(filetype)

def rename_outs(in_dir,replicates):
	os.chdir(in_dir)
	fh_out = open('sumstats_combined.txt','w')
	fh_out.write('Replicate'+'\t'+'Private'+'\t'+'Sites'+'\t'+'Variant_Sites'+'\t'+"Polymorphic_Sites"+'\t'+"%Polymorphic_Loci"\
	+'\t'+"Num_Indv"+'\t'+"Var_N"+'\t'+"StdErr_N"+'\t'+"P"+'\t'+"Var_P"+'\t'+"StdErr_P"+'\t'+"Obs_Het"+'\t'+"Var_OHe"+'\t'+"StdErr_OHe"\
	+'\t'+"Obs_Hom"+'\t'+"Var_OHo"+'\t'+"StdErr_OHo"+'\t'+"Exp_Het"+'\t'+"Var_EHe"+'\t'+"StdErr_EHe"+'\t'+"Exp_Hom"+'\t'+"Var_EHo"+'\t'+"StdErr_EHo"\
	+'\t'+"Pi"+'\t'+"Var_Pi"+'\t'+"StdErr_Pi"+'\t'+"Fis"+'\t'+"Var_Fis"+'\t'+"StdErr_Fis"+"\n")
	
	rep_list = [i.strip() for i in open(replicates)]
	for rep in rep_list:
		temp = rep + '_sumstats_summary_combined.txt'
		fh_temp = open(temp,'r')
		for line in fh_temp:
			if line.startswith('#'):
				pass
			else:
				line=line.strip()
				line=line.split('\t')
				fh_out.write(rep+'\t')
				for i in line[1:]:
					fh_out.write(i+'\t')
				fh_out.write('\n')

def main():
	#define the arguments
	args = get_args()
	subsample_loop(args.in_dir,args.boots,args.min_ind,args.max_ind,args.in_vcf,args.replicates)
	rename_outs(args.in_dir,args.replicates)
	
if __name__ == '__main__':
    main()	
	
	
	
	
	