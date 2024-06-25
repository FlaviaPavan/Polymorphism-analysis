#! usr/bin/python3
import sys
import re
import numpy as np
import pandas as pd

## Calculate the pi statistics ##
def calculate_statistics(lstpi, dege, fout):
	filename="file_rawpi.txt"
	with open(filename, 'a') as f: 
    		for i in lstpi: 
        		name= dege
        		f.write('%s\t' % (i))
        		f.write(name+'\n')
	f.close()
	
	nbsite=len(lstpi)
	l= pd.Series(lstpi)
	l = l.apply(lambda x: float(x))
	average = l.mean(skipna=True)
	standeviation= np.std(l)
	fout.write(f'{dege}\t{average}\t{standeviation}\t{nbsite}\n')

## Do a dictionary with positions associated to pi value ##
def load_pi_values(pi_file):
	dico_pi={}
	with open(pi_file, 'r') as fin:
		for line in fin:			
			line=line.strip("\n")
			if line.startswith('CHROM'):
				continue
			else:
				tab=line.split('\t')
				ID=f'{tab[0]}:{tab[1]}'
				dico_pi[ID]=tab[2]
		return dico_pi


def process_positions(input_file, dico_pi, line_prefix):	
	lst=[]
	lst2=[]
	lst3=[]
	
	with open(input_file, 'r') as fin:
		for line in fin:			
			line=line.strip("\n")
			if line.startswith(line_prefix):
				tab=line.split(" ")
				position=f'{tab[0]}:{tab[1]}'
				degenerated=tab[5]
				
				match=re.search('3',degenerated)
				if match:
					if position in dico_pi:						
						pi=dico_pi[position]
						lst.append(pi)
					else:
						print(position+' no pi data')
				
				match2=re.search('2',degenerated)
				if match2:					
					if position in dico_pi:
						pi=dico_pi[position]
						lst2.append(pi)
					else:
						print(position+' no pi data')
				
				match3=re.search('4',degenerated)		
				if match3:
					if position in dico_pi:
						pi=dico_pi[position]
						lst3.append(pi)
					else:
						print(position+' no pi data')
				else:
				
					continue
	return lst, lst2, lst3

def main():
	degenerated_position_file = sys.argv[1] # Degenerated target position file 
	pi_file = sys.argv[2] # File with calculated pi for each positions
	line_prefix = sys.argv[3] # Common prefix of all chromosome names in degenerated position file
	
	dico_pi = load_pi_values(pi_file)
	
	with open('pi_degenerated_position_stat.csv', 'w') as fout:
		fout.write('degenerated_position\taverage_pi\tstandard_deviation\tnb_values\n')
		lst, lst2, lst3 = process_positions(degenerated_position_file, dico_pi, line_prefix)
		calculate_statistics(lst, '0fold', fout)
		calculate_statistics(lst2, '2fold', fout)
		calculate_statistics(lst3, '4fold', fout)

if __name__ == "__main__":
	main()
    	

