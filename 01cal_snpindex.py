#!/usr/bin/env python

#eg: python 01cal_snpindex.py lobed_unlobed_pool_A01.vcf confidence_interval.txt lobed_unlobed_A01_index.xls lobed unlobed > lobed_unlobed_A01_index.log
import sys

i1 = open(sys.argv[1]) #.vcf
i2 = open(sys.argv[2]) #confidence_interval_file
o=open(sys.argv[3],'w')
name1 = sys.argv[4] #dominant pool name ####lobed
name2 = sys.argv[5] #recessive pool name ####unlobed

confi = {}
for line in i2:
	if 'DEPTH' in line:continue
	line = line.strip().split()
	confi[int(line[0])] = float(line[2])

o.write('#CHRMO\tPOS\tREF\tALT\tAD-pool-'+name1+'\tpool-snpindex-'+name1+'\tAD-pool-'+name2+'\tpool-snpindx-'+name2+'\tdelta\tconfidence_interval\tAD-parent-'+name1+'\tAD-parent-'+name2+'\n')
for li in i1:
	li=li.strip()
	if '#' not in li:
		line = li.split()
		p_domi = line[-3].split(':')
		p_rece = line[-4].split(':')
		f_domi = line[-1].split(':')
		f_rece = line[-2].split(':')
		
		
		if '|' not in line[-1] and '|' not in line[-2] and '|' not in line[-3] and '|' not in line[-4] and ',' not in line[3] and ',' not in line[4] and float(line[5]) >= 20 and p_domi[0]!='./.' and p_rece[0]!='./.' and f_domi[0]!='./.' and f_rece[0]!='./.' and 3 <= float(p_rece[2]) <= 75 and 3 <= float(p_domi[2]) <= 75 and 3 <= float(f_rece[2]) <= 75 and  3 <= float(f_domi[2]) <= 75 and ((p_rece[0] == '0/0' and p_domi[0] == '1/1') or (p_domi[0] == '0/0' and p_rece[0] == '1/1')):
			F2_domi_AD=f_domi[1]
			f_domi_AD=F2_domi_AD.split(',')
			F2_rece_AD=f_rece[1]
			f_rece_AD=F2_rece_AD.split(',')
			if p_rece[0] == '1/1':
				F2_domi_index=(float(f_domi_AD[0]))/(float(f_domi_AD[0])+float(f_domi_AD[1]))
				F2_rece_index=(float(f_rece_AD[0]))/(float(f_rece_AD[0])+float(f_rece_AD[1]))
			if p_rece[0] == '0/0':
				F2_domi_index=(float(f_domi_AD[1]))/(float(f_domi_AD[0])+float(f_domi_AD[1]))
				F2_rece_index=(float(f_rece_AD[1]))/(float(f_rece_AD[0])+float(f_rece_AD[1]))
			delta=F2_domi_index-F2_rece_index
			depth=int(f_domi[2])+int(f_rece[2])/2
			con=confi[depth]
			P_domi_AD=p_domi[1]
			P_rece_AD=p_rece[1]
			o.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\n'.format(line[0],line[1],line[3],line[4],F2_domi_AD,F2_domi_index,F2_rece_AD,F2_rece_index,delta,con,P_domi_AD,P_rece_AD))

		else:
			print li

i1.close()
i2.close()
o.close()
