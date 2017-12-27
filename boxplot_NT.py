#!/usr/bin/env python

'''
Author: Fei Jing
Date: 24/Nov/2017
Usage: python plot_paired.py expression.txt
Function: 
Version: python 2.7

'''
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from scipy.stats import ttest_1samp, wilcoxon, ttest_ind, mannwhitneyu
from math import log

def NT(expF):
        expF_handle=open(expF,'r')

        n = 0
        barcode01T_L =[]
	index01T_L =[]
        barcode11N_L =[]
	index11N_L =[]
	

        for line in expF_handle:
                n += 1
                if n==2:
                        barcode_L = line.split() #contains "barcode11N_L","barcode11N_L"
			
        for barcode in barcode_L:
		index = barcode_L.index(barcode)			#col: index+1
                if "ID" not in barcode:
                        if any (barcode.endswith (m) for m in ["01A","01B","01C","02A","06A"]):
				if barcode not in barcode01T_L:
		                        barcode01T_L.append(barcode)
					index01T_L.append(index)
                        else:
				if barcode not in barcode11N_L:
		                        barcode11N_L.append(barcode)
					index11N_L.append(index)

	return barcode01T_L, barcode11N_L, index01T_L,index11N_L
	
def getValue(expF,index01T_L,index11N_L,target):
	expF_handle=open(expF,'r')
	value01T_L=[]
	value11N_L=[]

	for Line in expF_handle:
		
		Line = Line.rstrip()
		
		ensembl = Line.split("\t",2)[0]
		gene = Line.split("\t",2)[1]
		value = Line.split("\t",2)[2]
		from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
		if ensembl in target:
			Line_L = Line.split()
			for t_id in index01T_L:
				value01T_L.append(Line_L[t_id])			

			for n_id in index11N_L:
				value11N_L.append(Line_L[n_id])
	return value01T_L,value11N_L
		

def Pvalue(L1,L2):
	t,p = mannwhitneyu(L1,L2)			# two-sample wilcoxon
	p="%.3g"%p
	return p


def Summary(barcode01T_L, barcode11N_L,value01T_L,value11N_L):
	print "Sample Barcode\tFPKM-UQ"

	for element in barcode01T_L:
		index = barcode01T_L.index(element)
		value_T = value01T_L[index]
		
		print element +"\t" + str(value_T)
	for element2 in barcode11N_L:
		index2 = barcode11N_L.index(element2)
		value_N = value11N_L[index2]
		print element2 +"\t" + str(value_N)	
	

def boxplot_2box(value01T_L,value11N_L,figureName,target,pvalue):
	## combine these different collections into a list  
	value01T_Log = [log(y,10) for y in value01T_L]
	value11N_Log = [log(x,10) for x in value11N_L]
	data_to_plot = [value01T_Log,value11N_Log]
	countT = len(value01T_L)
	countN = len(value11N_L)
	x_lable1 = "T (n=" + str(countT) + ")"
	x_lable2 = "N (n=" + str(countN) + ")" 
	# Create a figure instance
	fig = plt.figure(1,figsize=(9, 6))

	# Create an axes instance
	ax = fig.add_subplot(111)

	# Create the boxplot,patch_artist: fill with color,
	#bp = ax.boxplot(data_to_plot,labels=[x_lable1,x_lable2], patch_artist=True,widths=(0.5,0.5))
	
	bp = ax.boxplot(data_to_plot,labels=[x_lable1,x_lable2], widths=(0.4,0.4))

	# Labels
	Title = target + " Expression in Colorectal Cancer"
	ax.set_title(Title,fontsize=20)
	ax.set_ylabel('log10(FPKM-UQ)',color='black',fontsize=16)
	
	p_label = "P value="+ str(pvalue)
	at = AnchoredText(p_label,
                  prop=dict(size=12), frameon=False,
                  loc=1,)
	ax.add_artist(at)

	# Grid lines
	#ax.yaxis.grid(True)
	#ax.axvline(linewidth=4, color="r")


	for i in [1,2]:
		y = data_to_plot[i-1]
		x = np.random.normal(i,0.02,len(y))
		plt.plot(x,y,'r.',alpha=0.2)

	
	# Color
	#colors=['pink','lightblue']
	#for patch,color in zip(bp['boxes'],colors):
	#	patch.set_facecolor(color)
	
	# Save the figure
	plt.show()
	fig.savefig(figureName, bbox_inches='tight')
	


if __name__=="__main__":
        expF = sys.argv[1]
	gene = sys.argv[2]
	
	
	control_B2M = "ENSG00000166710"
	if "MAD1L1" in gene:
        	target = "ENSG00000002822.14"
	if "HSF1" in gene:
        	target = "ENSG00000185122.9"
	
	figureName = expF.split(".",1)[0]+"_boxplot_" + gene + ".png"

        barcode01T_L, barcode11N_L,index01T_L,index11N_L = NT(expF)
	value01T_L,value11N_L = getValue(expF,index01T_L,index11N_L,target)
	value01T_L = np.array(value01T_L).astype(np.float)
	value11N_L = np.array(value11N_L).astype(np.float)

	pvalue = Pvalue(value01T_L,value11N_L)
	boxplot_2box(value01T_L,value11N_L,figureName,gene,pvalue)
	Summary(barcode01T_L, barcode11N_L,value01T_L,value11N_L)




