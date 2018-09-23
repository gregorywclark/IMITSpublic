#!/usr/bin/env python

import numpy as np
import pandas as pd
import zlib,bz2
import sys
from math import log
import itertools
import seaborn as sns
from itertools import combinations
from itertools import product
from itertools import combinations_with_replacement
from itertools import islice
from collections import OrderedDict
from Bio.Seq import Seq
from itertools import islice
import statsmodels.api as sm
from statsmodels.formula.api import ols
from collections import defaultdict
import MySQLdb as mdb
import glob
import re

def Reverse(seq):
        orig=Seq(seq)
        rev=str(orig.reverse_complement())
        return rev

def window(seq, n):
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result    
    for elem in it:
        result =result[1:] + (elem,)
        yield result


matchpam=re.compile("[ACGTN]{21}GG",flags=re.X|re.I)
def trifonov(string):
	prod=1e9	##distribution looks normal/gaussian, but this scale gives better 'looking' numbers. 
			##e.g. from 1.248e-9 to 1.24. All numbers end up as e-9, rescaled to >0 and < 10.
	for n in range(2,len(string)+1):
		if n >=18: argmin=4
		else:	argmin=min(n,3)
			
		cd=set(map(lambda q: "".join(q),list(window(string,n))))
		gg=list(map(lambda q: "".join(q),list(window(string,argmin))))
		
		prod*=len(cd)/float(len(gg))
	return prod

def entropy(sequence):
	lsq=list(window(sequence,2))
	dNs=["".join(item) for item in lsq]
	uniquedN=list(set(dNs))
	if len(uniquedN) <= 1:
		return 0
	probs = [dNs.count(item)/float(len(lsq)) for item in uniquedN] 
	ent = 0.
	# Compute standard entropy.
	for i in probs:
		ent -= i * log(i,len(lsq))
	return ent

def getSequences(df):
	Sequences=defaultdict(list)
	for idx, lclrow in df.iterrows():
		#embryos=[lclrow['#Embryos Survived'],lclrow['#Embryos Survived to 2 cell stage'],lclrow['#Embryos Transfered']]
		if str(lclrow['gRNA Sequence (+ strand)']) != 'nan':
			Sequences[str(lclrow['Mi Attempt URL'])].append(str(lclrow['gRNA Sequence (+ strand)']))
		if str(lclrow['gRNA Sequence (+ strand).1']) != 'nan':
			Sequences[str(lclrow['Mi Attempt URL'])].append(str(lclrow['gRNA Sequence (+ strand).1']))
		if str(lclrow['gRNA Sequence (+ strand).2']) != 'nan':
			Sequences[str(lclrow['Mi Attempt URL'])].append(str(lclrow['gRNA Sequence (+ strand).2']))
		if str(lclrow['gRNA Sequence (+ strand).3']) != 'nan':
			Sequences[str(lclrow['Mi Attempt URL'])].append(str(lclrow['gRNA Sequence (+ strand).3']))
	return Sequences


def getTwoBits():
	import twobitreader as tbr
	basedir="/home/clarkg/MouseGenomes/MM10/"
	files=glob.glob(basedir+"*.2bit")
	genome={}
	for f in files:
		chrm=f.split("/")[-1].split(".")[4]
		genome[chrm]=tbr.TwoBitFile(f)
	return genome

def Clean_gRNA(df,genome):
	import twobitreader as tbr

	validgRNA=[]
	numbergRNA=[]
	missing=0
	found=0
	gFound=0
	gReverse=0
	gMissing=0
	MSEQ={}
	MSEQinfo={}
	writeGrna=open("crisporExtended.fasta",'w')
	for idx, lclrow in df.iterrows():
		Sequences=[]
		found=0
		#embryos=[lclrow['#Embryos Survived'],lclrow['#Embryos Survived to 2 cell stage'],lclrow['#Embryos Transfered']]
		if str(lclrow['gRNA Sequence (+ strand)']) != 'nan':
			Sequences.append(str(lclrow['gRNA Sequence (+ strand)']))
		if str(lclrow['gRNA Sequence (+ strand).1']) != 'nan':
			Sequences.append(str(lclrow['gRNA Sequence (+ strand).1']))
		if str(lclrow['gRNA Sequence (+ strand).2']) != 'nan':
			Sequences.append(str(lclrow['gRNA Sequence (+ strand).2']))
		if str(lclrow['gRNA Sequence (+ strand).3']) != 'nan':
			Sequences.append(str(lclrow['gRNA Sequence (+ strand).3']))
		Sequences=list(Sequences)
		othergRNA=[]
		URLID=str(lclrow['Mi Attempt URL']).split("/")[-1]
		ttlSeq=str(len(Sequences))
		for x in range(len(Sequences)):
			seqA=Sequences[x].strip()[3:]
			seqB=Sequences[x].strip()[:20]
			badlength=False
			if len(Sequences[x]) > 23:
				badlength=True
				print"ERROR,passing"
			if x == 0:
				suffix=""
			else:
				suffix="."+str(x)
			chromosome=str(lclrow['Chromosome (+ strand)'+suffix]).upper()
			start=int(lclrow['Start Co-od'+suffix])-1
			end=int(lclrow['End Co-od'+suffix])
			revseq=Reverse(Sequences[x])
			NamedGrna=chromosome+":"+str(start)+"-"+str(end)
			othergRNA.append(NamedGrna)
			try:
				genomicSeq=genome[chromosome][chromosome][start:end]
			except:
				print chromosome,chromosome,start,end	
				genomicSeq=""

			if Sequences[x] == genomicSeq and not badlength:
				#extendedSeq=genome[chromosome][chromosome][(start):(end)]
				grnaSeq=genome[chromosome][chromosome][(start):(end)]
				RgrnaSeq=Reverse(genome[chromosome][chromosome][(start):(end)])
				#writeGrna.write(",".join([URLID,chromosome,str(start),str(end),Sequences[x],"forward"])+"\n")
				if matchpam.search(grnaSeq):
					extendedSeq=genome[chromosome][chromosome][(start-30):(end+47)]
					writeGrna.write(">"+URLID+"_"+str(x+1)+"of"+ttlSeq+"_"+chromosome+":"+str(start-30)+"-"+str(end+47)+"_forward\n")
					writeGrna.write(extendedSeq+"\n")
				elif matchpam.search(RgrnaSeq):
					extendedSeq=Reverse(genome[chromosome][chromosome][(start-47):(end+30)])
					#grnaSeq=revesrse(grnaSeq)
					writeGrna.write(">"+URLID+"_"+str(x+1)+"of"+ttlSeq+"_"+chromosome+":"+str(start-47)+"-"+str(end+30)+"_forward\n")
					writeGrna.write(extendedSeq+"\n")
				found+=1
				gFound+=1
					
			elif revseq == genomicSeq and not badlength:
				#extendedSeq=genome[chromosome][chromosome][(start):(end)]
				#writeGrna.write(",".join([URLID,chromosome,str(start),str(end),Sequences[x],"reverse"])+"\n")
				#writeGrna.write(">"+URLID+"_"+chromosome+":"+str(start)+"-"+str(end)+"_reverse\n")
				#writeGrna.write(genomicSeq+"\t")
				grnaSeq=genome[chromosome][chromosome][(start):(end)]
				RgrnaSeq=Reverse(genome[chromosome][chromosome][(start):(end)])
				#writeGrna.write(",".join([URLID,chromosome,str(start),str(end),Sequences[x],"forward"])+"\n")
				if matchpam.search(grnaSeq):
					extendedSeq=genome[chromosome][chromosome][(start-30):(end+47)]
					writeGrna.write(">"+URLID+"_"+str(x+1)+"of"+ttlSeq+"_"+chromosome+":"+str(start-30)+"-"+str(end+47)+"_forward\n")
					writeGrna.write(extendedSeq+"\n")
				elif matchpam.search(RgrnaSeq):
					extendedSeq=Reverse(genome[chromosome][chromosome][(start-47):(end+30)])
					#grnaSeq=revesrse(grnaSeq)
					writeGrna.write(">"+URLID+"_"+str(x+1)+"of"+ttlSeq+"_"+chromosome+":"+str(start-47)+"-"+str(end+30)+"_forward\n")
					writeGrna.write(extendedSeq+"\n")
				gReverse+=1
				found+=1
			else:
				#writeGrna.write(",".join([URLID,chromosome,str(start),str(end),Sequences[x],"MISSING"])+"\n")
				MSEQ[chromosome+":"+str(start)+"-"+str(end)]=Sequences[x]
				othergRNA.remove(NamedGrna)

				MSEQinfo[chromosome+":"+str(start)+"-"+str(end)]=[othergRNA,str(lclrow['Gene Marker Symbol']),Sequences[x],str(lclrow['Mi Attempt URL']),str(lclrow['Production Centre'])]
				oops=[str(lclrow['Mi Attempt URL']),str(lclrow['Production Centre']),Sequences[x],chromosome+":"+str(start)+"-"+str(end)]
				#print ",".join(oops)
				gMissing+=1
				
		numbergRNA.append(len(Sequences))		
		if len(Sequences) == found:
			validgRNA.append(1)
		elif found > 0:
			validgRNA.append(2)
		elif found  ==0  and len(Sequences):
			validgRNA.append(0)
		else:
			validgRNA.append(3)
	#print missing,missing+found
	writeGrna.close()
	##We output as a fasta file
	kl=open("MISSING_SEQUENCES.fa",'w')
	for k,j in MSEQ.iteritems():
		kl.write(">"+k+"\n")
		kl.write(j+"\n")
	kl.close()
	##We output all info simply as a Pkl file, 
		## contains location,sequence,Mi attempt URL, Production centre, etc
	from cPickle import dump,load
	pl=open('INFO_MISSING.pkl','wb')
	dump(MSEQinfo,pl)
	pl.close()
	#print gFound,gReverse,gMissing,gFound+gReverse+gMissing
	df['#Found gRNA'] = pd.Series(validgRNA, index=df.index)
	df['#gRNA'] = pd.Series(numbergRNA, index=df.index)


def cutSizeDet(sorted_cuts):
	cutsize=[]
	if len(sorted_cuts) == 1:
		cutsize=1
		cuttype="1"
	elif len(sorted_cuts) == 2:
		cutsize=abs(sorted_cuts[0]-sorted_cuts[1])
		cuttype="2"
	elif len(sorted_cuts) == 4:
		Lft=[sorted_cuts[0],sorted_cuts[1]]
		Rgt=[sorted_cuts[2],sorted_cuts[3]]
		cutOuter=abs(sorted_cuts[3]-sorted_cuts[0])
		cutInner=abs(sorted_cuts[2]-sorted_cuts[1])
		###
		cutsize=np.mean([cutOuter,cutInner])
		cuttype="4"
	elif len(sorted_cuts) > 2:
		#print len(sorted_cuts),sorted_cuts
		leftD=[]
		leftC=[]
		rightC=[]
		rightD=[]
		for j in range(0,(len(sorted_cuts)/2)+1):
			leftD.append(sorted_cuts[j+1]-sorted_cuts[j])
			leftC.append(str(sorted_cuts[j])+"_"+str(sorted_cuts[j+1]))
		for j in reversed(range(len(sorted_cuts)/2+1,len(sorted_cuts))):
			rightD.append(sorted_cuts[j]-sorted_cuts[j-1])
			rightC.append(str(sorted_cuts[j-1])+"_"+str(sorted_cuts[j]))

		jntD=leftD+rightD
		jntC=leftC+rightC
		mDist=max(jntD)
		
		for a,b in zip(jntD,jntC):
			if a == mDist:
				flankLeft,flankRight=map(lambda i: int(i),b.split("_"))
				break
		Li=sorted_cuts.index(flankLeft)
		Ri=sorted_cuts.index(flankRight)
		Lft=sorted_cuts[:(Li+1)]
		Rgt=sorted_cuts[Ri:]
		cutOuter=abs(Rgt[-1]-Lft[0])
		cutInner=abs(flankRight - flankLeft)
		###
		cutsize=np.mean([cutOuter,cutInner])
		cuttype="other"
	else:
		print 20*"\tWTF\n"

	return cutsize,cuttype
		

if __name__ == "__main__":
	import warnings
	warnings.filterwarnings("ignore")
	# We get a future warning about a pandas core util that we don't even use

	
	#data=pd.ExcelFile("AllCas9IMPC_20180225.xlsx")
	#data=pd.ExcelFile("iMITS_Cas9_Data_20180503.xlsx")
	data=pd.ExcelFile("IMITS_data_2018_09_07.xlsx")
	Attempts=data.parse("Mouse Production")

	mm10=getTwoBits()

	Clean_gRNA(Attempts,mm10)
