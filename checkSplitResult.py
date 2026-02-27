#!/usr/bin/env python
#coding: utf8
###### Import Modules
from __future__ import print_function
import sys,os
import glob
import gzip
import time
from collections import defaultdict
from multiprocessing import Pool, freeze_support

###### Version and Date
prog_version = '0.1.0'
prog_date = '2020-04-21'

usage = '\n%s <fqDir>'%sys.argv[0]

def pairCheck(fq1, fq2):
	fh1 = gzip.open(fq1)
	fh2 = gzip.open(fq2)
	cu = 0
	flag = 0
	while True:
		rec1 = [fh1.readline().decode() for i in range(4)]
		rec2 = [fh2.readline().decode() for i in range(4)]
		if not rec1[0]:
			break
		if rec1[0].split('/')[0] != rec2[0].split('/')[0]:
			print('paired reads have different names: %s %s in %s'%(rec1[0].split('/')[0][1:], rec2[0].split('/')[0], os.path.basename(fq1)))
			flag = 1
			break
		cu += 1
	fh1.close()
	fh2.close()
	return [cu, flag]

def getReadNum(txtFile):
	cu = 0
	for k in open(txtFile):
		if k.startswith('#ReadNum'):
			cu = int(k.strip().split()[-1])
			break
	return cu

def checkPE(fqDir):
	#print (fqDir)
	#print (glob.glob(os.path.join(fqDir, '*.fq.gz')))
	fq = [k for k in glob.glob(os.path.join(fqDir, '*.fq.gz')) if os.path.getsize(k) > 10000][0]
	tag = gzip.open(fq).readline().decode().split('/')
	fqList = []
	PE = True
	if len(tag) > 1:
		fqList = glob.glob(os.path.join(fqDir, '*_1.fq.gz'))
	else:
		fqList = glob.glob(os.path.join(fqDir, '*.fq.gz'))
		PE = False
	return (fqList, PE)

def seCheck(fq):
	fh = gzip.open(fq)
	cu = 0
	while True:
		rec = [fh.readline().decode() for i in range(4)]
		if not rec[0]:
			break
		cu += 1
	return cu

def main():
	import argparse
	ArgParser = argparse.ArgumentParser(usage = usage)
	ArgParser.add_argument("-v", "--version", action="version", version=prog_version)
	ArgParser.add_argument("-r", "--rawFqDir", action="store", dest="rawFqDir", default=None, help="Input directory path of raw fq file.")

	(para, args) = ArgParser.parse_known_args()
	if len(args) != 1:
	    ArgParser.print_help()
	    print ("\nERROR: The parameters number is not correct!") #file=sys.stderr
	    sys.exit(1)

	fqDir = args[0]
	fqList,PE = checkPE(fqDir)
	if PE:
		rec = []
		p = Pool(10)
		for fq in fqList:
			rec.append(p.apply_async(pairCheck, args=(fq, fq.replace('_1.fq.gz', '_2.fq.gz'))))
		p.close()
		p.join()
		tmp = [0,0]
		for k in rec:
			k = k.get()
			tmp[0] += k[0]
			tmp[1] += k[1]
		if tmp[1] > 0:
			print ('ERROR!: paired reads have different name')
			sys.exit(1)
		if para.rawFqDir:
			fq = glob.glob(os.path.join(para.rawFqDir, '*_1.fq.gz'))[0]
			txtFile = fq.replace('fq.gz', 'fq.fqStat.txt')
			if txtFile:
				rawNum = getReadNum(txtFile)
				if rawNum != tmp[0]:
					print ('ERROR!: The total reads number after splitBarcode is not match the raw reads number: %d %d'%(tmp[0], rawNum))
					sys.exit(1)
		print ('PASS: no exceptions')
	else:
		if not para.rawFqDir:
			print ('For SE data inspection, you should add argument "-r" to appoint the rawFqDir')
			sys.exit(1)
		rec = []
		p = Pool(10)
		for fq in fqList:
			rec.append(p.apply_async(seCheck, args=(fq,)))
		p.close()
		p.join()
		num = 0
		for k in rec:
			num += k.get()
		fq = glob.glob(os.path.join(para.rawFqDir, '*.fq.gz'))[0]
		txtFile = fq.replace('fq.gz', 'fq.fqStat.txt')
		if txtFile:
			rawNum = getReadNum(txtFile)
			if rawNum != num:
				print ('ERROR!: The total reads number after splitBarcode is not match the raw reads number: %d %d'%(num, rawNum))
				sys.exit(1)
		print ('PASS: no exceptions')

if __name__ == "__main__":
	freeze_support()
	main()
