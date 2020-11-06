#!/usr/bin/python3

import argparse
import subprocess
import os
import json
import glob
import multiprocessing
import tempfile
import sklearn_json as skljson
from sklearn.linear_model import LogisticRegression
from itertools import repeat
import pkg_resources
import sys

insdir = os.path.expanduser('~/POSMM/')
os.system('mkdir -p ' + insdir)
os.system('rm -rf ' + insdir + 'tmp/')
os.system('mkdir -p ' + insdir + 'tmp/')
tmpdir = tempfile.mkdtemp(dir=insdir+'tmp/')
smm = insdir + '/smm'
smmsource = pkg_resources.resource_filename('POSMM','posmmsource/smm.cpp')

# Check if SMM binary exists, otherwise generate.

if os.path.isfile(smm) == False:
	print('SMM not found.  Compiling SMM')
	os.system('g++ -Ofast -std=c++11 ' + smmsource + ' -o ' + smm)
	
os.system('mkdir -p ' + tmpdir)

def linLoad():
	lin = {}
	with open('%s/lin.csv'%(tmpdir)) as infile:
		next(infile)
		for lines in infile:
			lines = lines.rstrip()
			values = lines.split(',')
			taxid = values[0]
			taxid = taxid.split('.')[0]
			lineage = values[1:9]
			lin[taxid] = lineage
	return lin

def chunkSplit(lst,sz):
	for x in range(0,len(lst),sz):
		yield lst[x:x+sz]
			
	
def genomeList(path,cpuc):
	gflist = []
	os.system('mkdir -p %s' %(tmpdir))
	fl = glob.glob(path + '/*.fna')
	if len(fl) == 0:
		print('No valid genomes detected for modelling.  Use --runmode setup to generate a set of genomes for analysis, or build your own using the manual.')
		sys.exit()
	print('Using ' + str(cpuc) + ' threads')
	print(str(len(fl)) + ' Genomes to Compare')
	splitcount = (len(fl) / cpuc) + 1
	if splitcount < 1.0:
		glf = open('%s/G0'%(tmpdir),'w')
		for f in fl:
			glf.write(f + '\n')
		gflist.append('%s/G0' %(tmpdir))
	else:
		splitcount = int(splitcount)
		fl_lst = chunkSplit(fl,splitcount)
		fc = 0
		for f in fl_lst:
			ofile = '%s/G'%(tmpdir) + str(fc)
			glf = open(ofile,'w')
			for f2 in f:
				glf.write(f2 + '\n')
			glf.close()
			gflist.append(ofile)
			fc += 1
	return gflist
	
def mpSet(gflist,fas,order):
	cmdlist = []
	for g in gflist:
		outg = g.replace('G','O')
		cmdlist.append(smm + ' %s %s %s raw %s' %(g,fas,order,outg))
	return cmdlist
		

def cmdFling(cmd):
	os.system(cmd)
	
def topResults():
	grez = []
	srez = []
	outfile = open('%s/ALL.markov'%(tmpdir),'w')
	resultfiles = glob.glob('%s/*.markov'%(tmpdir))
	if len(resultfiles) == 1:
		with open(resultfiles[0]) as infile:
			for lines in infile:
				outfile.write(lines)
	else:
		with open(resultfiles[0]) as infile:
			for lines in infile:
				lines = lines.rstrip()
				values = lines.split('\t')
				grez.append(values[1])
				srez.append(float(values[0]))
		for r in resultfiles[1:]:
			lc = 0
			with open(r) as infile:
				for lines in infile:
					lines = lines.rstrip()
					values = lines.split('\t')
					if float(values[0]) > srez[lc]:
						srez[lc] = float(values[0])
						grez[lc] = values[1]
					lc += 1
		for x in range(len(grez)):
			outfile.write(str(srez[x]) + '\t' + (grez[x].split('/')[-1]).split('.')[0] + '\n')
		outfile.close()
			
	
def rawOutFinal(allmarkovfile,output):
	a = open(output,'w')
	with open(allmarkovfile) as infile:
		for lines in infile:
			lines = lines.rstrip()
			values = lines.split('\t')
			taxid = values[1]
			lineage = glin[taxid][1:]
			a.write('\t'.join(lineage) + '\t' + str(values[0]) + '\n')
	a.close()
		
def countFas(fastafile):
	countlst = []
	nucleo = {'A','T','G','C'}
	count = 0
	with open(fastafile) as infile:
		for lines in infile:
			if lines.startswith('>'):
				count = 0
				continue
			else:
				lines = (lines.rstrip()).upper()
				for nuc in lines:
					if nuc in nucleo:
						count += 1
				countlst.append(count)
	return countlst
						

def conOutFinal(cfile,order):
	a = open(cfile.replace('.tocnf','.possum'),'w')
	try:
			
		specMODEL = skljson.from_json(pkg_resources.resource_filename('POSMM','models/species_' + str(order) + '.json'))
		genMODEL = skljson.from_json(pkg_resources.resource_filename('POSMM','models/genus_' + str(order) + '.json'))
		famMODEL = skljson.from_json(pkg_resources.resource_filename('POSMM','models/family_' + str(order) + '.json'))
		claMODEL = skljson.from_json(pkg_resources.resource_filename('POSMM','models/class_' + str(order) + '.json'))
		ordMODEL = skljson.from_json(pkg_resources.resource_filename('POSMM','models/order_' + str(order) + '.json'))
		phyMODEL = skljson.from_json(pkg_resources.resource_filename('POSMM','models/phylum_' + str(order) + '.json'))
	except:
		print('ERROR: Models Missing.  Aborting')
		sys.exit()
	modlist = [phyMODEL,claMODEL,ordMODEL,famMODEL,genMODEL,specMODEL]
	matchpos = []
	for m in modlist:
		mp = (list(m.classes_)).index('Match')
		matchpos.append(mp)
	with open(cfile) as infile:
		for lines in infile:
			lines = lines.rstrip()
			values = lines.split('\t')
			taxid = values[2]
			lineage = glin[taxid][1:]
			rawscore = float(values[1])
			readlen = int(values[0])
			wrstr = []
			for x in range(len(modlist)):
				cfs = modlist[x].predict_proba([[rawscore,readlen]])[0][matchpos[x]]
				wrstr.append(lineage[x] + ':::' + str(cfs))
			a.write('\t'.join(wrstr) + '\n')
	a.close()

def topResultSplitter(fasfile,cpu):
	print('Splitting Markov output for Confidence Score Generation')
	cntlst = countFas(fasfile)
	fc = 0
	lc = 0
	try:
		with open('%s/ALL.markov'%(tmpdir)) as infile:
			for lines in infile:
				lc += 1
	except OSError:
		print('ERROR: Combined SMM Results Missing.  Aborting.')
		sys.exit()
	if len(cntlst) != lc:
		print('Index Error. Aborting.')
		sys.exit()
	split = int(lc/cpu) + 1
	lc = 0
	with open('%s/ALL.markov'%(tmpdir)) as infile:
		for lines in infile:
			if (lc%split) == False:
				splitout = open('%s/%s.tocnf'%(tmpdir,str(fc)),'w')
				fc += 1
			splitout.write(str(cntlst[lc]) + '\t' + lines)
			lc += 1
			
def confFinalOut(cpu,outfile):
	a = open(outfile,'w')
	for x in range(cpu):
		try:
			with open('%s/%s.possum'%(tmpdir,str(x))) as infile:
				for lines in infile:
					a.write(lines)
		except:
			print('ERROR: Mismatch Results from Split.  Abort')
			sys.exit()
	a.close()
		
def taxLoad(taxmap):
	print('Generating custom taxmap from %s'%(taxmap))
	try:
		taxlist = []
		with open(taxmap) as infile:
			for lines in infile:
				taxlist.append(lines.rstrip())
		return taxlist	
	except OSError:
		print('Problem with custom taxmap file.  Exiting...')
		sys.exit()
		
def taxLoadBest(domain,lin,afile):
	taxlist = []
	specs = {}
	print('Identifying best ' + domain + ' genome RefSeq entries')
	with open(afile) as infile:
		for lines in infile:
			if lines.startswith('#'):
				continue
			lines = lines.rstrip()
			values = lines.split('\t')
			taxid = values[5]
			gcf = values[0]
			if taxid not in lin:
				continue
			if lin[taxid][0] != domain:
				continue
			else:
				species = lin[taxid][6]
				curscore = 0
				if species not in specs:
					specs[species] = [-1,gcf]
				if values[4] == 'representative genome':
					curscore += 1
				if values[10] == 'complete genome':
					curscore += 3
				if values[12] == 'Full':
					curscore += 1
				if curscore > specs[species][0]:
					specs[species][0] = curscore
					specs[species][1] = gcf
	print('Building taxlist')
	for s in specs:
		taxlist.append(specs[s][1])
		print(s)
		print(specs[s][1])
	return taxlist
	

def main():
	# Set up arguments

	parser = argparse.ArgumentParser(prog='POSMM')

	parser.add_argument('--runmode','-r',type=str,default='confidence',choices=['confidence','raw','setup'],help='Run-mode of POSMM.  Using setup will setup a default group of genomes for model generation.  Confidence will output taxa classifications with confidence scores.  Raw will output taxa classifications with raw Markov model / read scores.')
	parser.add_argument('--genomes','-g',type=str,default='Genomes.txt',help='List of Genomes for Read Comparison. Use --runmode setup to download default index or use your own.')	
	parser.add_argument('--order','-o',type=int,default=12,choices=range(10,13),help='Set Markov model order to value between 10-12.  Higher results in better classification accuracy, but longer run-time. Default is 12')
	parser.add_argument('--threads','-t',type=int,default=len(os.sched_getaffinity(0)),help='Set to number of desired threads.  Default uses all available.  Currently Reported: ' + str(len(os.sched_getaffinity(0))))
	parser.add_argument('--fasta','-f',type=str,help='Input Metagenomic Fasta File.')
	parser.add_argument('--gdir','-gd',type=str,help='Location of default genome and taxmap directory if using setup mode.  Defaults to ~/POSMM/Genomes',default='~/POSMM/Genomes')
	parser.add_argument('--gtype','-gt',type=str,help='Type of RefSeq genome set to download.  Default is bacteria, which attempts to download the highest quality genome for all unique bacterial species.',choices=['bacteria','virus','archaea','custom'],default='bacteria')
	parser.add_argument('--output','-out',type=str,help='Output file name')
	parser.add_argument('--taxlist','-tx',type=str,default=insdir+'taxmap.txt',help='Location of taxmap (List of prokaryotic Refseq GCF accessions, one per line) to build database from when using --runmode setup and --gtype custom. Default is ~/POSMM/taxmap.txt.')
	


	args = parser.parse_args()

	if ((args.fasta == None) and (args.runmode != 'setup')):
		parser.error("Input FASTA required for analysis (--fasta/-f) unless using --runmode setup")
		
	if ((args.output == None) and (args.runmode != 'setup')):
		parser.error("Output filename required for analysis (--output/-out) unless using --runmode setup")


	if (args.runmode == 'setup'):
		print('First-Run Mode: Downloading RefSeq Assembly Summary')
		os.system('mkdir -p %s'%(tmpdir))
		os.system('wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt -O %s/assembly_summary_refseq.txt'%(tmpdir))
		os.system('wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -O %s/taxdump.tar.gz'%(tmpdir))
		os.system('tar -xzvf ' + tmpdir + '/taxdump.tar.gz -C ' + tmpdir)
		try:
			cmd = 'ncbitax2lin %s/nodes.dmp %s/names.dmp --output %s/lin.csv.gz' %(tmpdir,tmpdir,tmpdir)
			process = subprocess.Popen(cmd,shell=True)
			process.wait()
			cmd2 = 'gunzip -f %s/lin.csv.gz'%(tmpdir)
			process = subprocess.Popen(cmd2,shell=True)
			process.wait()
		except OSError:
			print('ncbitax2lin not installed.  Try running pip install ncbitax2lin.')
			sys.exit()
		lin = linLoad()
		print('Saving Genomes/Taxmap to ' + args.gdir)
		os.system('mkdir -p ' + args.gdir)
		linjsonfile = args.gdir + '/lineage.json'
		jdf = open(os.path.expanduser(linjsonfile),'w')
		json.dump(lin,jdf)
		jdf.close()
		taxmap = []
		specs = {}
		if args.gtype == 'custom':
			taxmap = taxLoad(args.taxlist)
		elif args.gtype == 'bacteria':
			taxmap = taxLoadBest('Bacteria',lin,tmpdir+'/assembly_summary_refseq.txt')
		elif args.gtype == 'virus':
			taxmap = taxLoadBest('Viruses',lin,tmpdir+'/assembly_summary_refseq.txt')
		elif args.gtype == 'archaea':
			taxmap = taxLoadBest('Archaea',lin,tmpdir+'/assembly_summary_refseq.txt')
		else:
			print('Invalid Library Type (GTYPE ERROR)')
			sys.exit()
		with open('%s/assembly_summary_refseq.txt' %(tmpdir)) as infile:
			for lines in infile:
				if lines.startswith('#'):
					continue
				values = lines.split('\t')
				if args.gtype == 'custom':
					taxid = values[5]
					gcf = values[0]
					if gcf in taxmap:
						link = values[19]
						os.system('wget ' + link + '/' + link.split('/')[-1] + '_genomic.fna.gz' + ' -O ' + os.path.expanduser(args.gdir) + '/' + taxid + '.' + gcf + '.fna.gz')
				elif args.gtype == 'bacteria':
					taxid = values[5]
					gcf = values[0]
					if gcf in taxmap:
						link = values[19]
						os.system('wget ' + link + '/' + link.split('/')[-1] + '_genomic.fna.gz' + ' -O ' + os.path.expanduser(args.gdir) + '/' + taxid + '.' + gcf + '.fna.gz')
		os.system('gunzip -f ' + os.path.expanduser(args.gdir) + '/*.gz')
	else:
		print('Loading lineage from JSON')
		try:
			jdf = open(os.path.expanduser(args.gdir) + '/lineage.json')
		except OSError:
			print('No valid lineage file found.  Either generate your own lineage.json, or use --runmode setup to build directly from Refseq')
			sys.exit()
		global glin
		glin = json.load(jdf)
		jdf.close()
		print('Reading Genome Directory')
		gflist = genomeList(os.path.expanduser(args.gdir),args.threads)
		mplist = mpSet(gflist,args.fasta,str(args.order))
		pool = multiprocessing.Pool(args.threads)
		pool.map(cmdFling,mplist)
		topResults()
		
		if (args.runmode == 'raw'):
			rawOutFinal(tmpdir + '/ALL.markov',args.output)
		elif (args.runmode == 'confidence'):
			topResultSplitter(args.fasta,args.threads)
			cfiles = glob.glob('%s/*.tocnf'%(tmpdir))
			pool = multiprocessing.Pool(args.threads)
			pool.starmap(conOutFinal,zip(cfiles,repeat(args.order)))
			confFinalOut(args.threads,args.output)
	os.system('rm -rf ' + tmpdir)


if __name__=="__main__":
	main()

					
