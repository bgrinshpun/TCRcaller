'''
Module for CDR3 identification from TCR sequence data

Boris Grinshpun, 2017

'''

import sys
import time
import re # for parsing the cigar

outputdir=sys.argv[1]
chain=sys.argv[2]
coordinatefile=sys.argv[3]
Vmotiffile=sys.argv[4]
Jmotiffile=sys.argv[5]
codonmapfile=sys.argv[6]
refch=sys.argv[7]


discarded=open(outputdir+'/discarded.TCR'+chain+'.tsv','w')
onematch=open(outputdir+'/VorJonly.TCR'+chain+'.tsv','w')
completematch=open(outputdir+'/VJcalled.TCR'+chain+'.tsv','w')

onematch.write('readID'+'\t'+'CassetteID'+'\n')
discarded.write('readID'+'\t'+'reason'+'\n')
completematch.write('readID'+'\t'+'nucleotide'+'\t'+'cdr3'+'\t'+'Vgene'+'\t'+'Jgene'+'\n')



# Load References

# reference coordinates
coord=open(coordinatefile,'r').readlines()
coord=[i.strip().split('\t') for i in coord]
indices=[1,2,3]
coord=[[x[i] for i in indices] for x in coord][1:]
coord=[[i, int(j), int(k)] for i,j,k in coord]

codonmap=open(codonmapfile,'r').readlines()
codonmap=[line.strip().split() for line in codonmap]
codondict=(dict(codonmap))

Vmotifs=open(Vmotiffile,'r').readlines()
Vmotifs=[i.strip().split('\t') for i in Vmotifs]
indices=[0,1,3]
Vmotifs=[[x[i] for i in indices] for x in Vmotifs][1:]
Vmotifdict={i:(j,int(k)) for i,j,k in Vmotifs}

Jmotifs=open(Jmotiffile,'r').readlines()
Jmotifs=[i.strip().split('\t') for i in Jmotifs]
indices=[0,1,2]
Jmotifs=[[x[i] for i in indices] for x in Jmotifs][1:]
Jmotifdict={i:(j,k) for i,j,k in Jmotifs}

complementdict={'A':'T','T':'A','C':'G','G':'C','N':'N'}

# precompile cigar regex
patt = re.compile(r'(\d+)([A-Z])') # precompiled regex for parsing the CIGAR string
def splitcigar(C,direction="fwd"):
	'''
	Parse the cigar string

	'''
	splitstring=[]
	while C: # capture each instance of the regex in the cigar string
		m = re.match(patt,C)
		if m is None: break
		g1,g2 = m.group(1),m.group(2)
		C = C[len(g1)+len(g2):] # remove captured piece and continue
		if(direction=="rev"):
			splitstring.insert(0,[int(g1),g2]) # reverse cigar if reversed mapping
		else:
			splitstring.append([int(g1),g2]) # append all captured cigar pieces
	
	# identify the starting position of the mapped portion
	matchstring=[(i,clen,cid) for i,(clen,cid) in enumerate(splitstring) if cid=="M"]	
	minpos=0
	maxpos=0;
	if(matchstring[0][0]==0):
		readstart=0
	else:
		minpos=matchstring[0][0]
		readstart=sum([i for i,cid in splitstring[:minpos]])

	maxpos=matchstring[-1][0]
	maplength=sum([clen for clen,cid in splitstring[minpos:(maxpos+1)] if cid != "D"]) # allow for insertions but not deletions
	return readstart,maplength

# translate into all reading frames
def translate(seq,codondict):
	allframes=[[],[],[]]
	for rf in range(3):
		pos=rf
		aatrans=""
		while pos+3<=len(seq):
			s=seq[pos:pos+3]
			if re.search('N',s):
				nextaa='-'
			else:
				nextaa=codondict[s]
			aatrans=aatrans+nextaa
			pos=pos+3
		allframes[rf]=aatrans
	return allframes

def complement(seq):
	return ''.join([complementdict[x] for x in seq])

# translate and clip based on V and J reference to produce final CDR3
def processCDR3(TRV, TRJ, seq, readID):
	myseq=seq

	# Jcassette
	Jm=Jmotifdict[TRJ][1]
	Jm=Jm[0:2]+'.'+Jm[3] # allow a little wobble
	
	fwdAA=translate(myseq, codondict=codondict)
	revAA=translate(complement(myseq[::-1]), codondict=codondict)
	aa_inframe_fwd=[[i,aa] for i, aa in enumerate(fwdAA) if re.search(Jm,aa)]
	aa_inframe_rev=[[i,aa] for i, aa in enumerate(revAA) if re.search(Jm,aa)]
	if len(aa_inframe_fwd)>0 and len(aa_inframe_rev)==0:
		for match in re.finditer(Jm, aa_inframe_fwd[0][1]):
    			pass
		theind=match.start(0)+1
		thecdr3=aa_inframe_fwd[0][1][:theind]
		ref_ind=aa_inframe_fwd[0][0]
		myseq=myseq[ref_ind:((theind)*3+ref_ind)]

	elif len(aa_inframe_fwd)==0 and len(aa_inframe_rev)>0:
                for match in re.finditer(Jm, aa_inframe_rev[0][1]):
                        pass
		myseq=complement(myseq[::-1])
                theind=match.start(0)+1
                thecdr3=aa_inframe_rev[0][1][:theind]
                ref_ind=aa_inframe_rev[0][0]
                myseq=myseq[ref_ind:(theind)*3+ref_ind]
	else:
		
		discarded.write(readID+"\t"+'no_TRJ_RF'+'\n')
		#discarded.write(readID+"\t"+'no_TRJ_RF'+'\t'+TRV+'\t'+TRJ+'\t'+'\n')

		return None

	# Vcassette
	Vfull=Vmotifdict[TRV][0]
	Vpos=Vmotifdict[TRV][1]
	Vm=Vfull[(-Vpos-3):(-Vpos+1)]
        if re.search(Vm,thecdr3):
		for match in re.finditer(Vm, thecdr3):
                        pass
		theind=match.end(0)-1
		thecdr3=thecdr3[theind:]
		myseq=myseq[(theind*3):]
		return [myseq,thecdr3]
	else:
		discarded.write(readID+"\t"+'no_TRV_RF'+'\n')
		#discarded.write(readID+"\t"+'no_TRV_RF'+'\t'+TRV+'\t'+TRJ+'\t'+'\n')
		return None

#### Call TCRs #####
istart=0;
i=0;
start_time = time.time()
read_dict={};
for line in sys.stdin:
    i+=1;
    line=line.strip().split('\t') 
    readID=line[0]
    if(i==1):
    	readIDprev=line[0]
    ch=line[2]
    if(ch==refch):
	start_position=int(line[3])
	if readIDprev==readID:  # collect all read mappings
		read_dict[start_position]=line
	else:
		# find mapping cassettes
		matches=[]
		for key in read_dict.keys():
			matches.append([[key,cas] for cas, p1, p2 in coord if (key>=(p1-25) and key<=(p2+25))])  # include a small buffer in case mapping start is slightly outside the cassette region

		matches=[x for x in  matches for x in x if len(matches)>0]
		allcas=list(set([cas for pos,cas in matches]))

		if len(allcas)>2: # discard multiple matches to a V or a J but keep if one is much higher quality
			check_Vpos=[[pos,cas] for pos,cas in matches if re.match('TR.V',cas)]
			check_Jpos=[[pos,cas] for pos,cas in matches if re.match('TR.J',cas)]
			vmatch=[]; jmatch=[];
			discard=False
			if len(check_Vpos)>2:
				quality=[int(read_dict[k][4]) for k,c in check_Vpos]
				if max(quality)>0:
					vmatch=[[k,c] for (k,c),q in zip(check_Vpos,quality) if q/max(quality)>0.5]
					if len(vmatch)>1:
						discarded.write(readIDprev+'\t'+'too_many_calls'+'\n')
						discard=True
				else:
                                        discarded.write(readIDprev+'\t'+'too_many_calls'+'\n')
                                        discard=True
			if len(check_Jpos)>2:
				quality=[int(read_dict[k][4]) for k,c in check_Jpos]
				if max(quality)>0:
					jmatch=[[k,c] for (k,c),q in zip(check_Jpos,quality) if q/max(quality)>0.5]
					if len(jmatch)>1 and len(vmatch)==1:
						discarded.write(readIDprev+'\t'+'too_many_calls'+'\n')
						discard=True
				else:
                                        discarded.write(readIDprev+'\t'+'too_many_calls'+'\n')
                                        discard=True	

			if discard is False:
				matches=vmatch+jmatch
				allcas=list(set([cas for pos,cas in matches]))
			else:
                		readIDprev=readID;
                		read_dict={};
				continue;

		# definite discardation
		if len(allcas)==1: # single cassette only
			onematch.write(readIDprev+'\t'+matches[0][1]+'\n')
		elif len(allcas)==0: # unmapped
			discarded.write(readIDprev+"\t"+'no_calls'+'\n')
		else:
			if (allcas[0][0:4] != allcas[1][0:4]) :  # one is V while the other is J
				line1=read_dict[matches[0][0]]
				line2=read_dict[matches[1][0]]

				if line1[1] == "16" or line1[1] == "2064" : # check mapping flag
					sequence=complement(line1[9][::-1])
					readstart1,maplength1=splitcigar(line1[5],'rev')
					readstart2,maplength2=splitcigar(line2[5])
				else:
					sequence=line1[9]
					readstart1,maplength1=splitcigar(line1[5])

                                if line2[1] == "16" or line2[1] == "2064" : # check mapping flag
                                        readstart2,maplength2=splitcigar(line2[5],'rev')
                                else:
                                        readstart2,maplength2=splitcigar(line2[5])

				if readstart1<readstart2:
					startpos=readstart1
					endpos=readstart2+maplength2
				else:
					startpos=readstart2
                                        endpos=readstart1+maplength1

				CDR3=sequence[(startpos-1):endpos]

				if re.search('TR.V', allcas[0]):
					TRV=allcas[0]
					Vmap=sequence[readstart1:readstart1+maplength1]
					TRJ=allcas[1]
					Jmap=sequence[readstart2:readstart2+maplength2]
					
					result=processCDR3(TRV, TRJ, CDR3, readIDprev)
					if result is not None:
						completematch.write("\t".join([readIDprev,result[0], result[1], TRV, TRJ,'\n']))
				else:
					TRJ=allcas[0]
					Jmap=sequence[readstart1:readstart1+maplength1]
					TRV=allcas[1]
					Vmap=sequence[readstart2:readstart2+maplength2]
					result=processCDR3(TRV, TRJ, CDR3, readIDprev)
                                        if result is not None:
                                                completematch.write("\t".join([readIDprev,result[0], result[1], TRV, TRJ,'\n']))
			else: 
				discarded.write(readIDprev+'\t'+'same_cassette_type'+'\n')

		readIDprev=readID;
		read_dict={};


    if i%100000==0:
	print "    [Processed lines %s--%s, (%s sec)] " % (istart,i, "{0:.2f}".format(time.time()-start_time))
	istart=(i+1);
	start_time=time.time()
