import re
import subprocess
import ast
import pickle
from .constants import *
# Get condition lists
import itertools
def getranklist(jsonfile):
    def getrank(dt):
        listolist=[data[x] for x in data['rankorder']]
        return listolist
    data=json.load(open(jsonfile))
    ranktuple=list(itertools.product(*getrank(data)))
    ranklist=[list(x) for x in ranktuple]
    return ranklist

##### Helper functions############3
#  Saves & Load obj from pickled file. 
#######################
def pickle_dump(obj,picklefile):
    with open (picklefile,'wb') as f:
        pickle.dump(obj,f)
def pickle_load(picklefile):
    with open (picklefile,'rb') as f:
        obj = pickle.load(f)     
        return obj
def mapint(string):
    return list(map(int,string.split(",")))
# @profile
def subprocess_cmd(command):
    #  Run bash command line and capture the output
    process = subprocess.Popen(command,stdout=subprocess.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()
    return(proc_stdout)
def validate(rawstring):
    #  Check CHR strings
    string=re.sub('(?i)' + re.escape('chr'), '', rawstring)
    string=re.sub('(?i)' + re.escape('mt'), 'M', string)
    if string=="X" or string=="Y" or string=="M":
        return string
    else:
        assert int(string) in range(1,23),'{} is not a valid chromosome'.format(rawstring)
        return string
def checksequence(sequence):
    #  Check DNA sequence make sense
    valid_nuc=set('atcgATCGN')
    seqset=set(sequence)
    if seqset-valid_nuc!=set():
        raise ValueError("Invalid DNA sequence")
def masksequence(sequence):
    mask=sequence.translate(str.maketrans('atcg','NNNN'))
    checksequence(mask)
    return mask
def getsequence(CHR,Start,Length):
    # Get sequence from reference
    fasta_file=genomeDir+genomeBuild+"/SNP_masked_chr"+str(CHR)+".fa"
    f = open(fasta_file, 'r')
    f.seek(Start-1)
    sequence=f.read(Length)
    checksequence(sequence)
    f.close()
    return sequence
def rc(sequence):
    c=sequence.translate(str.maketrans('atcgATCG','tagcTAGC'))
    r_c=c[::-1]
    checksequence(r_c)
    return r_c
def overlap(primer1,primer2):
    # return True if overlap, False if not overlap
    # use strict criteria
    answer=True
    left1start,left1end=primer1.cood[0],primer1.cood[1]
    left2start,left2end=primer2.cood[0],primer2.cood[1]
    right1start,right1end=primer1.cood[2],primer1.cood[3]
    right2start,right2end=primer2.cood[2],primer2.cood[3]
    if left2start<left1start and left2start>left1end and left2end<left1start and left2end>left1end\
        and right2start<right1start and right2start>right1end and right2end<right1start and right2end>right1end:
        answer=False
    return answer
##########################
#####3 Dynamic P3 Settings
#####################3
def getP3setting(conditions,settingsfile=defaultfile):
    #  Get a P3 setting file using conditions
    #  Warning, this does not result in a valid p3settings file
    #  Additional missing information will be added in bpp3obj.getp3input()
    [maskstatus,procedure,optT,deltaT,size,minsize,optsize]=conditions
    mask=int(maskstatus=="Mask")
    minT=optT-deltaT
    maxT=optT+deltaT
    cmd='sed "s/optT/{}/g" {} |sed "s/minT/{}/g"|sed "s/maxT/{}/g"|sed "s/mask/{}/g"|sed "s/minS/{}/g" > P3tempSettings.temp'\
                              .format(optT,settingsfile,minT,maxT,mask,minsize)
    returncode=subprocess.call(cmd,shell=True)
    return returncode
 ########################################   
#### Primer 3 and Primer objects ########3
##################################3
class p3output():
    #  This contains & parse teh P3 output into a easy to use python object
    # @profile
    def __init__(self,rawstring):
        #  Parse the raw P3 output into a data structure, with relevant information organized into self.cleandict
        decodep3output=rawstring.decode("ascii")
        teststring=decodep3output[:-1]
        testlist=teststring.split("\n")
        cleanlist=['\''+re.sub('(?i)' + re.escape('='), '\':\'', s)+'\'' for s in testlist if s.count("=")==1]
        cleanstring=",".join(cleanlist)
        self.cleandict=ast.literal_eval('{'+cleanstring+'}')
        try:
            # Count the number of primers returned. 0 primer will return self.npair=0
            self.npair=int(self.cleandict['PRIMER_PAIR_NUM_RETURNED'])
        except:
            #  Catch when primer3 returns error, prints the error message and escaping with KeyError
            print(self.cleandict['PRIMER_ERROR'])
            raise KeyError("P3RunError")
    def __getitem__(self,index):
        #  self[i] will return the information for the ith primer 
        if not index<self.npair and not index>-1:
            raise IndexError('{} is not a valid index'.format(index))
        return {k:self.cleandict[k] for k in self.cleandict.keys() \
                if '_'+str(index)+'_' in k or k.endswith('_'+str(index))}
class primerobj():
    #  This object contains all relevant information for a single primer
    def __init__(self,rawp3input,coord,length,condition,ID,chrA,chrB):
        #  Initialize primer with information from p3output object and other metadata in self.data
        self.leftprimer=[rawp3input[k] for k in rawp3input.keys() if re.match("PRIMER_LEFT_[0-9]*_SEQUENCE",k)!=None]
        ### Paranoid checking, checks whether the somehow there is more than 1 sequence for left/right primer
        if len(self.leftprimer)>1:
            raise ValueError("Two or more left sequence")
        self.leftprimer=self.leftprimer[0]
        self.rightprimer=[rawp3input[k] for k in rawp3input.keys() if re.match("PRIMER_RIGHT_[0-9]*_SEQUENCE",k)!=None]
        if len(self.rightprimer)>1:
            raise ValueError("Two or more right sequence")
        self.rightprimer=self.rightprimer[0]
        #  Parses the left & right coordinates from the P3output
        self.left=mapint([rawp3input[k] for k in rawp3input.keys() if re.match("PRIMER_LEFT_[0-9]*$",k)!=None][0])
        self.right=mapint([rawp3input[k] for k in rawp3input.keys() if re.match("PRIMER_RIGHT_[0-9]*$",k)!=None][0])
        #  Returns a list of position [left_start,left_end,right_start,right_end]
        self.cood=[coord[0]+(self.left[0]-1),coord[0]+(self.left[0]-1)+(self.left[1]-1),coord[2]\
                   +(self.right[0]-length-1)-(self.right[1]-1),coord[2]+(self.right[0]-length-1)]
        self.prod_size=int([rawp3input[k] for k in rawp3input.keys() if re.match("PRIMER_PAIR_[0-9]*_PRODUCT_SIZE$",k)!=None][0])
        self.ltm=float([rawp3input[k] for k in rawp3input.keys() if re.match("PRIMER_LEFT_[0-9]*_TM$",k)!=None][0])  # Left Tm
        self.rtm=float([rawp3input[k] for k in rawp3input.keys() if re.match("PRIMER_RIGHT_[0-9]*_TM$",k)!=None][0])  # Right TM
        self.condition=condition
        self.ID=ID
        self.chrA=chrA
        self.chrB=chrB
        self.location=[self.chrA+":"+str(self.cood[0]),self.chrA+":"+str(self.cood[1]),self.chrB+":"+str(self.cood[2]),self.chrB+":"+str(self.cood[3])]
        #  Initialize a empty dictionary for QC information. 
        self.QC={}
        self.score=Score(self.QC)
    def __repr__(self):
        return '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.\
        format(self.ID,self.condition,self.leftprimer,self.rightprimer,self.prod_size,self.location,self.ltm,self.rtm,self.QC,self.score)
    # @profile
    def do_BLAT(self,port):
        #  Do one & two sided blat using the BLAT server setup in bpp3obj, parse the information returned
        #  Return True/False value for primerobj.QC['blat/oneblat']
        self.QC['oneblat']=False
        self.QC['blat']=False
        with open ("temp.fa",'w') as f:
            f.write('>{}:LEFT\n'.format(self.ID))
            f.write('{}\n'.format(self.leftprimer))
            f.write('>{}:RIGHT\n'.format(self.ID))
            f.write('{}\n'.format(self.rightprimer))
        subprocess.call("gfClient localhost {} '' temp.fa out.psl -minScore=20 -nohead".format(str(port)),shell=True,stdout=subprocess.PIPE)
        f=open("out.psl")
        string=f.read()
        [left,right]=[string.count("LEFT"),string.count("RIGHT")]
        if left == 1 or right ==1:
            self.QC['oneblat']=True
        if left == 1 and right ==1:
            self.QC['blat']=True
        self.score=Score(self.QC)
    # @profile
    def do_isPCR(self,port):
        #  Do isPCR with BLAT server on port
        #  Returns the number of isPCR hits into QC.['isPCR']
        proc=subprocess.Popen("{} localhost {} '' {} {} stdout -out=bed|wc -l".format(gfpcr,port, self.leftprimer,self.rightprimer),shell=True,stdout=subprocess.PIPE)
        output=int(proc.stdout.read())
        self.QC['isPCR']=output
        self.score=Score(self.QC)
        
##############################
class primerobj2():
    #  This object contains all relevant information for a single primer
    def __init__(self,rawp3input,coord,condition,ID,chrA,chrB):
        #  Initialize primer with information from p3output object and other metadata in self.data
        self.leftprimer=[rawp3input[k] for k in rawp3input.keys() if re.match("PRIMER_LEFT_[0-9]*_SEQUENCE",k)!=None]
        ### Paranoid checking, checks whether the somehow there is more than 1 sequence for left/right primer
        if len(self.leftprimer)>1:
            raise ValueError("Two or more left sequence")
        self.leftprimer=self.leftprimer[0]
        self.rightprimer=[rawp3input[k] for k in rawp3input.keys() if re.match("PRIMER_RIGHT_[0-9]*_SEQUENCE",k)!=None]
        if len(self.rightprimer)>1:
            raise ValueError("Two or more right sequence")
        self.rightprimer=self.rightprimer[0]
        #  Parses the left & right coordinates from the P3output
        self.left=mapint([rawp3input[k] for k in rawp3input.keys() if re.match("PRIMER_LEFT_[0-9]*$",k)!=None][0])
        self.right=mapint([rawp3input[k] for k in rawp3input.keys() if re.match("PRIMER_RIGHT_[0-9]*$",k)!=None][0])
        #  Returns a list of position [left_start,left_end,right_start,right_end]
        self.cood=[coord[0]+(self.left[0]-1),coord[0]+(self.left[0]-1)+(self.left[1]-1),
                  coord[0]+(self.right[0]-1)-(self.right[1]-1),coord[0]+(self.right[0]-1)]
        self.prod_size=int([rawp3input[k] for k in rawp3input.keys() if re.match("PRIMER_PAIR_[0-9]*_PRODUCT_SIZE$",k)!=None][0])
        self.ltm=float([rawp3input[k] for k in rawp3input.keys() if re.match("PRIMER_LEFT_[0-9]*_TM$",k)!=None][0])  # Left Tm
        self.rtm=float([rawp3input[k] for k in rawp3input.keys() if re.match("PRIMER_RIGHT_[0-9]*_TM$",k)!=None][0])  # Right TM
        self.condition=condition
        self.ID=ID
        self.chrA=chrA
        self.chrB=chrB
        self.location=[self.chrA+":"+str(self.cood[0]),self.chrA+":"+str(self.cood[1]),self.chrB+":"+str(self.cood[2]),self.chrB+":"+str(self.cood[3])]
        #  Initialize a empty dictionary for QC information. 
        self.QC={}
        self.score=Score(self.QC)
    def __repr__(self):
        return '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.\
        format(self.ID,self.condition,self.leftprimer,self.rightprimer,self.prod_size,self.location,self.ltm,self.rtm,self.QC,self.score)
    # @profile
    def do_BLAT(self,port):
        #  Do one & two sided blat using the BLAT server setup in bpp3obj, parse the information returned
        #  Return True/False value for primerobj.QC['blat/oneblat']
        self.QC['oneblat']=False
        self.QC['blat']=False
        with open ("temp.fa",'w') as f:
            f.write('>{}:LEFT\n'.format(self.ID))
            f.write('{}\n'.format(self.leftprimer))
            f.write('>{}:RIGHT\n'.format(self.ID))
            f.write('{}\n'.format(self.rightprimer))
        subprocess.call("gfClient localhost {} '' temp.fa out.psl -minScore=20 -nohead".format(str(port)),shell=True,stdout=subprocess.PIPE)
        f=open("out.psl")
        string=f.read()
        [left,right]=[string.count("LEFT"),string.count("RIGHT")]
        if left == 1 or right ==1:
            self.QC['oneblat']=True
        if left == 1 and right ==1:
            self.QC['blat']=True
        self.score=Score(self.QC)
    # @profile
    def do_isPCR(self,port):
        #  Do isPCR with BLAT server on port
        #  Returns the number of isPCR hits into QC.['isPCR']
        proc=subprocess.Popen("{} localhost {} '' {} {} stdout -out=bed|wc -l".format(gfpcr,port, self.leftprimer,self.rightprimer),shell=True,stdout=subprocess.PIPE)
        output=int(proc.stdout.read())
        self.QC['isPCR']=output
        self.score=Score(self.QC)
###############################################
######### Quality Control, BLAT & isPCR  ##########
################################################
def Score(QC):
    #  Scoring function for Qc
    #  Implemented with oneblat, blat & isPCR
    #  Does not care if ther is no QC information, or missing one QC criteria
    #  The score is defined as follows
    #  If blat succeeds, +1
    #  If isPCR hits = 0, +6
    #  If isPCR hit = 1, +3
    #  If oneblat succeeds, +1
    score=0
    for criteria in QC.keys():
        if criteria=="blat" and QC[criteria]==True:
            score+=1
        if criteria=="oneblat" and QC[criteria]==True:
            score+=1
        if criteria=="isPCR":
            if QC[criteria]==0:
                score+=6
            elif QC[criteria]==1:
                score+=3
    return score


        