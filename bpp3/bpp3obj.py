import bpp3.helper as helper
import subprocess
import time
import os
import pickle
#################
# This script calls Primer3 repeatedly, at decreasing order of stringency, until a set is found
# BLAT is performed on the resulting primer pairs
# Header required for input file
# The following are required CHR_A    START_A    DIR_A    CHR_B    START_B    DIR_B    CNVID    SV_Type
# These need to be "exactly" as is. The code will keep other columns, but there won't be analysis done with them
# The headers can be in any order. 
#################
class BpP3obj():   
    def __init__(self,input_file,port):
        #  Initialize obj by parsing the input file
        #  Port is for setting up the BLAT server for QC. Server does NOT start here.
        #  The "CNVID" column is especially important and must be unique. The IDs are used to identify CNVs throughout
        #  self.data is a dictionary with CNVID as keys. All info pertaining for a particular CNVID is contained within self.data[CNVID]
        self.req_keys=helper.req_keys
        self.sv_types=helper.sv_types
        self.data={}
        self.port=port
        with open (input_file,'r') as f:
            keys = f.readline().rstrip().split("\t")
            self.keys=keys
            if len(keys) > len(set(keys)):
                raise IndexError("Un-unique keys")
            if set(self.req_keys)-set(keys)!=set():
                raise KeyError("Wrong file header.")
            IDIndex=keys.index("CNVID")
            for line in f:
                lineDat=line.rstrip().split("\t")
                ID=lineDat[IDIndex]
                if ID in self.data.keys():
                    raise ValueError("Duplicate CNV ID")
                self.data[ID]=dict(zip(keys,lineDat))    
                # For SNV and INS the start end are the same. Add 1 to start_B so that the sequence grabbing does not overlap
                if self.data[ID]["START_A"]==self.data[ID]["START_B"] and self.data[ID]["CHR_A"]==self.data[ID]["CHR_B"]: 
                    self.data[ID]["START_B"]=str(int(self.data[ID]["START_B"])+1)
    def ids(self):
        return self.data.keys()
    def sanitate(self):
        # Checks for errors and correct idiosyncracies in the input. Throws exception if problem. 
        # TODO add sanity check for type of SV
        for ID in self.data.keys():
            # Check chromosome fix if neccessary
            self.data[ID]["CHR_A"]=helper.validate(self.data[ID]["CHR_A"])
            self.data[ID]["CHR_B"]=helper.validate(self.data[ID]["CHR_B"])    
            # Check orientation
            assert self.data[ID]["DIR_A"]=="+" or self.data[ID]["DIR_A"]=="-",\
            '{} is not valid orientation'.format( self.data[ID]["DIR_A"])
            assert self.data[ID]["DIR_B"]=="+" or self.data[ID]["DIR_B"]=="-",\
            '{} is not valid orientation'.format( self.data[ID]["DIR_B"])
            # Check position is number and make it int            
            try:
                self.data[ID]["START_A"]=int(self.data[ID]["START_A"])
                self.data[ID]["START_B"]=int(self.data[ID]["START_B"])
            except ValueError:
                raise ValueError('{} and/or {} is not valid numerical position'.format(self.data[ID]["START_A"],self.data[ID]["START_B"]))
            # check SV type make sense TODO
            assert self.data[ID]["SV_Type"] in self.sv_types,'{} is not a supported SV type'.format(self.data[ID]["SV_Type"])    
    #  For
    # @profile
    def getP3input(self,ID,llength,rlength,buffer,condition=['Mask', 'paired', 60, 3, 140, 50, 140]):
        #  With CNV ID, left/right length, buffer and conditions, obtain a reasonable P3 settings
        #  No DNA information used here. 
        #  Condition is (Mask status, Opt_Tm, Tm_wiggle room, opt_PROD_size)
        #  We assume that opt-prod-size=max-prod-size
        #  Output into self.data[ID]['p3input'], as well as modifying the P3 Settings file (MaxS and OptS)
        #  This function returns a return code. Any value of return code other than 0 indicates an error
        #  This function also handles cases where length>size for DUPs. 
        #  When that occurs, we use size//2 for rlength and llength, and also outputs a warning message
        # print('left:{}right:{}'.format(str(llength),str(rlength)))
        returncode=helper.getP3setting(condition)
        minS,optsize=condition[5],condition[6]
        if returncode!=0:
            raise ValueError("GetP3SettingFailed")
        p3length=llength+rlength
        # first implement as before 
        # We handle DUP cases with the -/+  orientations. If it is +/- orientation then assume doing normal direction PCR
        if self.data[ID]['SV_Type']=="DUP" and self.data[ID]['DIR_A']=="+" and self.data[ID]['DIR_B']=="-": 
            print("For SVID {} we treat DUP as DEL with normal direction PCR primers".format(ID))
        if self.data[ID]['SV_Type']=="DUP" and self.data[ID]['DIR_A']=="-" and self.data[ID]['DIR_B']=="+": 
            size = abs(int(self.data[ID]['START_A'])-int(self.data[ID]['START_B']))
            if size<p3length: # Need to change this
                rlength=size//2
                llength=size//2
                minS=max(size-30,55)
                condition_special=list(condition)
                condition_special[4]=minS
                returncode=helper.getP3setting(condition_special)
                if returncode!=0:    
                    raise ValueError("GetP3SettingFailed")
                print("For SVID {} DUP size is too small for specified length {}, using size {} as length instead".format(ID,p3length,size))
        if self.data[ID]["DIR_A"]=="+":
            startA=self.data[ID]["START_A"]-llength+1
            seqA=helper.getsequence(self.data[ID]["CHR_A"],startA,llength)
            endA=startA+llength-1
        else:
            startA=self.data[ID]["START_A"]    
            chrA=self.data[ID]["CHR_A"]
            seqA=helper.rc(helper.getsequence(chrA,startA,llength))
            endA=-startA
            startA=-(startA+llength-1)
        if self.data[ID]["DIR_B"]=="-":
            startB=self.data[ID]["START_B"]
            seqB=helper.getsequence(self.data[ID]["CHR_B"],startB,rlength)
            endB=startB+rlength-1
        else:
            startB=self.data[ID]["START_B"]-rlength+1
            seqB=helper.rc(helper.getsequence(self.data[ID]["CHR_B"],startB,rlength))
            endB=-startB
            startB=-(startB+rlength-1)
        
       
        length=llength+rlength
        self.data[ID]['p3inputCood']=[startA,endA,startB,endB]
        self.data[ID]['p3input']=seqA+seqB
        self.data[ID]['llength']=llength
        cmd='sed -i "s/optS/{}/;s/maxS/{}/g" P3tempSettings.temp'.format(str(optsize),str(length))
        return_code=subprocess.call(cmd,shell=True)
    # @profile
    def runP3(self,buffer,maxlength=500):
        #  This runs the entire P3 analysis, iterativly loosening the P3 constraints as indicated in self.helper.P3ranklist
        #  When P3 returns non-zero result for condition, stop and will not consider lesser conditions
        #  Getting the sequence dictated by the P3 settings, joining together, masking are done here
        #  The raw P3 command is saved in self.data[ID]['p3command']
        #  The P3 output is parsed by "p3outputobj" in helper into a p3outputobj class object. For more info see doc in helper.py       #
        #  Each primer in p3outputobj is parsed into primerobj. For more information see helper.py
        #  The primers are gathered into the list self.data[ID]['primers'] for each CNVID
        #  The CNVs that did not result in a single valid primer is collected in the list self.failed
        self.failed=[]
        for ID in self.ids():
            print("Running P3 on CNV ID {}".format(ID))
            self.data[ID]['primers']=[]            
            condition_number=0
            npairs=0
            while npairs==0 and condition_number<len(helper.P3ranklist):  # Iterate through the conditions until non zero result
                condition=helper.P3ranklist[condition_number]
                [maskstatus,procedure,optT,deltaT,length,minsize,optsize]=condition
                self.data[ID]['maskstatus']=maskstatus
                try:
                    self.data[ID]['maskstatus']
                except:
                    raise ValueError("!!!")
                # Parameters of either paired end or single  end procedures
                if procedure=="paired":
                    llength,rlength=length//2,length//2   # Initial length
                elif procedure=="left":
                    llength,rlength=length,maxlength
                elif procedure=="right":
                    llength,rlength=maxlength,length
                # Catch problems for INV
                if self.data[ID]['SV_Type']=="INV":
                    size = abs(int(self.data[ID]['START_A'])-int(self.data[ID]['START_B']))
                    if self.data[ID]['DIR_A']=="+" and self.data[ID]['DIR_B']=="+":
                        if rlength>size:
                            rlength=size-1
                            print('For CNV {} the right length is longer than the SV length, using rlength={} instead'.format(ID,str(rlength)))
                    elif self.data[ID]['DIR_A']=="-" and self.data[ID]['DIR_B']=="-":
                        if llength>size:
                            llength=size-1
                            print('For CNV {} the left length is longer than the SV length, using llength={} instead'.format(ID,str(llength)))
                    else:
                        raise ValueError("SV ID {} is marked as INV but the breakpoint orientation is incorrect.".format(ID))                
                self.getP3input(ID,llength,rlength,buffer,condition)   # In case of DUPs rlength/llength will be changed
                setting="P3tempSettings.temp"  
                maskstatus=condition[0]
                if maskstatus=="Mask":
                    DNA=helper.masksequence(self.data[ID]['p3input'])
                elif maskstatus=="Unmask":
                    DNA=self.data[ID]['p3input']
                elif maskstatus!="Unmask":
                    raise ValueError('{} is not valid masking status'.format(maskstatus))
                inputData='\'SEQUENCE_TEMPLATE='+DNA+'\nSEQUENCE_TARGET='+str(self.data[ID]['llength']-buffer)+','+str(2*buffer)+'\nPRIMER_THERMODYNAMIC_PARAMETERS_PATH='+helper.primer3ParametersLocation+'\n=\''
                p3Command='echo -e '+inputData+' | '+helper.primer3Location+' -p3_settings_file='+setting
                self.data[ID]['p3command']=p3Command
                # Primer 3 is run here. 
                rawp3output=helper.subprocess_cmd(p3Command)           
                p3outputobj=helper.p3output(rawp3output)
                self.data[ID]['p3command']=p3Command
                self.data[ID]['p3output']=p3outputobj
                npairs=p3outputobj.npair
                if npairs !=0:
                    npairs=0
                    for i in range(0,p3outputobj.npair):
                        primer=helper.primerobj(p3outputobj[i],self.data[ID]['p3inputCood'],llength,condition,ID,self.data[ID]['CHR_A'],self.data[ID]['CHR_B'])
                        primer.do_isPCR(self.port)
                        if primer.QC['isPCR'] ==0 or (primer.QC['isPCR']==1):# and self.data['DIR_A']=="+" and self.data[DIR_B]=="-"):
                            self.data[ID]['primers'].append(primer)
                            npairs+=1                        
                condition_number+=1
            if npairs==0:
                print("No primers found for CNV ID {}".format(ID))
                self.failed.append(ID)
            os.unlink("P3tempSettings.temp")
            self.failed.sort()
    def printprimers(self,outputfile):
        #  Print all primers
        #  Will print the CNV informations along with the primer informations
        with open(outputfile,'w') as g:
            header="\t".join(self.keys)+"\tID\tcondition\tleftprimer\trightprimer\tprod_size\tlocation\tltm\trtm\tQC\tScore\n"
            g.write(header)
            for ID in self.rankedIDlist:
                sampledata="\t".join([str(self.data[ID][k]) for k in self.keys])
                if self.data[ID]['primers']==[]:
                    print("No primers found for SVID {}".format(ID))
                else:
                    for primer in self.data[ID]['primers']:
                        g.write(sampledata+"\t"+repr(primer))
    def printHQdiff(self,outputfile):
        #  Print a pair of primers, that are of the highest quality, that are most different from each other
        #  Will print the CNV informations along with the primer informations
        with open(outputfile,'w') as g:
            header="\t".join(self.keys)+"\tID\tcondition\tleftprimer\trightprimer\tprod_size\tlocation\tltm\trtm\tQC\tScore\n"
            g.write(header)
            for ID in self.rankedIDlist:
                sampledata="\t".join([str(self.data[ID][k]) for k in self.keys])
                if self.data[ID]['primers']==[]:
                    print("No primers found for SVID {}".format(ID))
                else:
                    bestprimerlist=[l for l in self.data[ID]['primers'] if l.score==self.data[ID]['primers'][0].score]                  
                    top2list=[bestprimerlist[0],bestprimerlist[-1]]
                    if helper.overlap(top2list[0],top2list[1])==True:
                        bestprimerlist=[l for l in self.data[ID]['primers'] if (self.data[ID]['primers'][0].score-l.score)<2]
                        bestprimerlist.sort(key=lambda x: (x.cood), reverse=True)
                        try: 
                            top2list=[bestprimerlist[0],bestprimerlist[-1]]
                        except:
                            print(bestprimerlist)
                    for primer in top2list:
                        g.write(sampledata+"\t"+repr(primer))
                        
    def printbestprimers(self,number,outputfile):
        #  Print the top 3 primers for each CNV
        #  Will print the CNV informations along with the primer informations
        with open(outputfile,'w') as g:
            header="\t".join(self.keys)+"\tID\tcondition\tleftprimer\trightprimer\tprod_size\tlocation\tltm\trtm\tQC\tScore\n"
            g.write(header)
            for ID in self.rankedIDlist:
                sampledata="\t".join([str(self.data[ID][k]) for k in self.keys])
                if self.data[ID]['primers']==[]:
                    print("No primers found for SVID {}".format(ID))
                else:
                    n=0
                    while n<number:
                        try:
                            primer = self.data[ID]['primers'][n]
                        except IndexError:
                            print('Small Primer Set for CNVID {}'.format(ID))
                            n=number
                        g.write(sampledata+"\t"+repr(primer))
                        n+=1
    def printfailedprimers(self,outputfile):
        with open (outputfile,"w") as g:
            header="\t".join(self.keys)+"\tP3InputCoodinates\tP3InputSequence\n"
            g.write(header)
            for ID in self.failed:
                sampledata="\t".join([str(self.data[ID][k]) for k in self.keys])
                g.write('{}\t{}\t{}\n'.format(sampledata,str(self.data[ID]['p3inputCood']),self.data[ID]['p3input']))
    # @profile          
    def QC(self):
        #  Run BLAT & isPCR for all primers
        #  For more information on BLAT and isPCR see bpp3.helper
        #  This also sorts the primers by quality score
        if self.checkBLATserver()!=0:
            raise ValueError("The BLAT/PCR server has not been set up properlly")
        for ID in self.ids():
            print("Starting QC on SVID {}".format(ID))            
            for primer in self.data[ID]['primers']:
                  primer.do_BLAT(self.port)
                  # primer.do_isPCR(self.port)
            self.data[ID]['primers'].sort(key=lambda x: (x.score,x.cood[0],x.cood[2]), reverse=True)
            self.data[ID]['score']=0
            if self.data[ID]['primers']!=[]:            
                self.data[ID]['score']=self.data[ID]['primers'][0].score+9*(self.data[ID]['maskstatus']=="Mask")
        self.rankedIDlist=list(self.ids())
        self.rankedIDlist.sort(key=lambda x: self.data[x]['score'], reverse=True)
        try:
            os.unlink("out.psl")
            os.unlink("temp.fa")
        except OSError:
            pass
    def startBLATserver(self):
        #  This starts a local background gfServer for BLAT and isPCR on localhost:self.port 
        #  The server can be reused, running this code twice will show "Server already setup"
        #  When starting server, returns server ready status every 5 seconds. 
        return_code=self.checkBLATserver()
        if return_code != 0:
            blatservercmd="gfServer start localhost "+str(self.port)+" "+helper.blatReference+' -tileSize=11 -log=blat.log -stepSize=5 -canStop'
            print("Starting BLAT server on localhost:{}".format(str(self.port)))
            subprocess.Popen(blatservercmd,shell=True,stdout=subprocess.PIPE)
            timepassed=0
            while return_code != 0:
                if return_code==255:
                    print('Blat server not ready, {} seconds, return code {}'.format(timepassed,return_code))
                    return_code=self.checkBLATserver()
                    time.sleep(5)
                    timepassed+=5
                else:
                    raise ValueError('server error return code {}'.format(str(return_code)))
            print("Server Ready")
            return return_code
        else:
            print("Server already setup")
            return return_code
    def stopBLATserver(self):
        blatservercmd="gfServer stop localhost "+str(self.port)
        return_code=subprocess.call(blatservercmd,shell=True)
        return return_code
    def checkBLATserver(self):
        cmd='gfServer status localhost '+str(self.port)
        return_code=subprocess.call(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        return return_code
    def pickle_dump(self,picklefile):
        #  Dump the entire data obj into a file. 
        with open (picklefile,'wb') as f:
            pickle.dump(self,f)

