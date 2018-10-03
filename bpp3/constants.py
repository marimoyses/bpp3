import os
primer3Location = 'primer3_core'; # path for primer3_core command
primer3ParametersLocation = '/opt/bpp3/Software/primer3-2.3.7/src/primer3_config/'; # path for 'primer3_config' file
blatReference = '/hg19/hg19.lex.fa.2bit'; # path for BLAT genome reference files
genomeDir = '/'; # path for directory containing genome build directories (e.g. 'hg19')
genomeBuild = 'hg19'; # genome build (hg19 is default, but user can choose other builds if directories for them exist)
gfpcr='gfPcr'
req_keys=["CHR_A","START_A","DIR_A","CHR_B","START_B","DIR_B","CNVID","SV_Type"]
sv_types=["","DEL","SNV","DUP","INV","INS","CTX_qq","INS:ME:ALU","INS:ME:SVA","INS:ME:LINE1"]
defaultfile = os.path.join(os.path.dirname(__file__), "PD_P3Settings_Run1.txt")

# Need to modify this, so that can add variables easier.
# deltaT=[1,2]  # Tm wiggle room. If optT=60 & deltaT=1 then minT=59 maxT=60
# optT=[60]
# size=[120,200]
# minsize=[80]
# optsize=[100]
# procedure=["paired"]
# maskstatus=["Mask"]
 # P3ranklist will be a list of P3 constraints, from the most strict (Masked, 60,1,150) to the most lenient (Unmask, 60,3,300)
 # The conditions from most important to least is maskstatus, optT(constant as of now), deltaT,size
# P3ranktuple=list(itertools.product(maskstatus,procedure,optT,deltaT,size,minsize,optsize))
# P3ranklist=[list(x) for x in P3ranktuple]
# with open("conditionList.txt","w") as f:
    # for condition in P3ranklist:
        # f.write('{}\n'.format(str(condition)))