import sys
import re
import os.path
import shutil
import glob
import argparse
import subprocess
import json
from time import time, localtime, strftime                        
from Bio import codonalign

Codon_trees_lib_path = os.path.dirname(sys.argv[0]).replace("scripts", "lib")
sys.path.append(Codon_trees_lib_path)
import patric_api
import phylocode

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("genomeIdsFile", type=str, help="file with PATRIC genome IDs, one per line, optional content after tab delimiter ignored")
parser.add_argument("--genomeObjectFile", metavar="file", type=str, help="genome object (json file) to be added to ingroup")
#parser.add_argument("--genomeObjectName", metavar="name", type=str, help="name for genome object")
parser.add_argument("--outgroupIdsFile", metavar="file", type=str, help="ougroup genome ids, one per line (or first column of TSV)")
parser.add_argument("--maxGenes", metavar="#", type=int, default=50, help="number of genes in concatenated alignment")
parser.add_argument("--bootstrapReps", metavar="#", type=int, default=0, help="number of raxml 'fast boostrap' replicates")
parser.add_argument("--maxGenomesMissing", metavar="#", type=int, default=0, help="ingroup genomes allowed to lack a member of any homolog group")
parser.add_argument("--maxAllowedDups", metavar="maxDups", type=int, default=0, help="duplicated gene occurrences allowed within homolog group")
parser.add_argument("--endGapTrimThreshold", metavar="maxPropGaps", type=float, default=0.5, help="stringency of end-gap trimming, lower for less trimming")
parser.add_argument("--raxmlExecutable", metavar="program_name", type=str, default="raxml", help="program to call, possibly with path")
parser.add_argument("--rateModel", metavar="rateModel", type=str, choices = ['CAT', 'GAMMA'], default="CAT", help="variable rate category model CAT|GAMMA")
parser.add_argument("--proteinModel", metavar="substModel", type=str, default="WAGF", help="raxml protein substitution model")
parser.add_argument("--analyzeCodons", action='store_true', help="analyze only codons (ignore amino acids)")
parser.add_argument("--analyzeProteins", action='store_true', help="analyze only amino acids")
parser.add_argument("--threads", metavar="T", type=int, default=2, help="number of threads for raxml")
#parser.add_argument("--runRaxml", action='store_true', help="Deprecated: raxml run by default, use 'deferRaxml' to turn off")
parser.add_argument("--deferRaxml", action='store_true', help="set this flag if you do not want raxml to be run automatically (you can run it manually later using the command file provided)")
parser.add_argument("--outputDirectory", type=str, metavar="out_dir", help="directory for output, create if it does not exist")
parser.add_argument("--pathToFigtreeJar", type=str, metavar="jar_file", help="specify this to generate PDF graphic: java -jar pathToFigtreeJar -graphic PDF CodonTree.nex CodonTree.pdf")
parser.add_argument("--focusGenome", metavar="genome_id", type=str, help="genome to be highlighted in color in Figtree")
parser.add_argument("--debugMode", action='store_true', help="turns on more progress output to log file")
#parser.add_argument("--enableGenomeGenePgfamFileReuse", action='store_true', help="read genes and pgfams from stored file matching genomeIdsFile if it exists")
args = parser.parse_args()
starttime = time()
logfileName = os.path.basename(sys.argv[0])
logfileName = re.sub("\..*", "", logfileName)
logfileName += ".log"
global LOG 
LOG = open(logfileName, 'w')
LOG.write("starting %s\n"%sys.argv[0])
LOG.write(strftime("%a, %d %b %Y %H:%M:%S", localtime(starttime))+"\n")
LOG.write("args= "+str(args)+"\n\n")
LOG.flush()
phylocode.LOG = LOG
patric_api.LOG = LOG

subprocess.check_call(['which', args.raxmlExecutable])
phylocode.checkMuscle()

genomeIds = []
with open(args.genomeIdsFile) as F:
    for line in F:
        m = re.match(r"(\d+\.\d+)", line)
        if m:
            genomeIds.append(m.group(1))
LOG.write("elapsed seconds = %f\n"%(time()-starttime))
LOG.write("from %s got %d genomeIds\n%s\n"%(args.genomeIdsFile, len(genomeIds), "\t".join(genomeIds)))
if args.genomeObjectFile:
    LOG.write("genome object file: %s\n"%args.genomeObjectFile)
LOG.flush()

outgroupIds = []
if args.outgroupIdsFile:
    with open(args.outgroupIdsFile) as F:
        for line in F:
            m = re.match(r"(\d+\.\d+)", line)
            if m:
                outgroupIds.append(m.group(1))
LOG.write("got % outgroupIds\n%s\n"%(len(outgroupIds), "\t".join(outgroupIds)))
LOG.flush()

if len(genomeIds) + len(outgroupIds) < 4:
    LOG.write("too few genomeIds to build a tree with: %d"%(len(genomeIds)+len(outgroupIds)))
    sys.exit(1)

fileBase = os.path.basename(args.genomeIdsFile)
fileBase = re.sub("\..*", "", fileBase)

# if either codons or proteins is specified, analyze just that, otherwise analyze both
if (args.analyzeCodons or args.analyzeProteins):
    if args.analyzeCodons:
        LOG.write("analyzing just codon nucleotides, not proteins\n")
    else:
        LOG.write("analyzing just proteins, not nucleotides\n")
else:
    args.analyzeCodons = args.analyzeProteins = True
    LOG.write("analyzing both codons and proteins\n")
LOG.flush()

if not args.outputDirectory:
    args.outputDirectory=fileBase+"_dir/"
if not args.outputDirectory.endswith("/"):
    args.outputDirectory += "/"
if os.path.exists(args.outputDirectory):
    LOG.write("data directory %s exists\n"%args.outputDirectory)
else:
    os.mkdir(args.outputDirectory)

patric_api.debug = args.debugMode
phylocode.debug = args.debugMode

# this is where we gather the list of Pgfam genes for each ingroup genome ID
genomeGenePgfamList=[]
if False and args.genomeGenePgfamsFile: # reserve for future use (convenient for debugging)
    with open(args.genomeGenePgfamsFile) as F:
        for line in F:
            if "genome" in line:
                continue # header
            row = line.rstrip("\n").split("\t")
            if len(row) == 3:
                genomeGenePgfamList.append(row)
else:
    genomeGenePgfamList = patric_api.getPatricGenesPgfamsForGenomeList(genomeIds)

genomeObject=None
genomeObject_genomeId=None
genomeObject_name=None
if args.genomeObjectFile:
    #try:
    genomeObject = json.load(open(args.genomeObjectFile))
    genomeObjectGenePgfams = patric_api.getPatricGenesPgfamsForGenomeObject(genomeObject)
    genomeGenePgfamList.extend(genomeObjectGenePgfams)
    genomeObject_genomeId = genomeObject['id']
    genomeObject_name = genomeObject['scientific_name']
    genomeIds.append(genomeObject_genomeId)
    args.focusGenome = genomeObject_genomeId
    LOG.write("parsed json file %s, got PGFam genes=%d, total now is %d\n"%(args.genomeObjectFile, len(genomeObjectGenePgfams), len(genomeGenePgfamList)))
    LOG.flush()
    #except Exception as e:
    #LOG.write("Problem reading genome object json file.\n%s\n"%str(e))

# add outgroup genes+pgfams to list, get dynamically as the outgroup might change from run to run
if len(outgroupIds):
    genomeGenePgfamList.extend(patric_api.getPatricGenesPgfamsForGenomeList(outgroupIds))

with open(args.outputDirectory+fileBase+".genomeGenePgfams.txt", 'w') as F:
    for row in genomeGenePgfamList:
        F.write("\t".join(row)+"\n")

LOG.write("got genes and pgfams for genomes, len=%d\n"%len(genomeGenePgfamList))
for row in genomeGenePgfamList[:1]:
    LOG.write("\t".join(row)+"\n")
LOG.flush()
if not len(genomeGenePgfamList):
    LOG.write("got no genes and pgfams for genomes, exiting\n")
    sys.exit(1)

LOG.write("allowing %d genomes missing per PGfam of ingroup (out of %d total)\n"%(args.maxGenomesMissing, len(genomeIds)))
if args.maxGenomesMissing >= len(genomeIds):
    raise Exception("getSingleCopyPgfams: maxGenomesMissing too large: %d"%args.maxGenomesMissing)

# call to getSingleCopyPgfams uses ingroup taxa, outgroup is not involved in selecting single copy pgfams
singleCopyPgfams = phylocode.selectSingleCopyPgfams(genomeGenePgfamList, genomeIds, requiredGenome=args.focusGenome, maxGenomesMissing=args.maxGenomesMissing, maxAllowedDups=args.maxAllowedDups)

LOG.write("got single copy pgfams, num=%d\n"%len(singleCopyPgfams))
if len(singleCopyPgfams) > args.maxGenes:
    singleCopyPgfams=singleCopyPgfams[0:args.maxGenes]
    LOG.write("\tselecting top single-family genes: %d\n"%len(singleCopyPgfams))
LOG.flush()
if not len(singleCopyPgfams):
    LOG.write("got no single copy pgfams, exiting\n")
    sys.exit(1)
with open(args.outputDirectory+fileBase+".singlishCopyPgfams.txt", 'w') as F:
    for pgfam in singleCopyPgfams:
        F.write(pgfam+"\n")

allGenomeIds = genomeIds
allGenomeIds.extend(outgroupIds)
#genesForPgfams = phylocode.getGenesForPgfams(genomeGenePgfamList, allGenomeIds, singleCopyPgfams)
genesForPgfams={}
for pgfam in singleCopyPgfams:
    genesForPgfams[pgfam] = []
genomeObjectGeneDna={}
genomeObjectGenes=set()
numGenesAdded=0
for row in genomeGenePgfamList:
    genome, gene, pgfam = row
    if genome in allGenomeIds and pgfam in genesForPgfams:
        genesForPgfams[pgfam].append(gene)
        if genome == genomeObject_genomeId:
            genomeObjectGenes.add(gene)
        numGenesAdded += 1

LOG.write("got %d genes for %d pgfams\n"%(numGenesAdded, len(genesForPgfams)))
LOG.flush()
for pgfamId in singleCopyPgfams: #genesForPgfams:
    if pgfamId not in genesForPgfams:
        LOG.write("singleCopy pgfamId %s not in genesForPgfams\n"%pgfamId)
        continue
if genomeObject:
    genomeObjectProteins = patric_api.getGenomeObjectProteins(genomeObjectGenes, genomeObject)
    genomeObjectGeneDna = patric_api.getGenomeObjectGeneDna(genomeObjectGenes, genomeObject)

proteinAlignments = {}
codonAlignments = {}
alignedTaxa=set()
#phylocode.generateAlignmentsForCodonsAndProteins(genesForPgfams, proteinAlignments, codonAlignments)
for pgfamId in genesForPgfams: #genesForPgfams:
    proteinSeqRecords = patric_api.getProteinBioSeqRecordsForPatricIds(genesForPgfams[pgfamId])
    if args.genomeObjectFile:
        for geneId in genesForPgfams[pgfamId]:
            if geneId in genomeObjectProteins:
                proteinSeqRecords.append(genomeObjectProteins[geneId])
    proteinAlignment = phylocode.alignSeqRecordsMuscle(proteinSeqRecords)
    proteinAlignment = phylocode.resolveDuplicatesPerPatricGenome(proteinAlignment)
    proteinAlignment.sort()
    try:
        codonAlignment = phylocode.proteinToCodonAlignment(proteinAlignment, genomeObjectGeneDna)
        if codonAlignment: # if an error happened, we don't do next steps
            phylocode.relabelSequencesByGenomeId(codonAlignment)
            if codonAlignment.get_alignment_length() % 3:
                raise Exception("codon alignment length not multiple of 3 for %s\n"%pgfamId)
            if args.endGapTrimThreshold:
                codonAlignment = phylocode.trimEndGaps(codonAlignment, args.endGapTrimThreshold)
            codonAlignments[pgfamId] = codonAlignment
    except Exception as e:
        LOG.write("Exeption aligning codons: %s\n"%str(e))
    phylocode.relabelSequencesByGenomeId(proteinAlignment)
    for seqRecord in proteinAlignment:
        alignedTaxa.add(seqRecord.id)
    if args.endGapTrimThreshold:
        proteinAlignment = phylocode.trimEndGaps(proteinAlignment, args.endGapTrimThreshold)
    proteinAlignments[pgfamId] = proteinAlignment
numTaxa=len(alignedTaxa)

LOG.write("protein and codon alignments completed. num prot als = %d, num codon als = %d\n"%(len(proteinAlignments), len(codonAlignments)))
LOG.write("First prot alignment has %d elements\n"%len(proteinAlignments.values()[0]))
LOG.write("original_id of first prot: %s\n"%proteinAlignments.values()[0][0].annotations['original_id'])
LOG.flush()

# generate hopefully uniq output file name base
phyloFileBase = fileBase+"_%dtaxa"%(numTaxa)
if args.analyzeCodons:
    phyloFileBase += "_%dcds"%len(codonAlignments)
if args.analyzeProteins:
    phyloFileBase += "_%dproteins"%len(proteinAlignments)

proteinPositions=0
codonPositions = 0
# write the genes included in each homology group (and those paralogs excluded)
with open(args.outputDirectory+phyloFileBase+".pgfamsAndGenesIncludedInAlignment.txt", 'w') as F:
    for pgfamId in singleCopyPgfams:
        if pgfamId in genesForPgfams:
            if pgfamId in proteinAlignments:
                proteinPositions += proteinAlignments[pgfamId].get_alignment_length()
                genesNotIncluded = set(genesForPgfams[pgfamId])
                genesIncluded = set()
                F.write(pgfamId+"\tProteins\t")
                for seqRecord in proteinAlignments[pgfamId]:
                    originalId = seqRecord.annotations['original_id']
                    F.write("\t"+originalId)
                    if originalId not in genesNotIncluded:
                        LOG.write("Problem: originalId %s not in genesForPgfams for %s\n"%(originalId, pgfamId))
                    else:
                        genesNotIncluded.remove(originalId)
                    genesIncluded.add(originalId)
                if len(genesNotIncluded):
                    F.write("\tdeletedParalogs: \t"+"\t".join(genesNotIncluded))
                F.write("\n")
                if pgfamId in codonAlignments:
                    F.write(pgfamId+"\tCodons\t")
                    codonPositions += codonAlignments[pgfamId].get_alignment_length()
                    for seqRecord in codonAlignments[pgfamId]:
                        originalId = seqRecord.annotations['original_id']
                        F.write("\t"+originalId)
                        genesIncluded.remove(originalId)
                    if len(genesIncluded):
                        F.write("\tlackingCodonAlignment: \t"+"\t".join(genesIncluded))
                    F.write("\n")
                else:
                    F.write(pgfamId+"\nNo codon alignment\n")
            else:
                F.write(pgfamId+"\tNo protein alignment\n")
        else:
            F.write(pgfamId+" has no genes selected")
        F.write("\n")
    F.close()

# finally, output concatenated protein and/or DNA alignment and partitions and raxml command to appropriate files
raxmlCommand=''
if args.analyzeProteins and args.analyzeCodons:
    phyloFileBase += "_codonAndProteins"
    phylocode.outputCodonsProteinsPhylip(codonAlignments, proteinAlignments, args.outputDirectory+phyloFileBase+".phy")
    with open(args.outputDirectory+phyloFileBase+".partitions", 'w') as PartitionFile:
        for i in range(1,4):
            PartitionFile.write("DNA, codon%d = %d-%d\\3\n"%(i, i, codonPositions))
        PartitionFile.write("%s, proteins = %d-%d\n"%(args.proteinModel, codonPositions+1, codonPositions+proteinPositions))
    raxmlCommand = [args.raxmlExecutable, "-s", phyloFileBase+".phy", "-n", phyloFileBase, "-m",  "GTR%s"%args.rateModel, "-q",  phyloFileBase+".partitions",  "-p", "12345", "-T", str(args.threads)]

elif args.analyzeCodons:
    phyloFileBase += "_codonAlignment"
    phylocode.writeConcatenatedAlignmentsPhylip(codonAlignments, args.outputDirectory+phyloFileBase+".phy")
    with open(args.outputDirectory+phyloFileBase+".partitions", 'w') as PartitionFile:
        for i in range(1,4):
            PartitionFile.write("DNA, codon%d = %d-%d\\3\n"%(i, i, codonPositions))
    raxmlCommand = [args.raxmlExecutable, "-s", phyloFileBase+".phy", "-n", phyloFileBase, "-m",  "GTR%s"%args.rateModel, "-q",  phyloFileBase+".partitions",  "-p", "12345", "-T", str(args.threads)]

elif args.analyzeProteins:
    phyloFileBase += "_proteinAlignment"
    phylocode.writeConcatenatedAlignmentsPhylip(proteinAlignments, args.outputDirectory+phyloFileBase+".phy")
    raxmlCommand = [args.raxmlExecutable, "-s", phyloFileBase+".phy", "-n", phyloFileBase, "-m",  "PROT%s%s"%(args.rateModel, args.proteinModel), "-p", "12345", "-T", str(args.threads)]
raxmlCommand.extend(["-e", "1.0"]) # limit on precision, faster than default 0.1

if args.bootstrapReps > 0:
    raxmlCommand.extend(["-f", "a", "-x", "12345", "-N", str(args.bootstrapReps)]) 
with open(args.outputDirectory+phyloFileBase+".raxmlCommand.sh", 'w') as F:
    F.write(" ".join(raxmlCommand)+"\n")

if not args.deferRaxml:
    #remove RAxML files that clash in name, their existence blocks raxml from running
    for fl in glob.glob(args.outputDirectory+"RAxML_*"+phyloFileBase):
        os.remove(fl)
    proc = subprocess.Popen(raxmlCommand, cwd=args.outputDirectory)
    proc.wait()
    LOG.write("raxml completed: elapsed seconds = %f\n"%(time()-starttime))
    LOG.flush()
    genomeIdToName = {}
    for genomeId, genomeName in patric_api.getNamesForGenomeIds(allGenomeIds):
        genomeIdToName[genomeId] = genomeName+" "+genomeId
    if genomeObject:
        genomeIdToName[genomeObject_genomeId] = genomeObject_name+" "+genomeObject_genomeId
    originalNewick = ""
    raxmlNewickFileName = args.outputDirectory+"RAxML_bestTree."+phyloFileBase
    if args.bootstrapReps > 0:
        raxmlNewickFileName = args.outputDirectory+"RAxML_bipartitions."+phyloFileBase
    F = open(raxmlNewickFileName)
    originalNewick = F.read()
    F.close()
    renamedNewick = phylocode.relabelNewickTree(originalNewick, genomeIdToName)
    F = open(args.outputDirectory+"CodonTree.nwk", 'w')
    F.write(renamedNewick)
    F.close()
    LOG.write("codonTree output newick file saved to CodonTree.nwk\n")
    LOG.flush()

    # test to see if we can write a figtree nexus file
    #
    # We find the template file figtree.nex in the same directory
    # as our library code.
    #
    nexus_template_file = os.path.join(Codon_trees_lib_path, "figtree.nex")
    if os.path.exists(nexus_template_file):
        LOG.write("Found figtree template file: %s\n"%nexus_template_file)
        LOG.flush()
        figtreeParams = phylocode.readFigtreeParameters(nexus_template_file)
        nexusOutfileName = args.outputDirectory+phyloFileBase+".figtree.nex"
        nexusOut = open(nexusOutfileName, "w")
        phylocode.writeTranslatedNexusTree(nexusOut, originalNewick, genomeIdToName, figtreeParameters=figtreeParams, highlightGenome=args.focusGenome)
        nexusOut.close()
        LOG.write("nexus file written to %s\n"%nexusOutfileName)
        LOG.flush()
        #shutil.copy2(nexusOutfileName, args.outputDirectory+"CodonTree.nex")
        if not args.pathToFigtreeJar:
            if os.path.exists(os.path.join(Codon_tree_lib_path, "figtree.jar")):
                args.pathToFigtreeJar = os.path.join(Codon_tree_lib_path, "figtree.jar")
        if args.debugMode:
            LOG.write("found figtree.jar at %s\n"%args.pathToFigtreeJar)
        if os.path.exists(args.pathToFigtreeJar):
            figtreePdfName = args.outputDirectory+phyloFileBase+".figtree.pdf"
            LOG.write("run figtree to create tree figure as %s\n"%figtreePdfName)
            figtreeCommand = ['java',  '-jar', args.pathToFigtreeJar, '-graphic', 'PDF', nexusOutfileName, figtreePdfName]
            subprocess.Popen(figtreeCommand)
        
OUT = open(args.outputDirectory+"CodonTree.stats", 'w')
OUT.write("Statistics for CodonTree\n")
OUT.write("Num_genomes\t%s\n"%numTaxa)
OUT.write("Num_protein_alignments\t%s\n"%len(proteinAlignments))
OUT.write("Num_aligned_amino_acids\t%s\n"%proteinPositions)
OUT.write("Num_CDS_alignments\t%s\n"%len(proteinAlignments))
OUT.write("Num_aligned_nucleotides\t%s\n"%codonPositions)
OUT.write("PGFams\t%s\n"%",".join(sorted(proteinAlignments)))
OUT.write("command_line\t%s\n"%" ".join(raxmlCommand))
OUT.close()
LOG.write("output written to directory %s\n"%args.outputDirectory)
LOG.close()
sys.stdout.write("\n") 

