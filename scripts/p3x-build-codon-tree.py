import sys
import re
from math import sqrt
import os.path
import glob
import argparse
import subprocess
import json
import StringIO
from time import time, localtime, strftime                        
import patric_api
import phylocode
from Bio import SeqIO
from Bio.Alphabet import IUPAC

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--parametersJson", metavar="file.json", type=str, help="parameters in json format (command line overrides)")
parser.add_argument("--outputBase", metavar="filebase", type=str, help="base name for output files, def=codontree")
parser.add_argument("--genomeIdsFile", metavar="file", type=str, nargs="*", help="file with PATRIC genome IDs, one per line (or first column of TSV)")
parser.add_argument("--genomeGroupName", metavar="name", type=str, nargs="*", help="name of user's genome group at PATRIC")
parser.add_argument("--genomeObjectFile", metavar="file", type=str, help="genome object (json file)")
parser.add_argument("--genomePgfamGeneFile", metavar="file", type=str, help="read geneIDs per PGFam per genome from this file")
parser.add_argument("--optionalGenomeIdsFile", metavar="file", type=str, help="optional genome ids, one per line (or first column of TSV)")
parser.add_argument("--maxGenes", metavar="#", type=int, default=50, help="number of genes in concatenated alignment")
parser.add_argument("--excessGenesProp", metavar="prop", type=float, default=0.5, help="multiplier of maxGenes to add to filter out low-scoring alignments")
parser.add_argument("--excessGenesFixed", metavar="#", type=int, default=20, help="fixed excess genes to add to filter out low-scoring alignments")
parser.add_argument("--bootstrapReps", metavar="#", type=int, default=100, help="number of raxml 'fast boostrap' replicates")
parser.add_argument("--maxGenomesMissing", metavar="#", type=int, default=0, help="genomes allowed to lack a member of any homolog group")
parser.add_argument("--maxAllowedDups", metavar="maxDups", type=int, default=0, help="duplicated gene occurrences allowed within homolog group")
parser.add_argument("--endGapTrimThreshold", metavar="maxPropGaps", type=float, default=0.5, help="stringency of end-gap trimming, lower for less trimming")
parser.add_argument("--raxmlExecutable", metavar="program_name", type=str, default="raxml", help="program to call, possibly with path")
parser.add_argument("--rateModel", metavar="model", type=str, choices = ['CAT', 'GAMMA'], default="CAT", help="variable rate category model CAT|GAMMA")
parser.add_argument("--proteinModel", metavar="model", type=str, default="WAGF", help="raxml protein substitution model")
parser.add_argument("--analyzeCodons", action='store_true', help="analyze only codon nucleotides")
parser.add_argument("--analyzeProteins", action='store_true', help="analyze only amino acids")
parser.add_argument("--threads", metavar="T", type=int, default=2, help="threads for raxml")
parser.add_argument("--deferRaxml", action='store_true', help="does not run raxml")
parser.add_argument("--writePgfamAlignments", action='store_true', help="write fasta alignment per pgfam used for tree")
parser.add_argument("--writePgfamMatrix", action='store_true', help="write table of counts per pgfam per genome")
parser.add_argument("--outputDirectory", type=str, metavar="out_dir", help="for output, create if needed")
parser.add_argument("--pathToFigtreeJar", type=str, metavar="path", help="to generate PDF: java -jar path.jar -graphic PDF CodonTree.nex CodonTree.pdf")
parser.add_argument("--universalRolesFile", type=str, metavar="path", help="path to file with universal roles to select conserved genes")
parser.add_argument("--focusGenome", metavar="id", type=str, help="to be highlighted in color in Figtree")
parser.add_argument("--debugMode", action='store_true', help="more output to log file")
parser.add_argument("--authToken", metavar="STRING", type=str, help="patric authentication token")
parser.add_argument("--ignoreAuthEnv", action='store_true', help="turn off authorization by environmental variable")
parser.add_argument("--ignoreAuthRC", action='store_true', help="turn off authorization by file")
#parser.add_argument("--enableGenomeGenePgfamFileReuse", action='store_true', help="read genes and pgfams from stored file matching genomeIdsFile if it exists")
args = parser.parse_args()
starttime = time()
genomeIds = set() # list of genome IDs for tree building (enforcing maxGenomesMissing and maxAllowedDupes)
optionalGenomeIds = set()

if args.parametersJson:
    # read parameters in json file, avoiding overwriting any specified on command line
    params = json.load(open(args.parametersJson))
    if "output_path" in params and not args.outputDirectory:
        args.outputDirectory = os.path.abspath(params["output_path"])
    if "output_file" in params and not args.outputBase:
        args.outputBase = params["output_file"]
    if "genome_group" in params: #not necessary if UI flattens these out to ids
        if not args.genomeGroupName:
            args.genomeGroupName = []
        args.genomeGroupName.extend(params["genome_group"])
    if "genome_ids" in params: 
        genomeIds |= set(params["genome_ids"])
    if "optional_genome_ids" in params:
        optionalGenomeIds |= set(params["optional_genome_ids"])
    if "number_of_genes" in params:
        args.maxGenes = params["number_of_genes"]
    if "bootstraps" in params:
        if params["bootstraps"]:
            args.bootstrapReps = 100
    if "max_genomes_missing" in params:
        args.maxGenomesMissing = params["max_genomes_missing"]
    if "max_allowed_dups" in params:
        args.maxAllowedDups = params["max_allowed_dups"]
    
if not args.outputDirectory:
    args.outputDirectory="./"
if not args.outputDirectory.endswith("/"):
    args.outputDirectory += "/"
if not os.path.exists(args.outputDirectory):
    os.makedirs(args.outputDirectory)
if not args.outputBase:
    args.outputBase = "codontree"

logfileName = os.path.basename(sys.argv[0])
logfileName = re.sub("\..*", "", logfileName)
logfileName += ".log"
logfileName = os.path.join(args.outputDirectory, logfileName)

global LOG 
LOG = open(logfileName, 'w')
LOG.write("starting %s\n"%sys.argv[0])
LOG.write(strftime("%a, %d %b %Y %H:%M:%S", localtime(starttime))+"\n")
LOG.write("args= "+str(args)+"\n\n")
LOG.flush()
phylocode.LOG = LOG
patric_api.LOG = LOG
patric_api.PatricUser = None

if args.authToken:
    if args.debugMode:
        LOG.write("authorizing by passed token: %s\n"%args.authToken)
    patric_api.authenticateByString(args.authToken)
if not patric_api.PatricUser and not args.ignoreAuthEnv:
    if os.environ.has_key("KB_AUTH_TOKEN"):
        if args.debugMode:
            LOG.write("reading auth key from environment\n")
        patric_api.authenticateByString(os.environ.get('KB_AUTH_TOKEN'))
    elif args.debugMode:
        LOG.write("environment variable for authorization not found\n")
if not patric_api.PatricUser and not args.ignoreAuthRC:
    tokenFile = os.path.join(os.environ.get('HOME'), ".patric_token")
    if os.path.exists(tokenFile):
        LOG.write("reading auth key from file %s\n"%tokenFile)
        with open(tokenFile) as F:
            tokenString = F.read().rstrip()
            patric_api.authenticateByString(tokenString)
    elif args.debugMode:
        LOG.write("authorization file %s not found\n"%tokenFile)
LOG.write("Patric User = %s\n"%patric_api.PatricUser)

preflightTests = []
preflightTests.append(("phylocode.which(args.raxmlExecutable) != None", phylocode.which(args.raxmlExecutable) != None))
preflightTests.append(("phylocode.which('muscle') != None", phylocode.which('muscle') != None))

# update genomeIds set
if args.genomeIdsFile:
    for filename in args.genomeIdsFile:
        with open(filename) as F:
            for line in F:
                m = re.match(r"(\d+\.\d+)", line)
                if m:
                    genomeIds.add(m.group(1))
    LOG.write("elapsed seconds = %f\n"%(time()-starttime))
    LOG.write("from %s got %d genome IDs\n%s\n"%(",".join(args.genomeIdsFile), len(genomeIds), "\t".join(genomeIds)))
if args.genomeObjectFile:
    LOG.write("genome object file: %s\n"%args.genomeObjectFile)
LOG.flush()

if args.genomeGroupName:
    if not patric_api.PatricUser:
        raise Exception("No patric user is defined, required for using genome group.\n")
    for groupName in args.genomeGroupName:
        LOG.write("requesting genome IDs for user group %s\n"%groupName)
        for group in args.genomeGroupName:
            groupIds = patric_api.getGenomeGroupIds(groupName)
            genomeIds.update(set(groupIds))
            LOG.write("got %d ids for group %s, total IDs = %d\n"%(len(groupIds), groupName, len(genomeIds)))

if args.optionalGenomeIdsFile:
    with open(args.optionalGenomeIdsFile) as F:
        for line in F:
            m = re.match(r"(\d+\.\d+)", line)
            if m:
                optionalGenomeIds.add(m.group(1))
    LOG.write("got %d optionalGenomeIds\n%s\n"%(len(optionalGenomeIds), "\t".join(optionalGenomeIds)))
    LOG.flush()

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

if args.debugMode:
    patric_api.Debug = True
    phylocode.Debug = True
patric_api.LOG = LOG
phylocode.LOG = LOG


# this is where we gather the list of Pgfam genes for each genome ID
pgfamMatrix = {}
if args.genomePgfamGeneFile: # reserve for future use (convenient for debugging)
    if not os.path.exists(args.genomePgfamGeneFile):
        raise Exception("genomePgfamGeneFile specified as %s does not exist"%args.genomePgfamGeneFile)
    with open(args.genomePgfamGeneFile) as F:
        pgfamMatrix = patric_api.readPgfamGenomeMatrix(F)

genomeObject = None
genomeObject_genomeId = None
genomeObject_name = None
genomeObjectProteins = None
genomeObjectGeneDna = None
if args.genomeObjectFile:
    genomeObject = json.load(open(args.genomeObjectFile))
    genomeObjectProteins = patric_api.getGenomeObjectProteins(genomeObject)
    genomeObjectGeneDna = patric_api.getGenomeObjectGeneDna(genomeObject)
    genomeObjectGenePgfams = patric_api.getPatricGenesPgfamsForGenomeObject(genomeObject)
    for row in genomeObjectGenePgfams:
        genomeId, gene, pgfam = row
        if pgfam not in pgfamMatrix:
            pgfamMatrix[pgfam] = {}
        if genomeId not in pgfamMatrix[pgfam]:
            pgfamMatrix[pgfam][genomeId] = []
        pgfamMatrix[pgfam][genomeId].append(gene)
    
    genomeObject_genomeId = genomeObject['id']
    genomeObject_name = genomeObject['scientific_name']
    genomeIds.add(genomeObject_genomeId)
    args.focusGenome = genomeObject_genomeId
    LOG.write("parsed json file %s"%args.genomeObjectFile)
    LOG.flush()

allGenomeIds = genomeIds
allGenomeIds.update(optionalGenomeIds)
genomesWithData = set()
for pgfam in pgfamMatrix:
    genomesWithData.update(set(pgfamMatrix[pgfam].keys()))
genomesWithoutData = allGenomeIds - genomesWithData
if len(genomesWithoutData):
    if args.universalRolesFile:
        LOG.write("getPgfamGenomeMatrixFromUniversalRoles(%s)"%args.universalRolesFile)
        pgfamMatrix = patric_api.getPgfamMatrixFromUniversalRoles(genomesWithoutData, args.universalRolesFile, pgfamMatrix)
    else:
        LOG.write("getPgfamGenomeMatrix()\n")
        pgfamMatrix = patric_api.getPgfamGenomeMatrix(genomesWithoutData, pgfamMatrix)

if args.writePgfamMatrix:
    with open(os.path.join(args.outputDirectory, args.outputBase+".pgfamMatrix.txt"), 'w') as F:
       patric_api.writePgfamCountMatrix(pgfamMatrix, F)

genesPerGenome = {}
for genome in genomeIds:
    genesPerGenome[genome] = 0
for pgfam in pgfamMatrix:
    for genome in pgfamMatrix[pgfam]:
        genesPerGenome[genome] += 1

deletedGenomes = set()
for genome in genomeIds:
    if genesPerGenome[genome] == 0:
        deletedGenomes.add(genome)
        LOG.write("Deleting genome %s for lack of data.\n"%genome)
genomeIds -= deletedGenomes

LOG.write("got pgfams for genomes, len=%d\n"%len(pgfamMatrix))
LOG.flush()

preflightTests.append(("Need at least 4 genomes to build a tree", len(genomeIds) + len(optionalGenomeIds) >= 4))

if args.maxGenomesMissing:
    LOG.write("allowing %d genomes missing per PGfam (out of %d total)\n"%(args.maxGenomesMissing, len(genomeIds)))
preflightTests.append(("Num genomes - maxGenomesMissing >= 4", len(genomeIds) - args.maxGenomesMissing >= 4))

# call to getSingleCopyPgfams uses main genome, optional genomes are not involved in selecting single copy pgfams
singleCopyPgfams = phylocode.selectSingleCopyPgfams(pgfamMatrix, genomeIds, requiredGenome=args.focusGenome, maxGenomesMissing=args.maxGenomesMissing, maxAllowedDups=args.maxAllowedDups)

LOG.write("got single copy pgfams, num=%d\n"%len(singleCopyPgfams))

singleCopyGenesPerGenome = {}
for genome in genesPerGenome:
    singleCopyGenesPerGenome[genome] = 0
for pgfam in singleCopyPgfams:
    for genome in pgfamMatrix[pgfam]:
        singleCopyGenesPerGenome[genome] += 1

filesToMoveToDetailsFolder =  []
filesToDelete = []

preflightTests.append(("single copy pgfams > 0", len(singleCopyPgfams) > 0))
## perform preflight test
preflightFile = os.path.join(args.outputDirectory, args.outputBase+".preflight")
filesToDelete.append(args.outputBase+".preflight")

maxGenesWithExcess = int(args.maxGenes * (1+args.excessGenesProp) + args.excessGenesFixed)
if len(singleCopyPgfams) > maxGenesWithExcess:
    singleCopyPgfams=singleCopyPgfams[0:maxGenesWithExcess]
    LOG.write("\tevaluating alignments of %d single-family genes, excess of %d\n"%(len(singleCopyPgfams), maxGenesWithExcess - args.maxGenes))

proteinAlignments = {}
proteinAlignmentStats = {}
alignmentScore = {}
alignedTaxa=set()
for pgfamId in singleCopyPgfams:
    geneIdSet = set()
    for genome in pgfamMatrix[pgfamId]:
        geneIdSet.update(set(pgfamMatrix[pgfamId][genome]))
    proteinFasta = patric_api.getSequenceOfFeatures(geneIdSet, 'protein')
    proteinSeqDict = SeqIO.to_dict(SeqIO.parse(StringIO.StringIO(proteinFasta), "fasta", alphabet=IUPAC.extended_protein))
    for genomeId in pgfamMatrix[pgfamId]:
        for geneId in pgfamMatrix[pgfamId][genomeId]:
            if genomeId == genomeObject_genomeId:
                proteinSeqDict[geneId] = genomeObjectProteins[geneId]
        proteinSeqDict[geneId].annotations["genome_id"] = genomeId
        #proteinSeqRecords.append(proteinSeqDict[proteinId])
    if args.debugMode:
        LOG.write("protein set for %s has %d seqs\n"%(pgfamId, len(proteinSeqDict)))
    proteinAlignment = phylocode.alignSeqRecordsMuscle(proteinSeqDict.values())
    proteinAlignment = phylocode.resolveDuplicatesPerPatricGenome(proteinAlignment)
    proteinAlignment.sort()
    proteinAlignments[pgfamId] = proteinAlignment
    alignmentStats = phylocode.calcAlignmentStats(proteinAlignment)
    # dividing by sqrt(alignment length) yields result intermediate between sum_squared_freq and mean_squared_freq (where squared_freq is the within-column sum of squared state (letter, amino acid) frequency
    alignmentScore[pgfamId] = alignmentStats['sum_squared_freq'] / sqrt(alignmentStats['num_pos'])
    proteinAlignmentStats[pgfamId] = alignmentStats

LOG.write("protein alignments completed. num prot als = %d\n"%(len(proteinAlignments)))
# select top alignments by score
singleCopyPgfams = sorted(alignmentScore, key=alignmentScore.get, reverse=True)

if len(singleCopyPgfams) > args.maxGenes:
    singleCopyPgfams=singleCopyPgfams[0:args.maxGenes]
    LOG.write("\tselecting top %d single-family genes based on alignment score\n"%len(singleCopyPgfams))

filteredSingleCopyGenesPerGenome = {}
for genome in genesPerGenome:
    filteredSingleCopyGenesPerGenome[genome] = 0
for pgfam in singleCopyPgfams:
    for genome in pgfamMatrix[pgfam]:
        filteredSingleCopyGenesPerGenome[genome] += 1

genesPerGenomeFile = os.path.join(args.outputDirectory, args.outputBase+".genesPerGenome.txt")
filesToMoveToDetailsFolder.append(args.outputBase+".genesPerGenome.txt")
with open(genesPerGenomeFile, 'w') as F:
    F.write("Genome\tAllGenes\tSingleCopy\tFiltered_SingleCopy\n")
    for genome in sorted(genesPerGenome, key=genesPerGenome.get):
        F.write("%s\t%d\t%d\t%d\n"%(genome, genesPerGenome[genome], singleCopyGenesPerGenome[genome], filteredSingleCopyGenesPerGenome[genome]))

with open(preflightFile, 'w') as F:
    passed = True
    for test in preflightTests:
        passed &= test[1]
        F.write(test[0] + "\t" + str(test[1]) + "\n")
    F.write("All tests passed\t"+str(passed)+"\n")
    #if not passed:
    #    F.write("exiting\n")
    #    F.close()
    #    sys.exit(1)

codonAlignments = {}
for pgfamId in singleCopyPgfams:
    proteinAlignment = proteinAlignments[pgfamId]
    if args.debugMode:
        LOG.write("alignment for %s has %d seqs\n"%(pgfamId, len(proteinAlignment)))
    try:
        codonAlignment = phylocode.proteinToCodonAlignment(proteinAlignment, genomeObjectGeneDna)
        if codonAlignment: # if an error happened, we don't do next steps
            phylocode.relabelSequencesByGenomeId(codonAlignment)
            if codonAlignment.get_alignment_length() % 3:
                raise Exception("codon alignment length not multiple of 3 for %s\n"%pgfamId)
            if args.endGapTrimThreshold:
                codonAlignment = phylocode.trimEndGaps(codonAlignment, args.endGapTrimThreshold)
            codonAlignments[pgfamId] = codonAlignment
            if args.debugMode:
                LOG.write("dna alignment for %s has %d seqs\n"%(pgfamId, len(codonAlignment)))
                #SeqIO.write(codonAlignment[pgfamId][:2], LOG, "fasta")
    except Exception as e:
        LOG.write("Exception aligning codons: %s\n"%str(e))
    phylocode.relabelSequencesByGenomeId(proteinAlignment)
    for seqRecord in proteinAlignment:
        alignedTaxa.add(seqRecord.id)
    # trim protein alignment after codon alignment if trimming enabled
    if args.endGapTrimThreshold:
        proteinAlignment = phylocode.trimEndGaps(proteinAlignment, args.endGapTrimThreshold)
        proteinAlignments[pgfamId] = proteinAlignment 

numTaxa=len(alignedTaxa)

LOG.write("codon alignments completed. num codon als = %d\n"%(len(codonAlignments)))
LOG.flush()

# generate hopefully unique output file name base
phyloFileBase = args.outputBase # + "_%dtaxa"%(numTaxa)

# change to output directory to simplify file naming
os.chdir(args.outputDirectory)

proteinPositions=0
codonPositions = 0

raxmlCommand=''
if len(proteinAlignments):
# write the genes included in each homology group (and those paralogs excluded)
    pgfamsAndGenesIncludedInAlignmentFile = phyloFileBase+".pgfamsAndGenesIncludedInAlignment.txt"
    filesToMoveToDetailsFolder.append(pgfamsAndGenesIncludedInAlignmentFile)
    with open(pgfamsAndGenesIncludedInAlignmentFile, 'w') as F:
        for pgfamId in singleCopyPgfams:
            if pgfamId in proteinAlignments:
                proteinPositions += proteinAlignments[pgfamId].get_alignment_length()
                genesNotIncluded = set()
                for genome in pgfamMatrix[pgfamId]:
                    genesNotIncluded.update(set(pgfamMatrix[pgfamId][genome]))
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
        F.write("\n")
    F.close()

    if args.writePgfamAlignments:
        for pgfam in proteinAlignments:
            SeqIO.write(proteinAlignments[pgfam], pgfam+".faa", "fasta")

    alignmentStatsFile = phyloFileBase+".pgfamAlignmentStats.txt"
    filesToMoveToDetailsFolder.append(alignmentStatsFile)
    with open(alignmentStatsFile, "w") as F:
        first = True
        for pgfam in sorted(alignmentScore, key=alignmentScore.get, reverse=True): #proteinAlignments:
            stats = proteinAlignmentStats[pgfam]
            if first:
                F.write("PGFam\t"+"\t".join(sorted(stats.keys()))+"\tUsedInAnalysis\n")
                first = False
            F.write(pgfam)
            for key in sorted(stats):
                val = stats[key]
                if isinstance(val,int):
                    F.write("\t%d"%stats[key])
                else:
                    F.write("\t%.6f"%stats[key])
            F.write("\t"+str(pgfam in singleCopyPgfams))
            F.write("\n")

# change proteinAlignments to only include selected ones
    selectedAlignments = {}
    for pgfam in singleCopyPgfams:
        selectedAlignments[pgfam] = proteinAlignments[pgfam]
    proteinAlignments = selectedAlignments

# it is possible all codon alignments failed, remove from intent to analyze
    if len(codonAlignments) == 0:
        args.analyzeCodons = False

# finally, output concatenated protein and/or DNA alignment and partitions and raxml command to appropriate files
    if args.analyzeProteins and args.analyzeCodons:
        alignmentFile = phyloFileBase+".phy"
        filesToMoveToDetailsFolder.append(alignmentFile)
        phylocode.outputCodonsProteinsPhylip(codonAlignments, proteinAlignments, alignmentFile)
        with open(phyloFileBase+".partitions", 'w') as PartitionFile:
            for i in range(1,4):
                PartitionFile.write("DNA, codon%d = %d-%d\\3\n"%(i, i, codonPositions))
            PartitionFile.write("%s, proteins = %d-%d\n"%(args.proteinModel, codonPositions+1, codonPositions+proteinPositions))
        raxmlCommand = [args.raxmlExecutable, "-s", phyloFileBase+".phy", "-n", phyloFileBase, "-m",  "GTR%s"%args.rateModel, "-q",  phyloFileBase+".partitions",  "-p", "12345", "-T", str(args.threads)]
        filesToMoveToDetailsFolder.append(phyloFileBase+".partitions")

    elif args.analyzeCodons:
        alignmentFile = phyloFileBase+".phy"
        filesToMoveToDetailsFolder.append(alignmentFile)
        phylocode.writeConcatenatedAlignmentsPhylip(codonAlignments, alignmentFile)
        with open(phyloFileBase+".partitions", 'w') as PartitionFile:
            for i in range(1,4):
                PartitionFile.write("DNA, codon%d = %d-%d\\3\n"%(i, i, codonPositions))
        raxmlCommand = [args.raxmlExecutable, "-s", phyloFileBase+".phy", "-n", phyloFileBase, "-m",  "GTR%s"%args.rateModel, "-q",  phyloFileBase+".partitions",  "-p", "12345", "-T", str(args.threads)]
        filesToMoveToDetailsFolder.append(phyloFileBase+".partitions")

    elif args.analyzeProteins:
        alignmentFile = phyloFileBase+".phy"
        filesToMoveToDetailsFolder.append(alignmentFile)
        phylocode.writeConcatenatedAlignmentsPhylip(proteinAlignments, alignmentFile)
        raxmlCommand = [args.raxmlExecutable, "-s", phyloFileBase+".phy", "-n", phyloFileBase, "-m",  "PROT%s%s"%(args.rateModel, args.proteinModel), "-p", "12345", "-T", str(args.threads)]
    raxmlCommand.extend(["-e", "1.0"]) # limit on precision, faster than default 0.1

    if args.bootstrapReps > 0:
        raxmlCommand.extend(["-f", "a", "-x", "12345", "-N", str(args.bootstrapReps)]) 

    raxmlCommandFile = phyloFileBase+".raxmlCommand.sh"
    with open(raxmlCommandFile, 'w') as F:
        F.write(" ".join(raxmlCommand)+"\n")
    filesToMoveToDetailsFolder.append(phyloFileBase+".raxmlCommand.sh")

genomeIdToName = patric_api.getNamesForGenomeIdsByN(allGenomeIds)
svgTreeImage = None
if len(proteinAlignments) and not args.deferRaxml:
    #remove RAxML files that clash in name, their existence blocks raxml from running
    for fl in glob.glob("RAxML_*"+phyloFileBase):
        os.remove(fl)
    proc = subprocess.Popen(raxmlCommand, stdout=LOG)
    proc.wait()
    LOG.write("raxml completed: elapsed seconds = %f\n"%(time()-starttime))
    LOG.flush()
    if genomeObject:
        genomeIdToName[genomeObject_genomeId] = genomeObject_name+" "+genomeObject_genomeId
    originalNewick = ""
    raxmlNewickFileName = "RAxML_bestTree."+phyloFileBase
    if args.bootstrapReps > 0:
        raxmlNewickFileName = "RAxML_bipartitions."+phyloFileBase
    F = open(raxmlNewickFileName)
    originalNewick = F.read()
    F.close()
    treeWithGenomeIdsFile = phyloFileBase + "_treeWithGenomeIds.nwk"
    #filesToMoveToDetailsFolder.append(treeWithGenomeIdsFile)
    F = open(treeWithGenomeIdsFile, 'w')
    F.write(originalNewick)
    F.close()

    renamedNewick = phylocode.relabelNewickTree(originalNewick, genomeIdToName)
    renamedNewickFile = phyloFileBase+"_treeWithGenomeNames.nwk"
    filesToMoveToDetailsFolder.append(renamedNewickFile)
    F = open(renamedNewickFile, 'w')
    F.write(renamedNewick)
    F.close()
    LOG.write("codonTree newick relabeled with genome names written to "+renamedNewickFile+"\n")
    LOG.flush()

    if args.pathToFigtreeJar and os.path.exists(args.pathToFigtreeJar):
        nexusFilesWritten = phylocode.generateNexusFile(originalNewick, phyloFileBase, nexus_template = None, align_tips = "both", focus_genome = args.focusGenome, genomeIdToName=genomeIdToName)
        LOG.write("nexus file written to %s\n"%(", ".join(nexusFilesWritten)))
        filesToMoveToDetailsFolder.append(nexusFilesWritten[0])
        if len(nexusFilesWritten) > 1:
            filesToDelete.append(nexusFilesWritten[1]) # get rid of codontree_tipsAligned.nex
        for nexusFile in nexusFilesWritten:
            for imageFormat in ("PNG", "SVG", "PDF"):
                imageFile = phylocode.generateFigtreeImage(nexusFile, numTaxa=len(allGenomeIds), figtreeJar=args.pathToFigtreeJar, imageFormat=imageFormat)
                if "_tipsAligned" in nexusFile:
                    filesToMoveToDetailsFolder.append(imageFile)
                if imageFormat == "SVG" and not svgTreeImage:
                    svgTreeImage = imageFile
                LOG.write("created figtree figure: %s\n"%imageFile)

LOG.write(strftime("%a, %d %b %Y %H:%M:%S", localtime(time()))+"\n")
LOG.write("Total job duration %d seconds\n"%(time()-starttime))
        

analysisStatsFile = phyloFileBase+".analysisStats"
OUT = open(analysisStatsFile, 'w')
OUT.write("Statistics for CodonTree Analysis\n")
OUT.write("Num_genomes\t%s\n"%numTaxa)
OUT.write("Num_protein_alignments\t%s\n"%len(proteinAlignments))
OUT.write("Num_aligned_amino_acids\t%s\n"%proteinPositions)
OUT.write("Num_CDS_alignments\t%s\n"%len(proteinAlignments))
OUT.write("Num_aligned_nucleotides\t%s\n"%codonPositions)
OUT.write("PGFams\t%s\n"%",".join(sorted(proteinAlignments)))
OUT.write("command_line\t%s\n"%" ".join(raxmlCommand))
OUT.write("Total job duration %d seconds\n"%(time()-starttime))
OUT.close()
filesToMoveToDetailsFolder.append(analysisStatsFile)


htmlFile = os.path.abspath(args.outputBase+"_report.html")
HTML = open(htmlFile, 'w')
HTML.write("<h1>Phylogenetic Tree Report</h1>\n")
if svgTreeImage is not None and os.path.exists(svgTreeImage):
    HTML.write("<h2>Rendered Tree</h2>\n")
    HTML.write(open(svgTreeImage).read()+"\n\n")
alternateSvgFile = "%s_tipsAligned.svg"%phyloFileBase
if os.path.exists(alternateSvgFile):
    HTML.write("<p><a target='_parent' href='detail_files/%s'>Alternate View</a></p>"%alternateSvgFile)

HTML.write("<p><a target='_parent' href='detail_files/'>Files with more details</a></p>")


HTML.write("<h2>Tree Analysis Statistics</h2>\n<table border='1'>\n")
HTML.write("<tr><td><b>%s</b></td><td>%d</td></tr>\n"%("Num genomes", len(genomeIds)))
if len(deletedGenomes):
    HTML.write("<tr><td><b>%s</b></td><td>%d</td></tr>\n"%("Genomes lacking data", len(deletedGenomes)))
HTML.write("<tr><td><b>%s</b></td><td>%d</td></tr>\n"%("Max Allowed Deletions", args.maxGenomesMissing))
HTML.write("<tr><td><b>%s</b></td><td>%d</td></tr>\n"%("Max Allowed Duplications", args.maxAllowedDups))
HTML.write("<tr><td><b>%s</b></td><td>%d</td></tr>\n"%("Single-copy genes requested", args.maxGenes))
HTML.write("<tr><td><b>%s</b></td><td>%d</td></tr>\n"%("Single-copy genes found", len(singleCopyPgfams)))
HTML.write("<tr><td><b>%s</b></td><td>%d</td></tr>\n"%("Num protein alignments", len(proteinAlignments)))
HTML.write("<tr><td><b>%s</b></td><td>%d</td></tr>\n"%("Num aligned amino acids", proteinPositions))
HTML.write("<tr><td><b>%s</b></td><td>%d</td></tr>\n"%("Num CDS alignments", len(codonAlignments)))
HTML.write("<tr><td><b>%s</b></td><td>%d</td></tr>\n"%("Num aligned nucleotides", codonPositions))
m = re.search("/awe/work/../../([^_]+)_/", os.getcwd())
if m:
    patricJobId = m.group(1)
    HTML.write("<tr><td><b>PATRIC Job Idx=</b></td><td>"+patricJobId+"</td></tr>\n")
# tree analysis may or may not have completed, report only if appropriate
raxmlInfoFile = "RAxML_info."+phyloFileBase
if os.path.exists(raxmlInfoFile):
    filesToMoveToDetailsFolder.append(raxmlInfoFile)
    try:
        raxmlInfo = open(raxmlInfoFile).read()
        raxmlVersion = re.search("RAxML version ([\d\.]+)", raxmlInfo).group(1)
        HTML.write("<tr><td><b>%s</b></td><td>%s</td></tr>\n"%("RAxML Version", raxmlVersion))
        raxmlWarnings = []
        for m in re.finditer("IMPORTANT WARNING: (Sequences.*?identical)", raxmlInfo):
            raxmlWarnings.append(m.group(1))
        if len(raxmlWarnings):
            HTML.write("<tr><td><b>%s</b></td><td>%s</td></tr>\n"%("RAxML Warnings", "<br>\n".join(raxmlWarnings)))
        raxmlDuration = re.search("Overall execution time for full ML analysis: (.*) secs", raxmlInfo).group(1) 
        raxmlDuration = float(raxmlDuration)
        HTML.write("<tr><td><b>%s</b></td><td>%.1f seconds</td></tr>\n"%("RAxML Duration", raxmlDuration))
    except Exception as e:
        LOG.write("Exception parsing %s: %s\n"%(raxmlInfoFile, str(e))) 
HTML.write("<tr><td><b>%s</b></td><td>%.1f seconds</td></tr>\n"%("Total Job Duration", time()-starttime))
HTML.write("</table>\n\n")
if raxmlCommand and len(raxmlCommand):
    HTML.write("<h2>RAxML Command Line</h2>"+" ".join(raxmlCommand)+"\n")
else:
    HTML.write("<h2>RAxML Not Run</h2>\n")

HTML.write("<h2>Genome Statistics</h2>\n<table border='1'>\n")
HTML.write("<tr><th>GenomeId</th><th>Total Genes</th><th>Single Copy</th><th>Used</th><th>Name</th></tr>\n")
for genome in sorted(genesPerGenome, key=genesPerGenome.get):
    if genome not in genomeIds:
        LOG.write("%s not in genomeIds\n"%genome)
    elif genome not in singleCopyGenesPerGenome:
        LOG.write("%s not in singleCopyGenesPerGenome\n"%genome)
    elif genome not in genomeIdToName:
        LOG.write("%s not in genomeIdToName\n"%genome)
    else:
        HTML.write("<tr><td>%s</td><td>%d</td><td>%d</td><td>%d</td>"%(genome, genesPerGenome[genome], singleCopyGenesPerGenome[genome], filteredSingleCopyGenesPerGenome[genome]))
        if genome not in genomeIdToName:
            genomeIdToName[genome] = genome
        HTML.write("<td><a htarget='_blank' ref='https://patricbrc.org/view/Genome/%s'>%s</a></td></tr>\n"%(genome, genomeIdToName[genome]))
for genome in sorted(deletedGenomes):
    HTML.write("<tr><td>%s</td><td>0</td><td>0</td><td>0</td><td>Deleted</td></tr>\n"%(genome)) 

HTML.write("</table>\n\n")

if len(alignmentScore):
    HTML.write("<h2>Gene Family Statistics</h2>\n<table border='1'>\n")
    HTML.write("<p>Gene families are ranked by alignment score combining mean per-position variability, alignment length, and gappiness.</p>\n")
    pgfamProduct = patric_api.getProductsForPgfamsByN(alignmentScore.keys())
    statsToShow = ['num_pos', 'num_seqs', 'mean_squared_freq', 'prop_gaps']
    altHeadings = ['Align.<br>Length', 'Num<br>Seqs', 'Mean<br>Sqr Freq', 'Prop<br>Gaps']
    HTML.write("<tr><th>PGFam</th><th>Align.<br>Score</th><th>"+"</th><th>".join(altHeadings)+"</th><th>Used In<br>Analysis</th><th>Product</th></tr>\n")
    for pgfam in sorted(alignmentScore, key=alignmentScore.get, reverse=True): #proteinAlignments:
        stats = proteinAlignmentStats[pgfam]
        HTML.write("<tr><td>%s</td><td>%.2f</td>"%(pgfam, alignmentScore[pgfam]))
        for key in statsToShow:
            val = stats[key]
            if isinstance(val,int):
                HTML.write("<td>%d</td>"%stats[key])
            else:
                HTML.write("<td>%.3f</td>"%stats[key])
        HTML.write("<td>%s</td>\n"%str(pgfam in singleCopyPgfams))
        HTML.write("<td>%s</td></tr>\n"%pgfamProduct[pgfam])
    HTML.write("</table>\n")

if args.debugMode or len(singleCopyPgfams) < args.maxGenes:
    HTML.write("<h2>Strategies to Increase Single-Copy Gene Number</h2>\n")
    if len(singleCopyPgfams) < args.maxGenes:
        HTML.write("Number of single-copy genes (%d) was less than requested (%d).<br>\n"%(len(singleCopyPgfams), args.maxGenes))
    HTML.write("Examining the number of genes per genome in the 'Genome Statistics' table above may indicate incomplete or plasmid entries with few genes which may be removed.<br>\n")
    HTML.write("Criteria for calling single copy genes can be made more lenient by increasing Max Allowed Deletions and/or Max Allowed Duplications.<br>\n")
    #compute single copy counts for subsets of taxa, helps user identify optimal subsets for further analysis
    genomeSubsetSingleCopy = phylocode.countSingleCopyForGenomeSubsets(pgfamMatrix, genomeIds, maxAllowedDups = args.maxAllowedDups)
    if args.debugMode:
        LOG.write("len(genomeSubsetSingleCopy) = %d\n"%len(genomeSubsetSingleCopy))
    if len(genomeSubsetSingleCopy) > 0:
        HTML.write("<p>Omitting one of the following sets of genomes will provide the approximate boost in single-copy genes:</p>\n")
        HTML.write("<table border='1'>\n")
        HTML.write("<tr><th>scGene Boost</th><th>Genomes to Omit</th></tr>\n")
        maxToWrite = 5
        for subsetTuple in sorted(genomeSubsetSingleCopy, key=genomeSubsetSingleCopy.get, reverse=True):
            omissionSet = genomeIds - set(subsetTuple)
            HTML.write("<tr><td>%d</td><td>%s</td></tr>\n"%(genomeSubsetSingleCopy[subsetTuple], " ".join(sorted(omissionSet))))
            maxToWrite -= 1
            if not maxToWrite:
                break
        HTML.write("</table>\n")
HTML.write("</html>\n")
HTML.close()

LOG.write("output written to directory %s\n"%args.outputDirectory)
sys.stdout.write("\n")
detailsDirectory = "detail_files"
LOG.write("details files moved to subdirectory %s\n"%detailsDirectory)
if not os.path.isdir(detailsDirectory):
    if os.path.exists(detailsDirectory):
        os.remove(detailsDirectory)
    os.mkdir(detailsDirectory)
numMoved=0
for fn in filesToMoveToDetailsFolder:
    LOG.write("\t"+fn+"\t"+os.path.join(detailsDirectory, fn)+"\n")
    os.rename(fn, os.path.join(detailsDirectory, fn))
    numMoved += 1
LOG.write("files moved: %d\n"%numMoved)
filesToDelete.extend(glob.glob("RAxML*"))
if os.path.exists(phyloFileBase+".phy.reduced"):
    filesToDelete.append(phyloFileBase+".phy.reduced")
numDeleted = 0
for fn in filesToDelete:
    os.remove(fn)
    numDeleted += 1
    LOG.write("\t"+fn+"\n")
LOG.write("files deleted: %d\n"%numDeleted)
logfileName = os.path.basename(logfileName)
LOG.write("finally, will move this file, %s, to %s\n"%(logfileName, detailsDirectory))
LOG.close()
# finally, move the log file into the detailsDirectory
os.rename(logfileName, os.path.join(detailsDirectory, logfileName))
