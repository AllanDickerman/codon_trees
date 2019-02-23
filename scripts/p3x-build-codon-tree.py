import sys
import re
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
parser.add_argument("--bootstrapReps", metavar="#", type=int, default=0, help="number of raxml 'fast boostrap' replicates")
parser.add_argument("--maxGenomesMissing", metavar="#", type=int, default=0, help="genomes allowed to lack a member of any homolog group")
parser.add_argument("--maxAllowedDups", metavar="maxDups", type=int, default=0, help="duplicated gene occurrences allowed within homolog group")
parser.add_argument("--endGapTrimThreshold", metavar="maxPropGaps", type=float, default=0.5, help="stringency of end-gap trimming, lower for less trimming")
parser.add_argument("--raxmlExecutable", metavar="program_name", type=str, default="raxml", help="program to call, possibly with path")
parser.add_argument("--rateModel", metavar="rateModel", type=str, choices = ['CAT', 'GAMMA'], default="CAT", help="variable rate category model CAT|GAMMA")
parser.add_argument("--proteinModel", metavar="substModel", type=str, default="WAGF", help="raxml protein substitution model")
parser.add_argument("--analyzeCodons", action='store_true', help="analyze only codons (ignore amino acids)")
parser.add_argument("--analyzeProteins", action='store_true', help="analyze only amino acids")
parser.add_argument("--threads", metavar="T", type=int, default=2, help="number of threads for raxml")
parser.add_argument("--deferRaxml", action='store_true', help="does not raxml but provides command file")
parser.add_argument("--writePgfamAlignments", action='store_true', help="create fasta alignment file per pgfam used for tree")
parser.add_argument("--outputDirectory", type=str, metavar="out_dir", help="directory for output, create if it does not exist")
parser.add_argument("--pathToFigtreeJar", type=str, metavar="jar_file", help="specify this to generate PDF graphic: java -jar pathToFigtreeJar -graphic PDF CodonTree.nex CodonTree.pdf")
parser.add_argument("--focusGenome", metavar="genome_id", type=str, help="genome to be highlighted in color in Figtree")
parser.add_argument("--debugMode", action='store_true', help="turns on more progress output to log file")
parser.add_argument("--authenticateFile", type=str, metavar="file(optional)", nargs="?", default=os.path.expanduser("~/.patric_token"), help="authenticate patric user by token file")
parser.add_argument("--authenticateEnv", action='store_true', help="authenticate using environment variable KB_AUTH_TOKEN")
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
if args.authenticateFile:
    patric_api.authenticateByFile(args.authenticateFile)
else:
    patric_api.authenticateByEnv()

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
    for groupName in args.genomeGroupName:
        LOG.write("requesting genome IDs for user group %s\n"%args.genomeGroupName)
        ids = patric_api.getGenomeGroupIds(args.genomeGroupName)
        LOG.write("got %d ids for %s\n"%(len(ids), args.genomeGroupName))
        genomeIds.update(set(ids))

if args.optionalGenomeIdsFile:
    with open(args.optionalGenomeIdsFile) as F:
        for line in F:
            m = re.match(r"(\d+\.\d+)", line)
            if m:
                optionalGenomeIds.add(m.group(1))
    LOG.write("got %d optionalGenomeIds\n%s\n"%(len(optionalGenomeIds), "\t".join(optionalGenomeIds)))
    LOG.flush()

preflightTests.append(("Need at least 4 genomes to build a tree", len(genomeIds) + len(optionalGenomeIds) >= 4))

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
    pgfamMatrix = patric_api.getPgfamGenomeMatrix(genomesWithoutData, pgfamMatrix)

with open(os.path.join(args.outputDirectory, args.outputBase+".pgfamMatrix.txt"), 'w') as F:
   patric_api.writePgfamGenomeCountMatrix(pgfamMatrix, F)

LOG.write("got pgfams for genomes, len=%d\n"%len(pgfamMatrix))
LOG.flush()

if args.maxGenomesMissing:
    LOG.write("allowing %d genomes missing per PGfam (out of %d total)\n"%(args.maxGenomesMissing, len(genomeIds)))
preflightTests.append(("Num genomes - maxGenomesMissing >= 4", len(genomeIds) - args.maxGenomesMissing >= 4))

# call to getSingleCopyPgfams uses main genome, optional genomes are not involved in selecting single copy pgfams
singleCopyPgfams = phylocode.selectSingleCopyPgfams(pgfamMatrix, genomeIds, requiredGenome=args.focusGenome, maxGenomesMissing=args.maxGenomesMissing, maxAllowedDups=args.maxAllowedDups)

LOG.write("got single copy pgfams, num=%d\n"%len(singleCopyPgfams))
if len(singleCopyPgfams) > args.maxGenes:
    singleCopyPgfams=singleCopyPgfams[0:args.maxGenes]
    LOG.write("\tselecting top single-family genes: %d\n"%len(singleCopyPgfams))
LOG.flush()
with open(os.path.join(args.outputDirectory, args.outputBase+".singleCopyPgfams.txt"), 'w') as F:
    for pgfam in singleCopyPgfams:
        F.write(pgfam+"\n")

preflightTests.append(("single copy pgfams > 0", len(singleCopyPgfams) > 0))
## perform preflight test
with open(os.path.join(args.outputDirectory, args.outputBase+".preflight"), 'w') as F:
    passed = True
    for test in preflightTests:
        passed &= test[1]
        F.write(test[0] + "\t" + str(test[1]) + "\n")
    F.write("All tests passed\t"+str(passed)+"\n")
    if not passed:
        F.write("exiting\n")
        F.close()
        sys.exit(1)

proteinAlignments = {}
codonAlignments = {}
alignedTaxa=set()
for pgfamId in singleCopyPgfams:
    geneIdSet = set()
    for genome in pgfamMatrix[pgfamId]:
        geneIdSet.update(set(pgfamMatrix[pgfamId][genome]))
    proteinFasta = patric_api.getProteinFastaForPatricIds(geneIdSet)
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
                SeqIO.write(codonAlignment[pgfamId][:2], LOG, "fasta")
    except Exception as e:
        LOG.write("Exception aligning codons: %s\n"%str(e))
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

# generate hopefully unique output file name base
phyloFileBase = args.outputBase + "_%dtaxa"%(numTaxa)
if args.analyzeCodons:
    phyloFileBase += "_%dcds"%len(codonAlignments)
if args.analyzeProteins:
    phyloFileBase += "_%dproteins"%len(proteinAlignments)

# change to output directory to simplify file naming
os.chdir(args.outputDirectory)

proteinPositions=0
codonPositions = 0
# write the genes included in each homology group (and those paralogs excluded)
with open(phyloFileBase+".pgfamsAndGenesIncludedInAlignment.txt", 'w') as F:
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

with open("PgfamAlignmentStats.txt", "w") as F:
    first = True
    for pgfam in proteinAlignments:
        stats = phylocode.calcAlignmentStats(proteinAlignments[pgfam])
        if first:
            F.write("PGFam\t"+"\t".join(sorted(stats.keys()))+"\n")
            first = False
        F.write(pgfam)
        for key in sorted(stats):
            val = stats[key]
            if isinstance(val,int):
                F.write("\t%d"%stats[key])
            else:
                F.write("\t%.6f"%stats[key])
        F.write("\n")

# finally, output concatenated protein and/or DNA alignment and partitions and raxml command to appropriate files
raxmlCommand=''
if args.analyzeProteins and args.analyzeCodons:
    phylocode.outputCodonsProteinsPhylip(codonAlignments, proteinAlignments, phyloFileBase+".phy")
    with open(phyloFileBase+".partitions", 'w') as PartitionFile:
        for i in range(1,4):
            PartitionFile.write("DNA, codon%d = %d-%d\\3\n"%(i, i, codonPositions))
        PartitionFile.write("%s, proteins = %d-%d\n"%(args.proteinModel, codonPositions+1, codonPositions+proteinPositions))
    raxmlCommand = [args.raxmlExecutable, "-s", phyloFileBase+".phy", "-n", phyloFileBase, "-m",  "GTR%s"%args.rateModel, "-q",  phyloFileBase+".partitions",  "-p", "12345", "-T", str(args.threads)]

elif args.analyzeCodons:
    phylocode.writeConcatenatedAlignmentsPhylip(codonAlignments, phyloFileBase+".phy")
    with open(phyloFileBase+".partitions", 'w') as PartitionFile:
        for i in range(1,4):
            PartitionFile.write("DNA, codon%d = %d-%d\\3\n"%(i, i, codonPositions))
    raxmlCommand = [args.raxmlExecutable, "-s", phyloFileBase+".phy", "-n", phyloFileBase, "-m",  "GTR%s"%args.rateModel, "-q",  phyloFileBase+".partitions",  "-p", "12345", "-T", str(args.threads)]

elif args.analyzeProteins:
    phylocode.writeConcatenatedAlignmentsPhylip(proteinAlignments, phyloFileBase+".phy")
    raxmlCommand = [args.raxmlExecutable, "-s", phyloFileBase+".phy", "-n", phyloFileBase, "-m",  "PROT%s%s"%(args.rateModel, args.proteinModel), "-p", "12345", "-T", str(args.threads)]
raxmlCommand.extend(["-e", "1.0"]) # limit on precision, faster than default 0.1

if args.bootstrapReps > 0:
    raxmlCommand.extend(["-f", "a", "-x", "12345", "-N", str(args.bootstrapReps)]) 
with open(phyloFileBase+".raxmlCommand.sh", 'w') as F:
    F.write(" ".join(raxmlCommand)+"\n")

if not args.deferRaxml:
    #remove RAxML files that clash in name, their existence blocks raxml from running
    for fl in glob.glob("RAxML_*"+phyloFileBase):
        os.remove(fl)
    proc = subprocess.Popen(raxmlCommand)
    proc.wait()
    LOG.write("raxml completed: elapsed seconds = %f\n"%(time()-starttime))
    LOG.flush()
    genomeIdToName = {}
    for genomeId, genomeName in patric_api.getNamesForGenomeIds(allGenomeIds):
        genomeIdToName[genomeId] = genomeName+" "+genomeId
    if genomeObject:
        genomeIdToName[genomeObject_genomeId] = genomeObject_name+" "+genomeObject_genomeId
    originalNewick = ""
    raxmlNewickFileName = "RAxML_bestTree."+phyloFileBase
    if args.bootstrapReps > 0:
        raxmlNewickFileName = "RAxML_bipartitions."+phyloFileBase
    F = open(raxmlNewickFileName)
    originalNewick = F.read()
    F.close()
    renamedNewick = phylocode.relabelNewickTree(originalNewick, genomeIdToName)
    renamedNewickFile = phyloFileBase+"_treeWithGenomeNames.nwk"
    F = open(renamedNewickFile, 'w')
    F.write(renamedNewick)
    F.close()
    LOG.write("codonTree newick relabeled with genome names written to "+renamedNewickFile+"\n")
    LOG.flush()

    # Search for the template file figtree.nex in the same directories
    # as our library code, somewhere in sys.path.
    nexus_template_file = None
    for dirname in sys.path:
        if os.path.isfile(os.path.join(dirname, "figtree.nex")):
            nexus_template_file = os.path.join(dirname, "figtree.nex")
    if os.path.exists(nexus_template_file):
        figtreeParams = phylocode.readFigtreeParameters(nexus_template_file)
        LOG.write("Found figtree template file: %s\n"%nexus_template_file)
    else:
        figtreeParams = {}
        LOG.write("Could not find valid template nexus file.\n")
        LOG.flush()
    nexusFilesWritten = phylocode.generateNexusFile(originalNewick, phyloFileBase, nexus_template = nexus_template_file, align_tips = "both", focus_genome = args.focusGenome, genomeIdToName=genomeIdToName)
    LOG.write("nexus file written to %s\n"%(", ".join(nexusFilesWritten)))

    if not (args.pathToFigtreeJar and os.path.exists(args.pathToFigtreeJar)):
        LOG.write("Could not find valid path to figtree.jar\n")
        args.pathToFigtreeJar = None
    if args.pathToFigtreeJar:
        if os.path.exists(args.pathToFigtreeJar):
            LOG.write("found figtree.jar at %s\n"%args.pathToFigtreeJar)
            for nexusFile in nexusFilesWritten:
                figtreePdfName = re.sub(".nex", ".pdf", nexusFile)
                phylocode.generateFigtreeImage(nexusFile, figtreePdfName, len(allGenomeIds), args.pathToFigtreeJar)
                LOG.write("created figtree figure: %s\n"%figtreePdfName)
        else:
            message = "specified figtree.jar does not exist: %s\n"%args.pathToFigtreeJar
            LOG.write(message)

LOG.write(strftime("%a, %d %b %Y %H:%M:%S", localtime(time()))+"\n")
LOG.write("Total job duration %d seconds\n"%(time()-starttime))
        
OUT = open(phyloFileBase+"_codontree_analysis.stats", 'w')
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
LOG.write("output written to directory %s\n"%args.outputDirectory)
LOG.close()
sys.stdout.write("\n") 

