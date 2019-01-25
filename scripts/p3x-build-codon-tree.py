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
parser.add_argument("--genomeIdsFile", metavar="file", type=str, nargs="*", help="file with PATRIC genome IDs, one per line (or first column of TSV)")
parser.add_argument("--genomeGroupName", metavar="name", type=str, nargs="*", help="name of user's genome group at PATRIC")
parser.add_argument("--genomeObjectFile", metavar="file", type=str, help="genome object (json file) to be added to ingroup")
parser.add_argument("--genomePgfamGeneFile", metavar="file", type=str, help="read geneIDs per PGFam per genome from this file")
parser.add_argument("--optionalGenomeIdsFile", metavar="file", type=str, help="optional genome ids, one per line (or first column of TSV)")
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
parser.add_argument("--deferRaxml", action='store_true', help="does not raxml but provides command file")
parser.add_argument("--writePgfamAlignments", action='store_true', help="create fasta alignment file per pgfam used for tree")
parser.add_argument("--outputDirectory", type=str, default=".", metavar="out_dir", help="directory for output, create if it does not exist")
parser.add_argument("--pathToFigtreeJar", type=str, metavar="jar_file", help="specify this to generate PDF graphic: java -jar pathToFigtreeJar -graphic PDF CodonTree.nex CodonTree.pdf")
parser.add_argument("--focusGenome", metavar="genome_id", type=str, help="genome to be highlighted in color in Figtree")
parser.add_argument("--debugMode", action='store_true', help="turns on more progress output to log file")
#parser.add_argument("--enableGenomeGenePgfamFileReuse", action='store_true', help="read genes and pgfams from stored file matching genomeIdsFile if it exists")
args = parser.parse_args()
starttime = time()

if not args.outputDirectory:
    args.outputDirectory=fileBase+"_dir/"
if not args.outputDirectory.endswith("/"):
    args.outputDirectory += "/"
if not os.path.exists(args.outputDirectory):
    os.makedirs(args.outputDirectory)

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

subprocess.check_call(['which', args.raxmlExecutable])
phylocode.checkMuscle()

ingroupIds = set() # keep unique
if args.genomeIdsFile:
    for filename in args.genomeIdsFile
        with open(filename) as F:
            for line in F:
                m = re.match(r"(\d+\.\d+)", line)
                if m:
                    ingroupIds.add(m.group(1))
    LOG.write("elapsed seconds = %f\n"%(time()-starttime))
    LOG.write("from %s got %d ingroupIds\n%s\n"%(",".join(args.genomeIdsFile), len(ingroupIds), "\t".join(ingroupIds)))
if args.genomeObjectFile:
    LOG.write("genome object file: %s\n"%args.genomeObjectFile)
LOG.flush()

if args.genomeGroupName:
    for groupName in args.genomeGroupName:
        LOG.write("requesting genome IDs for user group %s\n"%args.genomeGroupName)
        ids = patric_api.getGenomeGroupIds(args.genomeGroupName)
        LOG.write("got %d ids for %s\n"%(len(ids), args.genomeGroupName))
        ingroupIds.update(set(ids))

optionalGenomeIds = set()
if args.optionalGenomeIdsFile:
    with open(args.optionalGenomeIdsFile) as F:
        for line in F:
            m = re.match(r"(\d+\.\d+)", line)
            if m:
                optionalGenomeIds.add(m.group(1))
LOG.write("got %d optionalGenomeIds\n%s\n"%(len(optionalGenomeIds), "\t".join(optionalGenomeIds)))
LOG.flush()

if len(ingroupIds) + len(optionalGenomeIds) < 4:
    LOG.write("too few genomes to build a tree with: %d"%(len(ingroupIds)+len(optionalGenomeIds)))
    sys.exit(1)

# figure out a base name for ouput files
if args.genomeIdsFile:
    fileBase = os.path.basename(args.genomeIdsFile)
    fileBase = re.sub("\..*", "", fileBase)
elif args.genomeGroupName:
    fileBase = args.genomeGroupName
else:
    fileBase = "codon_tree"


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


# this is where we gather the list of Pgfam genes for each ingroup genome ID
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
    ingroupIds.add(genomeObject_genomeId)
    args.focusGenome = genomeObject_genomeId
    LOG.write("parsed json file %s"%args.genomeObjectFile)
    LOG.flush()

allGenomeIds = ingroupIds
allGenomeIds.update(optionalGenomeIds)
genomesWithData = set()
for pgfam in pgfamMatrix:
    genomesWithData.update(set(pgfamMatrix[pgfam].keys()))
genomesWithoutData = allGenomeIds - genomesWithData
if len(genomesWithoutData):
    pgfamMatrix = patric_api.getPgfamGenomeMatrix(genomesWithoutData, pgfamMatrix)

with open(args.outputDirectory+fileBase+".pgfamMatrix.txt", 'w') as F:
   patric_api.writePgfamGenomeMatrix(pgfamMatrix, F)

LOG.write("got pgfams for genomes, len=%d\n"%len(pgfamMatrix))
LOG.flush()

LOG.write("allowing %d genomes missing per PGfam of ingroup (out of %d total)\n"%(args.maxGenomesMissing, len(ingroupIds)))
if args.maxGenomesMissing >= len(ingroupIds)-4:
    raise Exception("getSingleCopyPgfams: maxGenomesMissing too large: %d"%args.maxGenomesMissing)

# call to getSingleCopyPgfams uses ingroup taxa, optional genomes are not involved in selecting single copy pgfams
singleCopyPgfams = phylocode.selectSingleCopyPgfams(pgfamMatrix, ingroupIds, requiredGenome=args.focusGenome, maxGenomesMissing=args.maxGenomesMissing, maxAllowedDups=args.maxAllowedDups)

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
                SeqIO.write(codonSeqRecords[:2], LOG, "fasta")
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
        SeqIO.write(proteinAlignments[pgfam], args.outputDirectory+pgfam+".faa", "fasta")

with open(args.outputDirectory + "PgfamAlignmentStats.txt", "w") as F:
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
    nexusOutfileBase = os.path.join(args.outputDirectory, "CodonTree")
    nexusFilesWritten = phylocode.generateNexusFile(originalNewick, nexusOutfileBase, nexus_template = nexus_template_file, align_tips = "both", focus_genome = args.focusGenome, genomeIdToName=genomeIdToName)
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
        
OUT = open(args.outputDirectory+"CodonTree.stats", 'w')
OUT.write("Statistics for CodonTree\n")
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

