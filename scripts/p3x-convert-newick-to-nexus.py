import sys
import re
import os.path
import argparse
import patric_api
import phylocode

parser = argparse.ArgumentParser()
parser.add_argument("newickTreeFile", metavar="tree file", type=str, help="newick tree file labeled with PATRIC genome IDs")
parser.add_argument("--nexusTemplateFile", metavar="nexus", type=str, help="default values for Figtree visualization")
parser.add_argument("--focusGenome", metavar="genome_id", type=str, help="genome to be highlighted in color in Figtree")
parser.add_argument("--pathToFigtreeJar", type=str, metavar="jar_file", help="specify this to generate PDF graphic")
parser.add_argument("--outputDirectory", type=str, default=".", metavar="DIR", help="where output should go")
parser.add_argument("--alignTips", choices=['yes', 'no', 'both'], type=str, default='both', help="how to draw tips")
parser.add_argument("--debugMode", action='store_true', help="more progress output")
args = parser.parse_args()

logfile = os.path.basename(sys.argv[0])
logfile = re.sub('.py', '.nex', logfile)
LOG = open(logfile, 'w')
phylocode.LOG = LOG
if not args.nexusTemplateFile:
    for dirname in sys.path: # should be in .../codon_trees/lib
        if os.path.isfile(os.path.join(dirname, "figtree.nex")):
            args.nexusTemplateFile = os.path.join(dirname, "figtree.nex")
# read a model figtree nexus file
figtreeParams={}
if os.path.exists(args.nexusTemplateFile):
    LOG.write("Found figtree template file: %s\n"%args.nexusTemplateFile)
    LOG.flush()
    figtreeParams = phylocode.readFigtreeParameters(args.nexusTemplateFile)

newick = open(args.newickTreeFile).read()
genomeIds = re.findall("[(,]([^(,):]+)[,:)]", newick)
LOG.write("\n".join(genomeIds))

nexusOutfileBase = os.path.basename(args.newickTreeFile)
nexusOutfileBase = re.sub(".n[ewick]+$", "", nexusOutfileBase)
nexusOutfileBase = nexusOutfileBase+"_figtree"
nexusOutfileBase = os.path.join(args.outputDirectory, nexusOutfileBase)

nexusFilesWritten = phylocode.generateNexusFile(newick, nexusOutfileBase, nexus_template = args.nexusTemplateFile, align_tips = args.alignTips, focus_genome = args.focusGenome, genomeIdToName=None)
LOG.write("nexus file written to %s\n"%(", ".join(nexusFilesWritten)))

numTaxa = len(genomeIds)
if args.pathToFigtreeJar:
    if not os.path.exists(args.pathToFigtreeJar):
        message = "specified figtree.jar does not exist: %s\n"%args.pathToFigtreeJar
        LOG.write(message)
        raise Exception(message)
    for nexusFile in nexusFilesWritten:
        figtreePdfName = re.sub(".nex", "_figtree.pdf", nexusFile)
        phylocode.generateFigtreeImage(nexusFile, figtreePdfName, len(genomeIds), args.pathToFigtreeJar)
