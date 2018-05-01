import sys
import re
import os.path
import argparse
import inspect
from p3_allan import patric_api
from p3_allan import phylocode

parser = argparse.ArgumentParser()
parser.add_argument("newickTreeFile", metavar="tree file", type=str, help="newick tree file labeled with PATRIC genome IDs")
parser.add_argument("--modelNexusFile", metavar="nexus", type=str, help="genome object (json file) to be added to ingroup")
parser.add_argument("--focusGenome", metavar="genome_id", type=str, help="genome to be highlighted in color in Figtree")
args = parser.parse_args()

figtreeParams={}
if not args.modelNexusFile:
    pathToInstall = os.path.dirname(inspect.getfile(phylocode))
    args.modelNexusFile = pathToInstall+"/figtree.nex"
# read a model figtree nexus file
if os.path.exists(pathToInstall+"/figtree.nex"):
    figtreeParams = phylocode.readFigtreeParameters(args.modelNexusFile)

newick = open(args.newickTreeFile).read()
genomeIds = re.findall("[(,]([^(,):]+)[,:)]", newick)
sys.stderr.write("\n".join(genomeIds))

genomeIdToName = {}
for genomeId, genomeName in patric_api.getNamesForGenomeIds(genomeIds):
    genomeIdToName[genomeId] = genomeName+" "+genomeId

nexusOutfileName = re.sub(".n[ewick]+$", "", args.newickTreeFile)+"_figtree.nex"
nexusOut = open(nexusOutfileName, "w")
phylocode.writeTranslatedNexusTree(nexusOut, newick, genomeIdToName, figtreeParameters=figtreeParams, highlightGenome=args.focusGenome)
nexusOut.close()
sys.stderr.write("nexus file written to %s\n"%nexusOutfileName)
