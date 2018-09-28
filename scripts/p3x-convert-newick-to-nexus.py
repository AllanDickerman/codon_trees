import sys
import re
import os.path
import argparse
import subprocess
Codon_trees_lib_path = os.path.dirname(sys.argv[0]).replace("scripts", "lib")
sys.path.append(Codon_trees_lib_path)
import patric_api
import phylocode

parser = argparse.ArgumentParser()
parser.add_argument("newickTreeFile", metavar="tree file", type=str, help="newick tree file labeled with PATRIC genome IDs")
parser.add_argument("--modelNexusFile", metavar="nexus", type=str, help="genome object (json file) to be added to ingroup")
parser.add_argument("--focusGenome", metavar="genome_id", type=str, help="genome to be highlighted in color in Figtree")
parser.add_argument("--pathToFigtreeJar", type=str, metavar="jar_file", help="specify this to generate PDF graphic")
parser.add_argument("--debugMode", action='store_true', help="more progress output")
args = parser.parse_args()

figtreeParams={}
if not args.modelNexusFile:
    args.modelNexusFile = os.path.join(Codon_trees_lib_path, "figtree.nex")
# read a model figtree nexus file
if os.path.exists(args.modelNexusFile):
    figtreeParams = phylocode.readFigtreeParameters(args.modelNexusFile)

newick = open(args.newickTreeFile).read()
genomeIds = re.findall("[(,]([^(,):]+)[,:)]", newick)
sys.stderr.write("\n".join(genomeIds))

genomeIdToName = {}
for genomeId, genomeName in patric_api.getNamesForGenomeIds(genomeIds):
    genomeIdToName[genomeId] = genomeName+" "+genomeId

outfileBase = re.sub(".n[ewick]+$", "", args.newickTreeFile)
nexusOutfileName = outfileBase+"_figtree.nex"
nexusOut = open(nexusOutfileName, "w")
phylocode.writeTranslatedNexusTree(nexusOut, newick, genomeIdToName, figtreeParameters=figtreeParams, highlightGenome=args.focusGenome)
nexusOut.close()
sys.stderr.write("nexus file written to %s\n"%nexusOutfileName)

numTaxa = len(genomeIds)
if not args.pathToFigtreeJar:
    if os.path.exists(os.path.join(Codon_trees_lib_path, "figtree.jar")):
        args.pathToFigtreeJar = os.path.join(Codon_trees_lib_path, "figtree.jar")
if args.pathToFigtreeJar:
    if os.path.exists(args.pathToFigtreeJar):
        if args.debugMode:
            sys.stderr.write("found figtree.jar at %s\n"%args.pathToFigtreeJar)
    else:
        sys.stderr.write("Could not find figtree.jar")
        args.pathToFigtreeJar = None
if args.pathToFigtreeJar:
        figtreePdfName = outfileBase+"_figtree.pdf"
        #def generateFigtreeImage(nexusFile, outfileName, numTaxa, figtreeJarFile, imageFormat="PDF")
        phylocode.generateFigtreeImage(nexusOutfileName, figtreePdfName, len(genomeIds), args.pathToFigtreeJar)
        # Now write a version of tree with tip labels aligned (shows bootstap support better)
        if True: # possibly gate this by a parameter
            figtreeParams['rectilinearLayout.alignTipLabels'] = 'true'
            nexusOut = open(nexusOutfileName, "w")
            phylocode.writeTranslatedNexusTree(nexusOut, newick, genomeIdToName, figtreeParameters=figtreeParams, highlightGenome=args.focusGenome)
            nexusOut.close()
            figtreePdfName = outfileBase+"_figtree_tipLabelsAligned.pdf"
            phylocode.generateFigtreeImage(nexusOutfileName, figtreePdfName, len(genomeIds), args.pathToFigtreeJar)
            
