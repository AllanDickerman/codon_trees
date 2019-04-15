import sys
import re
import subprocess
import os.path
import warnings
from Bio import BiopythonWarning
from Bio.Alphabet import IUPAC
from Bio import AlignIO
from Bio import SeqIO
from Bio import Alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from collections import defaultdict
import patric_api
import StringIO
from Bio import BiopythonExperimentalWarning
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    from Bio import codonalign

Debug = False #shared across functions defined here
LOG = sys.stderr

def which(program):
    """ Can use to test for existence of needed files on path
    based on code from https://stackoverflow.com/users/20840/jay """
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

def getPgfamDistribution(genomeGenePgfamList):
    """ Given list of genome-gene-pgfam tuples, 
        tabulate counts per genome per pgfam, 
        num duplications per pgfam, 
        num single-copy.
        This information can be used to identify PGFams suitable for phylogenetic inference for these genomes.
        Alternatively, can use to identify genomes with sparse single-copy coverage to remove to create denser phylogenetic matrix. 
    """
    ggpMat = {} # genome-gene-pgfam matrix (really just a dictionary)
    for row in genomeGenePgfamList:
        genome, gene, pgfam = row
        if pgfam not in ggpMat:
            ggpMat[pgfam] = {}
        if genome not in ggpMat[pgfam]:
            ggpMat[pgfam][genome] = []
        ggpMat[pgfam][genome].append(gene)
    return ggpMat

def selectSingleCopyPgfams(pgfamMatrix, genomeIdList, requiredGenome=None, maxGenomesMissing=0, maxAllowedDups=0):
    # given a genome-gene-pgfam matrix
    # find the set of pgfam_ids which satisfy the maxGenomesMissing & maxAllowedDups criteria for specified genomes
    pgfamScore = {}
    scPgfamsIfGenomeOmitted = {}
    if not len(genomeIdList) > 3:
        raise Exception("getSingleCopyPgfams: number of genome IDs too low: %d"%len(genomeIdList))
    minGenomesPresent = len(genomeIdList) - maxGenomesMissing
    if minGenomesPresent <= 0:
        raise Exception("getSingleCopyPgfams: maxGenomesMissing too large: %d"%maxGenomesMissing)
    for pgfam in pgfamMatrix:
        if requiredGenome and requiredGenome not in pgfamMatrix[pgfam]:
            continue # skip pgfam if it does not include required genome
        pgfamGenomeCount = 0
        numDups = 0
        for genome in genomeIdList:
            if genome in pgfamMatrix[pgfam]:
                pgfamGenomeCount += 1
                numDups += len(pgfamMatrix[pgfam][genome]) - 1
        if pgfamGenomeCount >= minGenomesPresent and numDups <= maxAllowedDups:
            pgfamScore[pgfam] = pgfamGenomeCount - numDups # combined score of both factors, for sorting
    suitablePgfamList = sorted(pgfamScore, key=pgfamScore.get, reverse=True)
    return suitablePgfamList

def countSingleCopyForGenomeSubsets(pgfamMatrix, genomeIds, maxAllowedDups=0):
    """ Return dict of how many single copy genes would be found if 
    certain taxa were omitted from the analysis. 
    Takes into consideration maxAllowedDups, but not maxGenomesMissing.
    Return dict of genome tuples to count of single copy genes IF that set is omitted.
    Record subsets down to half the original set size.
    """
    scForSubset = {}
    minSetSize = int(len(genomeIds) * 0.75)
    for pgfam in pgfamMatrix:
        numDups = 0
        present = []
        missing = 0
        for genome in sorted(pgfamMatrix[pgfam]):
            x = len(pgfamMatrix[pgfam][genome])
            if x == 0:
                missing += 1
            else:
                present.append(genome)
                if x > 1:
                    numDups += (x - 1)
        if False and Debug:
            LOG.write("countSingleCopyForGenomeSubsets, pgfam=%s, missing=%d, present=%d, nd=%d\n"%(pgfam, missing, len(present), numDups))

        if numDups <= maxAllowedDups and (len(present) >= minSetSize) and (len(present) < len(genomeIds)):
            presentTuple = tuple(present)
            if presentTuple not in scForSubset:
                scForSubset[presentTuple] = 0
            scForSubset[presentTuple] += 1
    return scForSubset

def getGenesForPgfams(genomeGenePgfam, genomeIdList, singleCopyPgfams):
    genesForPgfams={}
    for pgfam in singleCopyPgfams:
        genesForPgfams[pgfam] = []
    for row in genomeGenePgfam:
        genome, gene, pgfam = row
        if genome in genomeIdList and pgfam in genesForPgfams:
            genesForPgfams[pgfam].append(gene)
    return genesForPgfams

def relabelNewickTree(newick, labelDict):
# newick file has tip labels eg genomIDs
# return value has replaced these according to dictionary passed as labelDict
    retval = newick
    for oldLabel in labelDict:
        retval = re.sub("\\b"+oldLabel+"\\b", '"'+labelDict[oldLabel]+'"', retval)
    return retval

def writeTranslatedNexusTree(outFile, newick, labelDict, figtreeParameters=None, highlightGenome=None):
    # intended to use with FigTree to generate figures, but can be used without a figtreeParameters object
#write taxa block
    taxonIds = re.findall(r"[(,]([^:(),]+)[,:)]", newick)
    outFile.write("#NEXUS\nbegin taxa;\n\tdimensions ntax=%d;\n\ttaxlabels\n"%len(taxonIds))
    for tax in taxonIds:
        if tax not in labelDict:
            labelDict[tax] = tax
        outFile.write("\t\"%s %s\""%(labelDict[tax], tax))
        if tax == highlightGenome:
            outFile.write("[&!color=#ff0000]")
        outFile.write("\n")
    outFile.write(";\nend;\n\n")
# prefix support values with "[&label="
    newick = re.sub(r"\)(\d+):", r")[&support=\1]:", newick)
# write trees block
    outFile.write("begin trees;\n")
    if labelDict:
        outFile.write("\ttranslate\n")
        for tax in taxonIds:
            translate = tax
            if tax in labelDict:
                translate = labelDict[tax]
            outFile.write("\t\t%s \"%s %s\",\n"%(tax, translate, tax))
        outFile.write("\t;\n")
    outFile.write("\ttree one = [&U] %s\n"%newick)
    outFile.write("\nend;\n\n")
# write figtree block
    if figtreeParameters:
        figtreeParameters["nodeLabels.isShown"]="true"
        figtreeParameters["nodeLabels.displayAttribute"]="\"support\""
        outFile.write("begin figtree;\n")
        for param in sorted(figtreeParameters):
            outFile.write("\tset %s=%s;\n"%(param, str(figtreeParameters[param])))
        outFile.write("\nend;\n\n")
    return

def readFigtreeParameters(filename):
    infile = open(filename)
    retval = {}
    inFigtreeBlock = False
    for line in infile:
        if re.search("begin figtree", line, re.IGNORECASE):
            inFigtreeBlock = True
        if re.search(r"^end;", line):
            inFigtreeBlock = False
        if inFigtreeBlock:
            m = re.search(r"set\s+(\S.*)=(\S.*);", line, re.IGNORECASE)
            if m:
                retval[m.group(1)] = m.group(2)
    return retval

def generateNexusFile(newick, outfileBase, nexus_template = None, align_tips = "both", focus_genome = None, genomeIdToName=None):
    figtreeParams={}
    if not nexus_template:
        LOG.write("Look for figtree.nex template in sys.path directories\n")
        for dirname in sys.path: # should be in .../codon_trees/lib
            if os.path.isfile(os.path.join(dirname, "figtree.nex")):
                nexus_template_file = os.path.join(dirname, "figtree.nex")
                LOG.write("Found figtree template file: %s\n"%nexus_template)
    # read a model figtree nexus file
    if nexus_template and os.path.exists(nexus_template):
        LOG.write("Read figtree template file: %s\n"%nexus_template)
        figtreeParams = readFigtreeParameters(nexus_template)
    genomeIds = re.findall("[(,]([^(,):]+)[,:)]", newick)
    if not genomeIdToName:
        genomeIdToName = patric_api.getNamesForGenomeIds(genomeIds)
    figtreeParams["trees.rooting"]="true"
    figtreeParams["trees.rootingType"]="Midpoint"
    figtreeParams["trees.order"]="true"
    figtreeParams["trees.orderType"]="increasing"
    figtreeParams["tipLabels.fontName"]="sanserif"
    figtreeParams["tipLabels.fontSize"]=14
    figtreeParams["tipLabels.fontStyle"]=0
    figtreeParams["tipLabels.isShown"]="true"
    filesWritten=[]
    if align_tips in ("no", "both"):
        figtreeParams['rectilinearLayout.alignTipLabels'] = 'false'
        nexusOut = open(outfileBase+".nex", "w")
        writeTranslatedNexusTree(nexusOut, newick, genomeIdToName, figtreeParameters=figtreeParams, highlightGenome=focus_genome)
        nexusOut.close()
        filesWritten.append(outfileBase+".nex")
    if align_tips in ("yes", "both"):
        figtreeParams['rectilinearLayout.alignTipLabels'] = 'true'
        nexusOut = open(outfileBase+"_tipsAligned.nex", "w")
        writeTranslatedNexusTree(nexusOut, newick, genomeIdToName, figtreeParameters=figtreeParams, highlightGenome=focus_genome)
        nexusOut.close()
        filesWritten.append(outfileBase+"_tipsAligned.nex")
    return filesWritten

def generateFigtreeImage(nexusFile, numTaxa=0, figtreeJar=None, imageFormat="PDF"):
    if Debug:
        LOG.write("generateTreeFigure(%s, %d, figtreeJarFile=%s, imageFormat=%s)\n"%(nexusFile, numTaxa, figtreeJar, imageFormat))
    if imageFormat not in ('PDF', 'SVG', 'PNG', 'JPEG'):
        raise Exception("imageFormat %s not in ('PDF', 'SVG', 'PNG', 'JPEG')"%imageFormat)
    imageFileName = nexusFile
    imageFileName = nexusFile.replace(".nex", ".")
    imageFileName += imageFormat.lower()
    if figtreeJar:
        figtreeCommand = ['java',  '-jar', figtreeJar, '-graphic', imageFormat]
    else:
        # assume a figtree executable is on the path
        figtreeCommand = ['figtree', '-graphic', imageFormat]
    if numTaxa > 40:
        height = 600 + 15 * (numTaxa - 40) # this is an empirical correction factor to avoid taxon name overlap
        figtreeCommand.extend(['-height', str(int(height))])
    figtreeCommand.extend([nexusFile, imageFileName])
    if Debug:
        LOG.write("running this command:\n%s\n"%" ".join(figtreeCommand))
    subprocess.call(figtreeCommand, stdout=LOG)
    return imageFileName

def checkMuscle():
    subprocess.check_call(['which', 'muscle'])

def alignSeqRecordsMuscle(seqRecords):
    muscleProcess = subprocess.Popen(['muscle', '-quiet'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    SeqIO.write(seqRecords, muscleProcess.stdin, 'fasta')
    muscleProcess.stdin.close()
    alignment = AlignIO.read(muscleProcess.stdout, "fasta", alphabet=seqRecords[0].seq.alphabet)
    alignment.sort()
    return(alignment)

def calcAlignmentStats(alignment):
    # analyze a BioPython MultipleSeqAlignment object to describe conservation levels
    stats = {}
    stats['num_pos'] = alignment.get_alignment_length()
    stats['num_seqs'] = len(alignment)
    numGaps = 0
    numNonGaps = 0
    sumSquaredFreq = 0
    for pos in range(0,stats['num_pos']):
        stateCount = {}
        states = alignment[: , pos]
        for residue in states:
            if residue == '-':
                numGaps += 1
                continue
            numNonGaps += 1
            if residue not in stateCount:
                stateCount[residue] = 0
            stateCount[residue] += 1
        for residue in stateCount:
            freq = stateCount[residue]/float(stats['num_seqs'])
            sumSquaredFreq += freq*freq
    stats['sum_squared_freq'] = sumSquaredFreq
    stats['mean_squared_freq'] = sumSquaredFreq/stats['num_pos']
    stats['gaps'] = numGaps
    stats['prop_gaps'] = numGaps/(float(numGaps+numNonGaps))
    #stats['non_gaps'] = numNonGaps
    return stats

def suggestAlignmentDeletions(alignment):
    """ 
    analyze a BioPython MultipleSeqAlignment object to calculate how many gaps could be avoided by omitting sequences
    return dict of tuple of seqIds to improvement score of deleting that set of seqs
    """
    delsets = {}
    for pos in range(0, alignment.get_alignment_length()):
        states = alignment[: , pos]
        gaps = 0.0
        nongap_seqs = []
        for i, residue in enumerate(states):
            if residue == '-':
                gaps += 1
            else:
                nongap_seqs.append(alignment[i].id)
        if gaps and gaps/len(alignment) >= 0.75: # a gap-heavy position
            key = tuple(nongap_seqs)
            if key not in delsets:
                delsets[key] = 0
            delsets[key] += 1
    return delsets

def calcSumAlignmentDistance(alignment, querySeq):
    sumDist = 0
    for record in alignment:
        for posPair in zip(querySeq, record.seq):
          if not posPair[0] == posPair[1]:
            sumDist += 1
        #sumDist += Levenshtein.distance(str(querySeq), str(record.seq))
    return(sumDist)

def trimEndGaps(alignment, trimThreshold=0.5):
    # trim leading and trailing gaps from alignment
    # trim columns where proportion of gaps is > trimThreshold, stop at first column below threshold
    # higher threshold means more trimming
    # return tuple with number of columns trimmed front and back
    if trimThreshold < 0 or trimThreshold >= 1.0:
        raise Exception("trimEndGaps: trimThreshold (%.2f) out of range 0.0-1.0"%trimThreshold)
    trimmableNumberOfGaps = int(trimThreshold * len(alignment)) # number of sequences with gaps allowed
    if trimmableNumberOfGaps == 0:
        if Debug:
            LOG.write("trimEndGaps: trimThreshold (%.2f) so lenient no trimming possible with %d sequences\n"%(trimThreshold, len(alignment)))
        #alignment.annotation["endgaps_trimmed"] = (0,0)
        return(alignment)
    leadGaps={}
    endGaps={}
    for record in alignment:
        leadGapLen = 0
        m = re.match("-+", str(record.seq))
        if m:
            leadGapLen=len(m.group(0))
        if leadGapLen not in leadGaps:
            leadGaps[leadGapLen] = 0
        leadGaps[leadGapLen] += 1
        endGapLen = 0
        m = re.search("-+$", str(record.seq))
        if m:
            endGapLen=len(m.group(0))
        if endGapLen not in endGaps:
            endGaps[endGapLen] = 0
        endGaps[endGapLen] += 1
    leadTrimLen = 0
    gappedSeqs = len(alignment) # start considering all sequence as gaps (at pos -1)
    for l in sorted(leadGaps):
        gappedSeqs -= leadGaps[l]
        if gappedSeqs <= trimmableNumberOfGaps:
            leadTrimLen = l
            break
    gappedSeqs = len(alignment) # start considering all sequence as gaps (at pos len+1)
    endTrimLen=0
    for l in sorted(endGaps):
        gappedSeqs -= endGaps[l]
        if gappedSeqs <= trimmableNumberOfGaps:
            endTrimLen = l
            break
    if leadTrimLen > 0 or endTrimLen > 0:
        endPos = alignment.get_alignment_length() - endTrimLen
        trimmedAlignment = alignment[:,leadTrimLen:endPos]
        for i in range(0, len(alignment)):
            trimmedAlignment[i].annotations = alignment[i].annotations
        alignment = trimmedAlignment
    #alignment.annotation["endgaps_trimmed"] = (leadTrimLen, endTrimLen)
    return(alignment)

def resolveDuplicatesPerPatricGenome(alignment):
# accepts Bio.MultipleSeqAlignment, returns list of seqIds to keep
# calculate average similarity/distance of each seqToResolve to entire alignment
# identify best seq per genome (highest avg similarity to rest of alignment) and save that one, remove others from seqIds list
# return list of seqIds to keep
    seqIds=list()
    genomesToResolve=set()
    seqsPerGenome = {}
    initialNumRecords=len(alignment)
    for record in alignment:
        seqIds.append(record.id)
        # assume PATRIC gene identifier like this: fig|1399771.3.peg.1094
        genomeId = ".".join(record.id.split(".")[:2]).split("|")[1]
        #LOG.write("%s\t%s\n"%(record.id, genomeId))
        if genomeId not in seqsPerGenome:
            seqsPerGenome[genomeId] = []
        else:
            genomesToResolve.add(genomeId)
        seqsPerGenome[genomeId].append(record.id)

    for genomeId in genomesToResolve:
        setToResolve = seqsPerGenome[genomeId]
        bestAlignDist = alignment.get_alignment_length()
        seqToKeep = None
        for seqId in setToResolve:
            dist = calcSumAlignmentDistance(alignment, record.seq)
            if seqToKeep == None:
                bestAlignDist=dist
                seqToKeep=seqId
            else:
                if dist < bestAlignDist:
                    bestAlignDist = dist
                    seqToKeep = seqId
        for seqId in setToResolve:
            if not seqId == seqToKeep:
                seqIds.remove(seqId)
    if len(seqIds) < initialNumRecords:
        if Debug:
            LOG.write("after resolveDups num seqIds is %d, versus prev %d\n"%(len(seqIds), initialNumRecords))
        reducedSeqs = []
        for record in alignment:
            if record.id in seqIds:
                record.seq = Seq(str(record.seq).replace('-', ''), record.seq.alphabet) # remove gaps
                reducedSeqs.append(record)
        alignment = alignSeqRecordsMuscle(reducedSeqs)
    return alignment

def proteinToCodonAlignment(proteinAlignment, extraDnaSeqs = None):
    protSeqDict = {}
    for seqRecord in proteinAlignment:
        protSeqDict[seqRecord.id] = seqRecord
    dnaFasta = patric_api.getDnaFastaForPatricIds(protSeqDict.keys())
    #if Debug:
    #     LOG.write("dnaFasta sample: %s\n"%dnaFasta[:100])

    dnaSeqDict = SeqIO.to_dict(SeqIO.parse(StringIO.StringIO(dnaFasta), "fasta", alphabet=IUPAC.ambiguous_dna))
    for seqId in protSeqDict:
        if extraDnaSeqs and seqId in extraDnaSeqs:
            dnaSeqDict[seqId] = extraDnaSeqs[seqId]
            if Debug:
                LOG.write("appending extra DNA seq %s\n"%seqId)
    if set(dnaSeqDict.keys()) != set(protSeqDict.keys()):
        raise Exception("Protein and DNA sets differ:\nProteins: %s\nDNA: %s\n"%(", ".join(sorted(protSeqDict)), ", ".join(sorted(dnaSeqDict))))
    allGood = True
    for seqId in dnaSeqDict:
        if not len(dnaSeqDict[seqId].seq):
            allGood = False
            #del(dnaSeqDict[seqId])
            LOG.write("warning: seqId %s length of dna was zero\n"%seqId)
    dnaSeqRecords=[]
    for proteinSeq in proteinAlignment:
        dnaSeqRecords.append(dnaSeqDict[proteinSeq.id])

    if Debug:
        LOG.write("dna seqs has %d seqs\n"%(len(dnaSeqRecords)))
        #LOG.write("DNA seq ids: %s\n"%(", ".join(sorted(dnaSeqDict))))
        #LOG.write("pro seq ids: %s\n"%(", ".join(sorted(protSeqDict))))
        #LOG.write("first two aligned DNA seqs:\n")
        #SeqIO.write(dnaSeqRecords[:2], LOG, "fasta")
        #LOG.flush()
   
    """
    # now check length of protein vs dna sequences, extend dna if needed to make match in numbers of codons
    for i, protRec in enumerate(proteinAlignment):
        protSeq = str(protRec.seq)
        protSeq.replace('-','')
        protLen = len(protSeq)
        if len(dnaSeqs[i].seq) < protLen*3:
            shortfall = (protLen*3) - len(dnaSeqs[i].seq)
            if Debug:
                LOG.write("DNA seq for %s is too short for protein, shortfall = %d\n"%(protRec.id, shortfall))
            # extend on both ends to be safe
            dnaSeqs[i].seq = "N"*shortfall + dnaSeqs[i].seq + "N"*shortfall
    """
    returnValue = None
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', BiopythonWarning)
        try:
            returnValue = codonalign.build(proteinAlignment, dnaSeqRecords, max_score=1000)
            for dnaSeq in returnValue:
                proteinRecord = protSeqDict[dnaSeq.id]
                if proteinRecord.annotations:
                    dnaSeq.annotations = proteinRecord.annotations.copy()

        except Exception as e:
            LOG.write("problem in codonalign, skipping\n%s\n"%str(e))
            #raise
    return returnValue
    
def relabelSequencesByGenomeId(seqRecordSet):
    #rename sequences by genome instead of sequence Id
    for seqRecord in seqRecordSet:
        originalId = seqRecord.id
        genomeId = ".".join(seqRecord.id.split(".")[:2]).split("|")[1]
        seqRecord.id = genomeId
        #stash original ID in the annotations dictionary of the seqRecord
        if not seqRecord.annotations:
            seqRecord.annotations={}
        seqRecord.annotations['original_id'] = originalId

def trimAlignments(alignmentDict, endGapTrimThreshold=0.5):
    if endGapTrimThreshold == 0:
        return
    for alignmentId in alignmentDict:
        alignmentDict[alignmentId] = trimEndGaps(alignmentDict[alignmentId], endGapTrimThreshold)
    return

def writeOneAlignmentPhylip(alignment, destination, idList, outputIds=True):
    nameFieldWidth = 10
    if outputIds:
        for id in idList:
            nameFieldWidth = max(nameFieldWidth, len(id))
    theseIds = set(rec.id for rec in alignment)
    if len(theseIds) < len(idList):
        #nullSeq = UnknownSeq(alignment.get_alignment_length(), alphabet=alignment._alphabet)
        nullSeq = '-'*alignment.get_alignment_length()
    seqDict = {}
    for record in alignment:
        seqDict[record.id] = record.seq
    for id in idList:
        if outputIds:
            destination.write("{:{width}}  ".format(id, width = nameFieldWidth))
        if id in theseIds:
            destination.write(str(seqDict[id])+"\n")
        else:
            destination.write(str(nullSeq)+"\n")

def writeConcatenatedAlignmentsPhylip(alignments, destination):
# alignments is dictionary of Bio.multipleSeqAlignments
    if type(destination) == str:
        destination = open(destination, "w")
    taxonIdSet = set(seq.id for aln in alignments for seq in alignments[aln])
    totalLength = 0
    for alignmentId in alignments:
        totalLength += alignments[alignmentId].get_alignment_length()
    destination.write("%d\t%d\n"%(len(taxonIdSet), totalLength))
    taxonIdList=sorted(taxonIdSet)
    alignmentIdList = sorted(alignments)
    writeOneAlignmentPhylip(alignments[alignmentIdList[0]], destination, taxonIdList, outputIds=True)
    for alignmentId in alignmentIdList[1:]:
        destination.write("\n")
        writeOneAlignmentPhylip(alignments[alignmentId], destination, taxonIdList, outputIds=False)

def outputCodonsProteinsPhylip(codonAlignments, proteinAlignments, destination):
    if Debug:
        LOG.write("outputCodonsProteinsPhylip, ncodonAln=%d, nprotAln=%d, destType=%s\n"%(len(codonAlignments), len(proteinAlignments), type(destination)))
    if len(codonAlignments) == 0:
        LOG.write("outputCodonsProteinsPhylip() called with zero codonAlignments()\n")
        return
    if type(destination) == str or type(destination) == unicode:
        if Debug:
            LOG.write("outputCodonsProteinsPhylip opening file %s\n"%destination)
        destination = open(destination, "w")
    codonPositions = 0
    taxonSet=set()
    for alignmentId in codonAlignments:
        codonPositions += codonAlignments[alignmentId].get_alignment_length()
        for seqRecord in codonAlignments[alignmentId]:
            taxonSet.add(seqRecord.id)
    proteinPositions = 0
    for alignmentId in proteinAlignments:
        proteinPositions += proteinAlignments[alignmentId].get_alignment_length()
        for seqRecord in proteinAlignments[alignmentId]:
            taxonSet.add(seqRecord.id)
    taxonIdList = sorted(taxonSet)
    #destination = open(directory+fileBase+"+codonsAndProteins.phy", 'w')
    destination.write("%d\t%d\n"%(len(taxonIdList), codonPositions+proteinPositions))
    alignmentIdList = sorted(codonAlignments)
    writeOneAlignmentPhylip(codonAlignments[alignmentIdList[0]], destination, taxonIdList, outputIds=True)
    for alignmentId in alignmentIdList[1:]:
        destination.write("\n")
        writeOneAlignmentPhylip(codonAlignments[alignmentId], destination, taxonIdList, outputIds=False)
    for alignmentId in sorted(proteinAlignments):
        destination.write("\n")
        writeOneAlignmentPhylip(proteinAlignments[alignmentId], destination, taxonIdList, outputIds=False)
    destination.write("\n")
    destination.close()

    return()

