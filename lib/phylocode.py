import sys
import re
import subprocess
import warnings
from Bio import BiopythonExperimentalWarning
warnings.simplefilter('ignore', BiopythonExperimentalWarning) # importing codonalign raises warnings
from Bio.Alphabet import IUPAC
from Bio import AlignIO
from Bio import SeqIO
from Bio import Alphabet
from Bio import codonalign
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import codonalign
from collections import defaultdict
import patric_api
import StringIO

Debug = False #shared across functions defined here
LOG = sys.stderr

def selectSingleCopyPgfams(genomeGenePgfamList, genomeIdList, requiredGenome=None, maxGenomesMissing=0, maxAllowedDups=0):
    # given a list of genome_ids, gene_ids, and pgfam_ids
    # find the set of pgfam_ids which satisfy the minPropPresent and maxAllowedDups criteria for specified genomes
    pgfamGenomeCount={}
    pgfamCount={}
    genomes=set()
    numDups = {}
    if not len(genomeIdList) > 3:
        raise Exception("getSingleCopyPgfams: number of genome IDs too low: %d"%len(genomeIdList))
    minGenomesPresent = len(genomeIdList) - maxGenomesMissing
    if minGenomesPresent <= 0:
        raise Exception("getSingleCopyPgfams: maxGenomesMissing too large: %d"%maxGenomesMissing)
    for row in genomeGenePgfamList:
        genome, gene, pgfam = row
        if genome not in genomeIdList:
            continue # only pay attention to genomes on list
        if pgfam not in pgfamGenomeCount:
            pgfamGenomeCount[pgfam] = {}
            pgfamCount[pgfam] = 0
            numDups[pgfam] = 0
        if genome not in pgfamGenomeCount[pgfam]:
            pgfamGenomeCount[pgfam][genome] = 0
        else:
            numDups[pgfam] += 1
        pgfamGenomeCount[pgfam][genome] += 1
        pgfamCount[pgfam] += 1
        genomes.add(genome)
    suitablePgfamList = []
    for pgfam in sorted(pgfamGenomeCount, key=lambda id: (pgfamGenomeCount[id], -numDups[id])):
        if requiredGenome and not requiredGenome in pgfamGenomeCount[pgfam]:
            continue # does not include a required genome
        if len(pgfamGenomeCount[pgfam]) >= minGenomesPresent and numDups[pgfam] <= maxAllowedDups:
            suitablePgfamList.append(pgfam)
    return suitablePgfamList

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
        outFile.write("\t\"%s\""%labelDict[tax])
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
            outFile.write("\t\t%s \"%s\",\n"%(tax, translate))
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
        if re.search("begin figtree", line):
            inFigtreeBlock = True
        if re.search(r"^end;", line):
            inFigtreeBlock = False
        if inFigtreeBlock:
            m = re.search(r"set\s+(\S.*)=(\S.*);", line)
            if m:
                retval[m.group(1)] = m.group(2)
    return retval

def generateFigtreeImage(nexusFile, outfileName, numTaxa, figtreeJarFile, imageFormat="PDF"):
    if Debug:
        LOG.write("generateTreeFigure(%s, %s, %d, %s, %s)\n"%(nexusFile, outfileName, numTaxa, figtreeJarFile, imageFormat))
    if imageFormat not in ('PDF', 'SVG', 'PNG', 'JPEG'):
        raise Exception("imageFormat %s not in ('PDF', 'SVG', 'PNG', 'JPEG')"%imageFormat)
    figtreeCommand = ['java',  '-jar', figtreeJarFile, '-graphic', imageFormat]
    if numTaxa > 40:
        height = 600 + 15 * (numTaxa - 40) # this is an empirical correction factor to avoid taxon name overlap
        figtreeCommand.extend(['-height', str(int(height))])
    figtreeCommand.extend([nexusFile, outfileName])
    subprocess.call(figtreeCommand)

def checkMuscle():
    subprocess.check_call(['which', 'muscle'])

def alignSeqRecordsMuscle(seqRecords):
    muscleProcess = subprocess.Popen(['muscle', '-quiet'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    SeqIO.write(seqRecords, muscleProcess.stdin, 'fasta')
    muscleProcess.stdin.close()
    alignment = AlignIO.read(muscleProcess.stdout, "fasta", alphabet=seqRecords[0].seq.alphabet)
    alignment.sort()
    return(alignment)

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
    seqIds=[]
    protSeqDict={}
    for seqRecord in proteinAlignment:
        seqIds.append(seqRecord.id)
        protSeqDict[seqRecord.id] = seqRecord
    dnaFasta = patric_api.getDnaFastaForPatricIds(seqIds)
    if Debug:
        LOG.write("dnaFasta sample: %s\n"%dnaFasta[:100])

    dnaSeqRecords = list(SeqIO.parse(StringIO.StringIO(dnaFasta), "fasta", alphabet=IUPAC.ambiguous_dna))
    for seqRecord in proteinAlignment:
        if extraDnaSeqs and seqRecord.id in extraDnaSeqs:
            dnaSeqRecords.append(extraDnaSeqs[seqRecord.id])
    if Debug:
        LOG.write("proteinToCodonAlignment: last DNA seq record before aligning:\n%s"%str(dnaSeqRecords[-1]))
    missingIds=[]
    allGood = True
    for seqRecord in dnaSeqRecords:
        if not len(seqRecord.seq):
            allGood = False
            missingIds.append(seqRecord.id)
    if not allGood:
        tempProtAl = []
        for seqRecord in proteinAlignment:
            if seqRecord.id not in missingIds:
                tempProtAl.append(seqRecord)
        proteinAlignment = tempProtAl
    #force to be in same order (necessary?)
    dnaSeqs_ordered = len(seqIds)*[None] #pre-allocate length of array
    for dnaSeq in dnaSeqRecords:
        index = seqIds.index(dnaSeq.id)
        dnaSeqs_ordered[index] = dnaSeq
        if len(dnaSeq.seq) < 1:
            message = "proteinToCodonAlignment: length of sequence %s is %d\n"%(dnaSeq.id, len(dnaSeq.seq))
            LOG.write(message+"\n")
            LOG.flush()
            raise Exception(message)
    dnaSeqRecords = dnaSeqs_ordered
    if Debug:
        LOG.write("number of dnaseqs missing = %d\n"%len(missingIds))
        LOG.write("dna seqs has %d seqs\n"%(len(dnaSeqRecords)))
        SeqIO.write(dnaSeqRecords[:2], LOG, "fasta")
        LOG.flush()
   
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
    if type(destination) == str:
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

def getPatricGenesPgfamsForGenomeObject(genomeObject):
# parse a PATRIC genome object (read from json format) for PGFams
    retval = [] # a list of tupples of (genomeId, Pgfam, geneId)
    genomeId = genomeObject['id']
    for feature in genomeObject['features']:
        if 'family_assignments' in feature:
            for familyAssignment in feature['family_assignments']:
                if familyAssignment[0] == 'PGFAM':
                    retval.append((genomeId, feature['id'], familyAssignment[1]))
    return retval

def getGenomeObjectProteins(patricIds, genomeObject):
# return dictionary of patricId -> BioPython.SeqRecord
    genomeId = genomeObject['id']
    retval = {}
    for feature in genomeObject['features']:
        patricId, product, genomeId, aa_sequence = '', '', '', ''
        patricId = feature['id']
        if not patricId in patricIds:
            continue
        if "protein_translation" in feature:
            aa_sequence = feature["protein_translation"]
        if 'function' in feature:
            product = feature['function']
        simpleSeq = Seq(aa_sequence, IUPAC.extended_protein)
        seqRecord = SeqRecord(simpleSeq, id=patricId, description=product)
        seqRecord.annotations["genome_id"] = genomeId
        retval[patricId] = seqRecord
    return retval

def getGenomeObjectGeneDna(patricIds, genomeObject):
# return dictionary of patricId -> BioPython.SeqRecord
    genomeId = genomeObject['id']
    contigSeq = {}
    for contig in genomeObject['contigs']:
        contigSeq[contig['id']] = contig['dna']
    retval = {} # dict of SeqRecords
    for feature in genomeObject['features']:
        if not feature['id'] in patricIds:
            continue
        geneId = feature['id']
        if geneId not in patricIds:
            continue
        product = ''
        if 'product' in feature:
            product = feature['function']
        if not 'location' in feature:
            continue
        contig, start, ori, length = feature['location'][0] # this should be an array of (contig, start, orientation, length)
        start = int(float(start))
        length = int(float(length))
        if ori == '+':
            start -= 1
            simpleSeq = Seq(contigSeq[contig][start:start+length], IUPAC.ambiguous_dna)
        if ori == '-':
            simpleSeq = Seq(contigSeq[contig][start-length:start], IUPAC.ambiguous_dna)
            simpleSeq = simpleSeq.reverse_complement()

        seqRecord = SeqRecord(simpleSeq, id=geneId, description=product)
        seqRecord.annotations["genome_id"] = genomeId
        retval[geneId] = seqRecord
    return retval


