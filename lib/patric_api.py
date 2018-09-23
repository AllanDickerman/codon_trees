import os
import sys
import re
import requests
import urllib
from Bio.Alphabet import IUPAC
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import json
from requests.packages.urllib3.exceptions import InsecureRequestWarning
requests.packages.urllib3.disable_warnings(InsecureRequestWarning)

Debug = False #shared across functions defined here
LOG = sys.stderr
Base_url="https://www.patricbrc.org/api/"

Session = requests.Session()
UserAtPatric = None
if os.environ.has_key("KB_AUTH_TOKEN"):
    Session.headers.update({ 'Authorization' : os.environ.get('KB_AUTH_TOKEN') })
elif os.path.exists(os.path.join(os.environ.get('HOME'), ".patric_token")):
    F = open(os.path.join(os.environ.get('HOME'), ".patric_token"))
    Session.headers.update({ 'Authorization' : F.read() })
    F.close()
if "authorization" in Session.headers:
    UserAtPatric = Session.headers["Authorization"].split(r"|")[3].split("=")[1]
    LOG.write("Patric user = %s\n"%UserAtPatric)
Session.headers.update({ 'accept': "text/tsv" })
Session.headers.update({ "Content-Type": "application/rqlquery+x-www-form-urlencoded" })
LOG.write(str(Session.headers)+"\n")

def getGenomeIdsNamesByName(name, limit='10'):
    query = "eq(genome_name,%s)"%name
    query += "&select(genome_id,genome_name)"
    query += "&limit(%s)"%limit
    ret = Session.get(Base_url+"genome/", params=query)
    if Debug:
        LOG.write(ret.url+"\n")
    return(ret.text.replace('"', ''))

def getGenomeGroupIds(genomeGroupName):
    genomeGroupSpecifier = UserAtPatric+"/home/Genome Groups/"+genomeGroupName
    genomeGroupSpecifier = "/"+urllib.quote(genomeGroupSpecifier)
    genomeGroupSpecifier = genomeGroupSpecifier.replace("/", "%2f")
    query = "in(genome_id,GenomeGroup("+genomeGroupSpecifier+"))"
    query += "&select(genome_id)"
    query += "&limit(1000)"
    if Debug:
        LOG.write("requesting group %s for user %s\n"%(genomeGroupName, UserAtPatric))
        LOG.write("query =  %s\n"%(query))
    ret = Session.get(Base_url+"genome/", params=query)
    if Debug:
        LOG.write(ret.url+"\n")
    return(ret.text.replace('"', '').split("\n"))[1:-1]

def getNamesForGenomeIds(genomeIdList):
    return getDataForGenomes(genomeIdList, ["genome_id", "genome_name"])

def getGenomeIdsByFieldValue(queryField, queryValue):
    req = sesssion.get(Base_url+"genome/", params="in(%s,%s)"%(queryField, queryValue)) 
    retval = []
    if Debug:
        LOG.write("getGenomeIdsByQuery: "+req.url+"\n")
        LOG.write(req.text+"\n")
    for line in req.text.split("\n"):
       retval.append(line)
    return retval

def getDataForGenomes(genomeIdList, fieldNames):
    query = "in(genome_id,(%s))"%",".join(genomeIdList)
    if fieldNames:
        query += "&select(%s)"%",".join(fieldNames)
    query += "&limit(%s)"%len(genomeIdList)

    response = Session.get(Base_url+"genome/", params=query)
    if not response.ok:
        LOG.write("Error code %d returned by %s in getDataForGenomes\n"%(response.status_code, req.url))
        LOG.write("length of query was %d\n"%len(query))
        LOG.write("url="+req.url+"\nquery="+query+"\n")
        raise Exception(errorMessage)
    data = response.text.replace('"','') #get rid of quotes
    rows = data.split("\n")[:-1] # leave off empty last element
    retval = []
    for row in rows:
        fields = row.split("\t")
        #if len(fields) != len(fieldNames):
         #   continue
        retval.append(fields)
    return(retval)

def getGenomeFeaturesByPatricIds(patricIdList, fieldNames=None):
    query="in(patric_id,("+",".join(map(urllib.quote, patricIdList))+"))"
    if fieldNames:
        query += "&select(%s)"%",".join(fieldNames)
    query += "&limit(%d)"%len(patricIdList)
    response=Session.get(Base_url+"genome_feature/", params=query)
    if not response.ok:
        LOG.write("Error code %d returned by %s in getGenomeFeaturesByPatricIds\n"%(response.status_code, Base_url))
        LOG.write("length of query was %d\n"%len(query))
        LOG.write("query="+req.url+"\n")
        errorMessage= "Error code %d returned by %s in getGenomeFeaturesByPatricIds\nlength of query was %d\n"%(response.status_code, Base_url, len(query))
        raise Exception(errorMessage)
    data = response.text.replace('"','')
    rows = data.split("\n")
    retval = []
    for row in rows[:-1]: # last line is empty (because of terminal line return)
        fields = row.split("\t")
        if len(fields) != len(fieldNames):
            LOG.write("getGenomeFeaturesByPatricIds: parsed fields (%d) is fewer than requested fields (%d):\n%s\n"%(len(fields), len(fieldNames), row))
            continue
        retval.append(fields)
    return(retval)

def getProteinBioSeqRecordsForPatricIds(patricIdList):
    data = getGenomeFeaturesByPatricIds(patricIdList, ["patric_id", "product", "genome_id", "aa_sequence"])
    #if Debug:
        #LOG.write("getProteinBioSeqRecordsForPatricIds, about to parse:\n%s\n"%data)
        #LOG.write("num rows = %d\n"%len(rows))
    retval = [] # will be list of Bio.SeqRecord
    for row in data[1:]: # first row is headers
        if len(row) < 4:
            LOG.write("problem in getProteinBioSeqRecordsForPatricIds, expected 4 fields and got %d: %s\n"%(len(row), "|".join(row)))
            continue
        patricId, product, genomeId, aa_sequence = row
        simpleSeq = Seq(aa_sequence, IUPAC.extended_protein)
        seqRecord = SeqRecord(simpleSeq, id=patricId, description=product)
        seqRecord.annotations["genome_id"] = genomeId
        retval.append(seqRecord)
    return(retval) # list of Bio.SeqRecord

def getDnaBioSeqRecordsForPatricIds(patricIdList):
    data = getGenomeFeaturesByPatricIds(patricIdList, ["patric_id", "product", "genome_id", "na_sequence"])
    #if Debug:
     #   LOG.write("getDnaBioSeqRecordsForPatricIds, first two rows of data:\n")
     #   LOG.write(rows[0]+"\n")
     #   LOG.write(rows[1]+"\n")
    retval = []
    for fields in data[1:]: #first row is headers
        if len(fields) < 4:
            LOG.write("problem in getDnaBioSeqRecordsForPatricIds, expected 4 fields and got %d: %s\n"%(len(fields), "|".join(fields)))
            continue
        patricId, product, genomeId, na_sequence = fields[:4]
        simpleSeq = Seq(na_sequence, IUPAC.ambiguous_dna)
        seqRecord = SeqRecord(simpleSeq, id=patricId, description=product)
        seqRecord.annotations["genome_id"] = genomeId
        retval.append(seqRecord)
    return(retval)

def getAllProteinsForGenomeId(genomeId):
    query="in(genome_id,("+genomeId+"))"
    fieldNames = ["patric_id", "product", "genome_name", "genome_id", "aa_sequence"]
    query += "&select(%s)"%",".join(fieldNames)
    query += "&limit(%d)"%10000
    """
    req = requests.Request('POST', Base_url+"genome_feature/", data=query)
    prepared = Session.prepare_request(req) #req.prepare()
    response=Session.send(prepared, verify=False)
    """
    response = Session.get(Base_url+"genome_feature/", params=query) #, 
    if not response.ok:
        LOG.write("Error code %d returned by %s in getAllProteinsForGenomeId\nlength of query was %d\n"%(response.status_code, Base_url, len(query)))
        LOG.write("query="+req.url+"\n")
        errorMessage= "Error code %d returned by %s\nlength of query was %d\n"%(response.status_code, Base_url, len(query))
        raise Exception(errorMessage)
    rows = response.text.split("\n")
    retval = []
    for row in rows[:-1]: # last line is empty (because of terminal line return)
        fields = row.replace('"', '').split("\t")
        if len(fields) != len(fieldNames):
            LOG.write("getGenomeFeaturesByPatricIds: parsed fields (%d) is fewer than requested fields (%d):\n%s\n"%(len(fields), len(fieldNames), row))
        retval.append(fields)
    return(retval)


def getPatricGenesPgfamsForGenomeList(genomeIdList):
    if Debug:
        LOG.write("getPatricGenesPgfamsForGenomeList() called for %d genomes\n"%len(genomeIdList))
        LOG.write("    Session headers=\n"+str(Session.headers)+"\n")
    retval = []
    # one genome at a time, so using 'get' should be fine
    for genomeId in genomeIdList:
        query="and(%s,%s,%s)"%("eq(genome_id,(%s))"%genomeId, "eq(feature_type,CDS)", "eq(pgfam_id,PGF*)")
        query += "&select(genome_id,patric_id,pgfam_id)"
        query += "&limit(10000)"
        response = Session.get(Base_url+"genome_feature/", params=query) #, 
        """
        req = requests.Request('POST', Base_url+"genome_feature/", data=query)
        prepared = Session.prepare_request(req) #req.prepare()
        response=Session.send(prepared, verify=False)
        """
        if Debug:
            LOG.write("    response URL: %s\n"%response.url)
            LOG.write("    len(response.text)= %d\n"%len(response.text))
        curLen = len(retval)
        for line in response.text.split("\n"):
            line = line.replace('"','')
            row = line.split("\t")
            if len(row) != 3:
                continue
            if not row[2].startswith("PGF"):
                continue
            retval.append(row)
        if Debug:
            LOG.write("    got %d pgfams for that genome\n"%(len(retval)-curLen))
    return(retval)

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


