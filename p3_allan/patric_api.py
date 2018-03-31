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

debug = False #shared across functions defined here
base_url="https://www.patricbrc.org/api/"

session = requests.Session()
if os.environ.has_key("KB_AUTH_TOKEN"):
    session.headers.update({ 'Authorization' : os.environ.get('KB_AUTH_TOKEN') })
    print(session.headers)


def getGenomeIdsNamesByName(name, limit='10'):
    query = "eq(genome_name,%s)"%name
    query += "&select(genome_id,genome_name)"
    query += "&limit(%s)"%limit
    ret = session.get(base_url+"genome/", params=query, headers={"accept":"text/tsv"})
    if debug:
        sys.stderr.write(ret.url+"\n")
    return(ret.text.replace('"', ''))

def getGenomeGroupIds(genomeGroupName):
    if not authorization:
        raise Exception("call to getGenomeGroupIds without an authorization string")
    userAtPatric = authorization.split(r"|")[3].split("=")[1]
    genomeGroupSpecifier = userAtPatric+"/home/Genome Groups/"+genomeGroupName
    genomeGroupSpecifier = "/"+urllib.quote(genomeGroupSpecifier)
    genomeGroupSpecifier = genomeGroupSpecifier.replace("/", "%2f")
    query = "in(genome_id,GenomeGroup("+genomeGroupSpecifier+"))"
    query += "&select(genome_id)"
    query += "&limit(1000)"
    #query = urllib.quote(query)
    headers = {"accept":"text/tsv"}
    headers["Authorization"] = authorization
    if debug:
        sys.stderr.write("requesting group %s for user %s\n"%(genomeGroupName, userAtPatric))
        sys.stderr.write("query =  %s\n"%(query))
    ret = requests.get(base_url+"genome/", params=query, headers=headers)
    if debug:
        sys.stderr.write(ret.url+"\n")
    return(ret.text.replace('"', '').split("\n"))[1:-1]

def getNamesForGenomeIds(genomeIdList):
    return getDataForGenomes(genomeIdList, ["genome_id", "genome_name"])

def getDataForGenomes(genomeIdList, fieldNames):
    query = "in(genome_id,(%s))"%",".join(genomeIdList)
    if fieldNames:
        query += "&select(%s)"%",".join(fieldNames)
    query += "&limit(%s)"%len(genomeIdList)

    headers={"Content-Type": "application/rqlquery+x-www-form-urlencoded", "accept":"text/tsv"}
    req = session.Request('POST', base_url+"genome/", headers=headers, data=query)
    prepared = req.prepare()
    #pretty_print_POST(prepared)
    s = requests.Session()
    response=s.send(prepared, verify=False)
    if not response.ok:
        sys.stderr.write("Error code %d returned by %s in getGenomeFeaturesByPatricIds\nlength of query was %d\n"%(response.status_code, req.url, len(query)))
        sys.stderr.write("url="+req.url+"\nquery="+query+"\n")
        errorMessage= "Error code %d returned by %s in getGenomeFeaturesByPatricIds\nlength of query was %d\n"%(response.status_code, req.url, len(query))
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
    #session.get(base_url+"genome", params=query, headers={"accept":"text/tsv"})
    #return(ret.text.replace('"', ''))

def getGenomeFeaturesByPatricIds(patricIdList, fieldNames=None):
    query="in(patric_id,("+",".join(map(urllib.quote, patricIdList))+"))"
    if fieldNames:
        query += "&select(%s)"%",".join(fieldNames)
    query += "&limit(%d)"%len(patricIdList)
    headers={"Content-Type": "application/rqlquery+x-www-form-urlencoded", "accept":"text/tsv"}
    req = requests.Request('POST', base_url+"genome_feature/", headers=headers, data=query)
    prepared = req.prepare()
    req = requests.Request('POST', base_url+"genome_feature/", headers=headers, data=query)
    prepared = req.prepare()
    #pretty_print_POST(prepared)

    response=session.send(prepared, verify=False)
    if not response.ok:
        sys.stderr.write("Error code %d returned by %s in getGenomeFeaturesByPatricIds\nlength of query was %d\n"%(response.status_code, base_url, len(query)))
        sys.stderr.write("query="+req.url+"\n")
        errorMessage= "Error code %d returned by %s in getGenomeFeaturesByPatricIds\nlength of query was %d\n"%(response.status_code, base_url, len(query))
        raise Exception(errorMessage)
    if not response.ok:
        sys.stderr.write("Error code %d returned by %s in getGenomeFeaturesByPatricIds\nlength of query was %d\n"%(response.status_code, base_url, len(query)))
        sys.stderr.write("query="+req.url+"\n")
        errorMessage= "Error code %d returned by %s in getGenomeFeaturesByPatricIds\nlength of query was %d\n"%(response.status_code, base_url, len(query))
        raise Exception(errorMessage)
    data = response.text.replace('"','')
    rows = data.split("\n")
    retval = []
    for row in rows[:-1]: # last line is empty (because of terminal line return)
        fields = row.split("\t")
        if len(fields) != len(fieldNames):
            sys.stderr.write("getGenomeFeaturesByPatricIds: parsed fields (%d) is fewer than requested fields (%d):\n%s\n"%(len(fields), len(fieldNames), row))
            continue
        retval.append(fields)
    return(retval)

def getProteinBioSeqRecordsForPatricIds(patricIdList):
    data = getGenomeFeaturesByPatricIds(patricIdList, ["patric_id", "product", "genome_id", "aa_sequence"])
    #if debug:
        #sys.stderr.write("getProteinBioSeqRecordsForPatricIds, about to parse:\n%s\n"%data)
        #sys.stderr.write("num rows = %d\n"%len(rows))
    retval = [] # will be list of Bio.SeqRecord
    for row in data[1:]: # first row is headers
        if len(row) < 4:
            sys.stderr.write("problem in getProteinBioSeqRecordsForPatricIds, expected 4 fields and got %d: %s\n"%(len(row), "|".join(row)))
            continue
        patricId, product, genomeId, aa_sequence = row
        simpleSeq = Seq(aa_sequence, IUPAC.extended_protein)
        seqRecord = SeqRecord(simpleSeq, id=patricId, description=product)
        seqRecord.annotations["genome_id"] = genomeId
        retval.append(seqRecord)
    return(retval) # list of Bio.SeqRecord

def getDnaBioSeqRecordsForPatricIds(patricIdList):
    data = getGenomeFeaturesByPatricIds(patricIdList, ["patric_id", "product", "genome_id", "na_sequence"])
    #if debug:
     #   sys.stderr.write("getDnaBioSeqRecordsForPatricIds, first two rows of data:\n")
     #   sys.stderr.write(rows[0]+"\n")
     #   sys.stderr.write(rows[1]+"\n")
    retval = []
    for fields in data[1:]: #first row is headers
        if len(fields) < 4:
            sys.stderr.write("problem in getDnaBioSeqRecordsForPatricIds, expected 4 fields and got %d: %s\n"%(len(fields), "|".join(fields)))
            continue
        patricId, product, genomeId, na_sequence = fields[:4]
        simpleSeq = Seq(na_sequence, IUPAC.ambiguous_dna)
        seqRecord = SeqRecord(simpleSeq, id=patricId, description=product)
        seqRecord.annotations["genome_id"] = genomeId
        retval.append(seqRecord)
    return(retval)

def getPatricGenesPgfamsForGenomeList(genomeIdList):
    retval = []
    headers={"accept":"text/tsv"}
    # one genome at a time, so using 'get' should be fine
    for genomeId in genomeIdList:
        query="and(%s,%s,%s)"%("eq(genome_id,(%s))"%genomeId, "eq(feature_type,CDS)", "eq(pgfam_id,PGF*)")
    # one genome at a time, so using 'get' should be fine
    for genomeId in genomeIdList:
        query="and(%s,%s,%s)"%("eq(genome_id,(%s))"%genomeId, "eq(feature_type,CDS)", "eq(pgfam_id,PGF*)")
        query += "&select(genome_id,patric_id,pgfam_id)"
        query += "&limit(10000)"
        if debug:
            sys.stderr.write("getPatricGenesPgfamsForGenomeList: about to request genes and pgfams for genome ids:\n%s\n%s\n"%(base_url, query))
        req = requests.get(base_url+"genome_feature/", params=query) #, 
        if debug:
            sys.stderr.write(req.url+"\n")
        for line in req.text.split("\n"):
            line = line.replace('"','')
            row = line.split(",")
            if len(row) != 3:
                continue
            if not row[2].startswith("PGF"):
                continue
            retval.append(row)
    return(retval)

def getPatricGenesPgfamsForGenomeObject(genomeObject):
# parse a PATRIC genome object (read from json format) for PGFams
    retval = [] # a list of tupples of (genomeId, Pgfam, geneId)
    genomeId = genomeObject['id']
    for feature in genomeObject['features']:
        if 'family_assignments' in feature:
            for famAss in feature['family_assignments']:
                if famAss[0] == 'PGFAM':
                    retval.append((genomeId, feature['id'], famAss[1]))
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


'''44
def getPatricGenesForGenomeList(genomeIdList):
    retval = []
    headers={"accept":"text/tsv"}
    if authorization:
        headers["Authorization"] = authorization
    # one genome at a time, so using 'get' should be fine
    for genomeId in genomeIdList:
        query="and(%s,%s)"%("eq(genome_id,(%s))"%genomeId, "eq(feature_type,CDS)")
        query += "&select(genome_id,patric_id,pgfam_id)"
        query += "&limit(10000)"
        if debug:
            sys.stderr.write("getPatricGenesPgfamsForGenomeList: about to request genes and pgfams for genome ids:\n%s\n%s\n"%(base_url, query))
        req = requests.get(base_url+"genome_feature/", params=query) #, 
        if debug:
            sys.stderr.write(req.url+"\n")
        for line in req.text.split("\n"):
            line = line.replace('"','')
            row = line.split(",")
            if len(row) != 3:
                continue
            if not row[2].startswith("PGF"):
                continue
            retval.append(row)
    return(retval)
'''
