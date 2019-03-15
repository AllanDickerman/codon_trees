import os
import sys
import re
import requests
import urllib
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
    Session.headers.update({ 'Authorization' : F.read().rstrip() })
    F.close()
if "authorization" in Session.headers:
    UserAtPatric = Session.headers["Authorization"].split(r"|")[3].split("=")[1]
    LOG.write("Patric user = %s\n"%UserAtPatric)
Session.headers.update({ 'accept': "text/tsv" })
Session.headers.update({ "Content-Type": "application/rqlquery+x-www-form-urlencoded" })
if Debug:
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

def getNamesForGenomeIds(genomeIdSet):
    return getDataForGenomes(genomeIdSet, ["genome_id", "genome_name"])

def getGenomeIdsByFieldValue(queryField, queryValue):
    req = sesssion.get(Base_url+"genome/", params="in(%s,%s)"%(queryField, queryValue)) 
    retval = []
    if Debug:
        LOG.write("getGenomeIdsByQuery: "+req.url+"\n")
        LOG.write(req.text+"\n")
    for line in req.text.split("\n"):
       retval.append(line)
    return retval

def getDataForGenomes(genomeIdSet, fieldNames):
    query = "in(genome_id,(%s))"%",".join(genomeIdSet)
    if fieldNames:
        query += "&select(%s)"%",".join(fieldNames)
    query += "&limit(%s)"%len(genomeIdSet)

    response = Session.get(Base_url+"genome/", params=query)
    if Debug:
        LOG.write("getDataForGenomes:\nurl="+response.url+"\nquery="+query+"\n")
    if not response.ok:
        LOG.write("Error code %d returned by %s in getDataForGenomes\n"%(response.status_code, response.url))
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
"""
# commented because this will no longer work with solr schema change, need to select my md5
def getGenomeFeaturesByPatricIds(patricIdList, fieldNames=None):
    query="in(patric_id,("+",".join(map(urllib.quote, patricIdList))+"))"
    if fieldNames:
        query += "&select(%s)"%",".join(fieldNames)
    query += "&limit(%d)"%len(patricIdList)
    response=Session.get(Base_url+"genome_feature/", params=query)
    if Debug:
        LOG.write("getGenomeFeaturesByPatricIds:\nurl="+response.url+"\nquery="+query+"\n")
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
"""
def getProteinFastaForPatricIds(patricIds):
    query="in(patric_id,("+",".join(map(urllib.quote, patricIds))+"))"
    query += "&limit(%d)"%len(patricIds)
    response=Session.get(Base_url+"genome_feature/", params=query, headers={'Accept': 'application/protein+fasta'})
    if Debug:
        LOG.write("getProteinFastaForByPatricIds:\nurl="+response.url+"\nquery="+query+"\n")
    if not response.ok:
        LOG.write("Error code %d returned by %s in getProteinFastaForPatricIds\n"%(response.status_code, Base_url))
        errorMessage= "Error code %d returned by %s in getGenomeFeaturesByPatricIds\nlength of query was %d\n"%(response.status_code, Base_url, len(query))
        LOG.write(errorMessage)
        LOG.flush()
        raise Exception(errorMessage)
    idsFixedFasta=""
    for line in response.text.split("\n"):
        if line.startswith(">"):
            parts = line.split("|")
            if len(parts) > 2:
                line = "|".join(parts[:2])
        idsFixedFasta += line+"\n"
    return idsFixedFasta
    
def getDnaFastaForPatricIds(patricIds):
    query="in(patric_id,("+",".join(map(urllib.quote, patricIds))+"))"
    query += "&limit(%d)"%len(patricIds)
    response=Session.get(Base_url+"genome_feature/", params=query, headers={'Accept': 'application/dna+fasta'})
    if Debug:
        LOG.write("getDnaFastaForByPatricIds:\nurl="+response.url+"\nquery="+query+"\n")
    if not response.ok:
        LOG.write("Error code %d returned by %s in getDnaFastaForPatricIds\n"%(response.status_code, Base_url))
        errorMessage= "Error code %d returned by %s in getGenomeFeaturesByPatricIds\nlength of query was %d\n"%(response.status_code, Base_url, len(query))
        LOG.write(errorMessage)
        LOG.flush()
        raise Exception(errorMessage)
    idsFixedFasta=""
    for line in response.text.split("\n"):
        if line.startswith(">"):
            parts = line.split("|")
            if len(parts) > 2:
                line = "|".join(parts[:2])
        idsFixedFasta += line+"\n"
    return idsFixedFasta
    
def getProteinsFastaForGenomeId(genomeId):
    query="in(genome_id,("+genomeId+"))"
    query += "&limit(25000)"
    response=Session.get(Base_url+"genome_feature/", params=query, headers={'Accept': 'application/protein+fasta'})
    if Debug:
        LOG.write("getProteinsFastaForGenomeId:\nurl="+response.url+"\nquery="+query+"\n")
    if not response.ok:
        LOG.write("Error code %d returned by %s in getProteinsFastaForGenomeId\n"%(response.status_code, Base_url))
        errorMessage= "Error code %d returned by %s in getProteinsFastaForGenomeId\nlength of query was %d\n"%(response.status_code, Base_url, len(query))
        LOG.write(errorMessage)
        LOG.flush()
        raise Exception(errorMessage)
    idsFixedFasta=""
    for line in response.text.split("\n"):
        if line.startswith(">"):
            parts = line.split("|")
            if len(parts) > 2:
                line = "|".join(parts[:2])+"\n"
        idsFixedFasta += line
    return idsFixedFasta

def getPatricGenesPgfamsForGenomeSet(genomeIdSet):
    if Debug:
        LOG.write("getPatricGenesPgfamsForGenomeSet() called for %d genomes\n"%len(genomeIdSet))
        LOG.write("    Session headers=\n"+str(Session.headers)+"\n")
    retval = []
    # one genome at a time, so using 'get' should be fine
    for genomeId in genomeIdSet:
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
