from flask import Flask, jsonify, abort, request, jsonify
from pyproteinsExt import uniprot as pExt
from . import collectionProxy as redisCollection
from decorator import decorator
import re
from functools import wraps





UNIPROT_COLLECTION=None
# temp hack
REDIS=False

def cleanup(rh=None, rp=None):
    redisCollection.bootstrap(host=rh, port=rp)
    redisCollection.cleanup()

def startup(xmlUniprot, redis=False, rh=None, rp=None):

    global UNIPROT_COLLECTION, REDIS
   # UNIPROT_COLLECTION = pExt.EntrySet(collectionXML=xmlUniprot)
   
    if xmlUniprot:
        print(f"Loading XML ressource {xmlUniprot} ...")
        _ = pExt.EntrySet(collectionXML=xmlUniprot)
        UNIPROT_COLLECTION = _    
    if redis:
        REDIS=True
        redisCollection.bootstrap(host=rh, port=rp)        
        if xmlUniprot:
            redisCollection.convert(_)
        UNIPROT_COLLECTION = redisCollection
        
    print("Uniprot storage service listening")
 
    app = Flask(__name__)

    app.add_url_rule('/ping', 'ping', ping)
    
    app.add_url_rule( '/api/entry/<uniprotID>', 'getProtein', getProtein,
                      methods = ['GET'] ) 
    
    app.add_url_rule( '/api/entries', view_func=getProteins,
                      methods = ['POST'] )

    app.add_url_rule( '/api/length', 'length', length,
                      methods = ['GET'] ) 
    
    app.add_url_rule( '/api/list/<interval>', 'listProtein', listProtein)

    app.add_url_rule( '/api/goterms',view_func=getGoTerms,
                      methods = ['POST'] )
    #app.add_url_rule('/unigo/<taxid>', 'view_unigo', view_unigo)
    
    return app

def unwrapUniprotList():
    if not request.is_json:
        print("Post data not json")
        abort(422)
    data = request.get_json()
    if not 'uniprotIDs' in data:
        print("Post data lacks uniprotIDs key")
        abort(422)
    _ = {k:v for k,v in data.items() if k != 'uniprotIDs'} 
   
    return (data['uniprotIDs'], _)

def getGoTerms():
    (uniprotList, opt) = unwrapUniprotList()
    
    entryIter = UNIPROT_COLLECTION.mget(uniprotList, raw=False)

    GoBag = { }

    for e in entryIter:
        if not e is None:
            print(type(e))
            for go in e.GO:
                print(go)
                if not go.id in GoBag:
                    GoBag[go.id] = go.asDict
                    GoBag[go.id]["members"] = []
                GoBag[go.id]["members"].append(e.id)
    
    return jsonify(GoBag)

def listProtein(interval=':1000'): 
    def parseInterval(string): 
        """Parse a slice-like expression"""       
        m1 = re.match("^:{0,1}([\d]+)$", string)
        if m1:           
            return (0, int(m1[1]))
        m2 = re.match("^([\d]+):$", string)
        if m2:
            return (int(m2[1]), None)
        m3 = re.match("^([\d]+):([\d]+)$", string)
        if m3:
            return (int(m3[1]), int(m3[2]))
    
    cstart, cstop = parseInterval(interval)
    #print(cstart, cstop)
    listIDs = redisCollection.getSliceIDs(cstart, cstop)
    
    return jsonify( {"entryIDs" : listIDs} )

def length():
    global UNIPROT_COLLECTION
    if REDIS:
        length = UNIPROT_COLLECTION.length()
    else:
        length = len(UNIPROT_COLLECTION)
    print(f"Current Collection size ${length}")
    return jsonify( {"totalEntry" : length} )

def getProtein(uniprotID):
    print(f"Seeking {uniprotID}")
    oProtein = UNIPROT_COLLECTION.get(uniprotID)
    if not oProtein:
        abort(404)
    return jsonify( { uniprotID : oProtein.toJSON() } )

def getProteins():
    uniprotList, opt = unwrapUniprotList()
    ## implements redis
    results = {}
    validCnt = 0
    # tmp hack  
    if not REDIS:
        for id in uniprotList:
            _ = UNIPROT_COLLECTION.get(id)
            results[id] = _.toJSON() if not _ is None else None
            validCnt = validCnt + 1 if results[id] else validCnt
    else:
        print("Fetching many via redis")
        for e in UNIPROT_COLLECTION.mget(uniprotList, raw=False):
            results[e.id] = e.toJSON() if not e is None else None
            validCnt = validCnt + 1 if results[e.id] else validCnt

    print(f"Returning { validCnt } valid elements")
    return jsonify(results)



def ping():
    #pyproteinsExt.model()
    return "pong"