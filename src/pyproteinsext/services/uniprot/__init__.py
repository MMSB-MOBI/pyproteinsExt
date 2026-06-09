from flask import Flask, jsonify, abort, request, jsonify
from pyproteinsext import uniprot as pExt
from . import collectionProxy as redisCollection
import re


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
        print("before", len(_))
        _.add('P27573')
        print('after', len(_))
        
        UNIPROT_COLLECTION = _    
    if redis:
        REDIS=True
        redisCollection.bootstrap(host=rh, port=rp)        
        if xmlUniprot:
            redisCollection.convert(_)
        UNIPROT_COLLECTION = redisCollection
        
    print("Uniprot storage service listening")
 
    app = Flask(__name__)

    app.add_url_rule('/model', 'model', model)
    
    app.add_url_rule( '/uniprot/<uniprotID>', 'getProtein', getProtein,
                      methods = ['GET'] ) 
    
    app.add_url_rule( '/uniprots', 'getProteins', getProteins,
                      methods = ['POST'] )

    app.add_url_rule( '/length', 'length', length,
                      methods = ['GET'] ) 
    
    app.add_url_rule( '/uniprot/list', 'list', listProtein)
    app.add_url_rule( '/uniprot/list/<interval>', 'listProtein', listProtein)
    app.add_url_rule('/missing', 'getMissing', get_missing, methods = ['POST'])
    #app.add_url_rule('/unigo/<taxid>', 'view_unigo', view_unigo)
    
    return app

def listProtein(interval=':1000'):
    def parseInterval(string): 
        """Parse a slice-like expression"""       
        m1 = re.match("^([\d]+):{0,1}$", string)
        if m1:
            return (0, int(m1[1]))
        m2 = re.match("^:([\d]+)$", string)
        if m2:
            return (int(m2[1]), None)
        m3 = re.match("^([\d]+):([\d]+)$", string)
        if m3:
            return (int(m3[1]), int(m3[2]))
    
    cstart, cstop = parseInterval(interval)

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

def _load_uniprot_ids_request():
    if not request.is_json:
        print("Post data not json")
        abort(422)
    data = request.get_json()
    if not 'uniprotIDs' in data:
        print("Post data lacks uniprotIDs key")
        abort(422)
    return data


def getProteins():
    data = _load_uniprot_ids_request()
    ## implements redis
    results = {}
    validCnt = 0
    # tmp hack
    if not REDIS:
        for id in data['uniprotIDs']:
            _ = UNIPROT_COLLECTION.get(id)
            results[id] = _.toJSON() if not _ is None else None
            validCnt = validCnt + 1 if results[id] else validCnt
    else:
        print("Fetching many via redis")
        print(data['uniprotIDs'][:10])
        #test = [e.id if e else None for e in UNIPROT_COLLECTION.mget(data['uniprotIDs'], raw=False)]
        #print(test[:10])
        #indices = [i for i, x in enumerate(test) if x == None]
        #print("Missing", indices)

        for uniprot_id, e in zip( data['uniprotIDs'],\
                                  UNIPROT_COLLECTION.mget(data['uniprotIDs'], raw=False) 
                                ):
            results[ uniprot_id if not e else e.id ] = e.toJSON() if not e is None else None
            validCnt = validCnt + 1 if e is not None else validCnt

    print(f"Returning { validCnt } valid elements")
    return jsonify(results)

def get_missing():
    '''Return missing proteins'''
    print("Check availability")
    data = _load_uniprot_ids_request()
    redis_proteins = [e.id if e else None for e in UNIPROT_COLLECTION.mget(data['uniprotIDs'], raw=False)]
    if None in redis_proteins:
        missing_indices = [i for i, x in enumerate(redis_proteins) if x == None]
        missing_proteins = [data['uniprotIDs'][i] for i in missing_indices]
        return jsonify(missing_proteins)
    return jsonify([])

def model():
    #pyproteinsext.model()
    return "Hello world"