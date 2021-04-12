from flask import Flask, jsonify, abort, request, jsonify
from pyproteinsExt import uniprot as pExt
from . import collectionProxy as redisCollection

UNIPROT_COLLECTION=None

def cleanup(rh=None, rp=None):
    redisCollection.bootstrap(host=rh, port=rp)
    redisCollection.cleanup()

def startup(xmlUniprot, redis=False, rh=None, rp=None):

    global UNIPROT_COLLECTION
   # UNIPROT_COLLECTION = pExt.EntrySet(collectionXML=xmlUniprot)
   
    if xmlUniprot:
        print(f"Loadind XML ressource {xmlUniprot}")
        _ = pExt.EntrySet(collectionXML=xmlUniprot)
        print(f"Loaded uniprot collection from {xmlUniprot} elements")

    if redis:
        redisCollection.bootstrap(host=rh, port=rp)
        
        if xmlUniprot:
            redisCollection.convert(_)
        UNIPROT_COLLECTION = redisCollection
    if not xmlUniprot:
        print("No data added to uniprot storage service")

    
    print("Uniprot storage service listening")
 
    app = Flask(__name__)

    app.add_url_rule('/model', 'model', model)
    
    app.add_url_rule( '/uniprot/<uniprotID>', 'getProtein', getProtein,
                      methods = ['GET'] ) 
    
    app.add_url_rule( '/uniprots', 'getProteins', getProteins,
                      methods = ['POST'] )

    app.add_url_rule( '/length', 'length', length,
                      methods = ['GET'] ) 
    
    #app.add_url_rule('/unigo/<taxid>', 'view_unigo', view_unigo)
    
    return app

def length():
    global UNIPROT_COLLECTION
    length = len(UNIPROT_COLLECTION)
    print(f"Current Collection size ${length}")
    return jsonify( {"totalEntry" : length } )

def getProtein(uniprotID):
    print(f"Seeking {uniprotID}")
    oProtein = UNIPROT_COLLECTION.get(uniprotID)
    if not oProtein:
        abort(404)
    return jsonify( { uniprotID : oProtein.toJSON() } )

def getProteins():
    if not request.is_json:
        print("Post data not json")
        abort(422)
    data = request.get_json()
    if not 'uniprotIDs' in data:
        print("Post data lacks uniprotIDs key")
        abort(422)
    
    results = {}
    validCnt = 0
    for id in data['uniprotIDs']:
        _ = UNIPROT_COLLECTION.get(id)
        results[id] = _.toJSON() if not _ is None else None
        validCnt = validCnt + 1 if results[id] else validCnt
    print(f"Returning { validCnt } valid elements")
    return jsonify(results)

def model():
    #pyproteinsExt.model()
    return "Hello world"