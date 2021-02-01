from flask import Flask, jsonify, abort, request, jsonify
from pyproteinsExt import uniprot as pExt

UNIPROT_COLLECTION=None

def startup(xmlUniprot):
    print("Loading")

    load(xmlUniprot)
 
    app = Flask(__name__)

    app.add_url_rule('/model', 'model', model)
    
    app.add_url_rule( '/uniprot/<uniprotID>', 'getProtein', getProtein,
                      methods = ['GET'] ) 
    
    app.add_url_rule( '/uniprots', 'getProteins', getProteins,
                      methods = ['POST'] )
    
    #app.add_url_rule('/unigo/<taxid>', 'view_unigo', view_unigo)
    
    return app

def load(xmlFile):
    global UNIPROT_COLLECTION
    UNIPROT_COLLECTION = pExt.EntrySet(collectionXML=xmlFile)
    print(f"Loaded uniprot collection from {xmlFile} elements")

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
    for id in data['uniprotIDs']:
        _ = UNIPROT_COLLECTION.get(id)
        results[id] = _.toJSON() if not _ is None else None
    
    print(f"Returning {results}")
    return jsonify(results)

def model():
    #pyproteinsExt.model()
    return "Hello world"