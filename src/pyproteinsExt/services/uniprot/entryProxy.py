from marshmallow import Schema, fields, post_load, INCLUDE
from ...uniprot import isValidID, GoKW# impoort unirpiotID validator

class GoTermSchema(Schema):
    id       = fields.Str()
    term     = fields.Str()
    evidence = fields.Str()
    @post_load
    def make_GoKW(self, data, **kwargs):
        return GoKW(**data)

class EntrySchema(Schema):
    #class Meta:
    #    unknown = INCLUDE
    name = fields.Str()
    id   = fields.Str(validate=isValidID)
    geneName = fields.Str(allow_none=True)
    fullName = fields.Str()
    taxid    = fields.Str() # Int ?
    GO = fields.List(fields.Nested(GoTermSchema))
    @post_load
    def make_uniprotEntryProxy(self, data, **kwargs):
        return ProxyEntry(**data)



class ProxyEntry():
    def __init__(self, id, name, geneName, fullName, taxid, GO):
        self.name     = name
        self.id       = id
        self.geneName = geneName
        self.fullName = fullName
        self.taxid    = taxid
        self.GO = GO
    def toJSON(self): # duplicate of pyproteins.Entry.toJSON ...
        container = {}
        for k, v in self.__dict__.items():
            if k == 'name':
                container[k] = v
            if k == 'GO':
                container[k] = [ go.__dict__ for go in v ]
            if k == 'id':
                container[k] = v
            if k == 'geneName':
                container[k] = v
            if k == 'fullName':
                container[k] = v
            if k == 'taxid':
                container[k] = v
        return container
    
    def hasGO(self):
        return keyword.upper() in (kw.id.upper() for kw in self.GO)
    
    def __repr__(self):
        #asStr = f"{self.id}:{self.AC}\n" 
        asStr = f"{self.id}:AC are not yet serialized\n" 
        asStr += f"{self.name}:{self.fullName}({self.geneName})\n"
        
        #if self.STRING_ID:
        #    asStr += f"STRING_ID:{self.STRING_ID}\n"
        asStr += f"taxid:{self.taxid}"
        #asStr += f"taxid:{self.taxid}:{self.lineage}\n"
        #asStr += f"KW:{self.KW}\n"
        asStr += f"GO:{self.GO}\n"
        
        return asStr