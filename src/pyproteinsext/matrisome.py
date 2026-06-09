import pyproteins.services.utils
class Matrisome:
    def __init__(self, masterFile=None):
        if not masterFile:
            raise ValueError ("You must provide a masterFile argument to Matrisome constructor")

        tsvBuffer = pyproteins.services.utils.tsvToDictList(fileName=masterFile)
        self.data = tsvBuffer['data']
        self.keymap = tsvBuffer['keymap']
        self.accessors = self._index()

    # index entry using relevant key found in their record
    def _index(self):
        accessors = {'uniprot' : {}, 'gene' : {}} # any other meaningfull accessor
        for d in self.data:
            #print "-->" + d['UniProt_IDs'] + '<--'
            if d['UniProt_IDs']:
                for k in d['UniProt_IDs'].split(':'):
                    if k in accessors['uniprot']:
                        print (k + ' found multiple time')
                    else :
                        accessors['uniprot'][k] = []
                    accessors['uniprot'][k].append(d)
            if d['Gene Symbol']:
                for k in d['Gene Symbol'].split(':'):
                    if k in accessors['gene']:
                        print (k + ' found multiple time')
                    else :
                        accessors['gene'][k] = []
                    accessors['gene'][k].append(d)
        return accessors

    def get(self, uniprotID=None):
        if uniprotID:
            if not uniprotID in self.accessors['uniprot']:
                return []
            return [{'Category' : d['Category'], 'Division' : d['Division']} for d in self.accessors['uniprot'][uniprotID]]
        return []
