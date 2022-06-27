from ..models.entryProxy import ProxyEntry
from typing import List

class EntrySetProxy():
    """
    A proxy for pyproteinsext.uniprot.EntrySet
    """
    def __init__(self, entry_list:List[ProxyEntry]):
        self.data = entry_list
        self._index = { e.id:e for e in self.data }

    def __len__(self):
        return len(self.data)

    def keys(self):
        for e in self.data:
            yield(e.id)
    
    def has(self, maybe_ac):
        return maybe_ac in self.index
    
    def __iter__(self):
        for entry in self.data:
            yield(entry)

    def get(self, uniprotID):
        if uniprotID in self.index:
            return self.index[uniprotID]
        return None
    
    @property
    def taxids(self):
        taxids = set()
        for e in self:
            taxids.add(e.taxid)
        return list(taxids)
