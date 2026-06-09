class Container(object):
    def __init__(self, fParsing=None, input=None):
        if not input:
            self.entries = {}
        else:
            self.entries = self.parsing(input, fParsing)

    def addParsing(self, other):
        for k in other.entries:
            if k not in self.entries:
                self.entries[k] = other.entries[k]
        return self

    def addEntry(self, new_entry):
        if new_entry.prot not in self.entries:
            self.entries[new_entry.prot] = new_entry

    def __len__(self):
        return len(self.entries)

    def __iter__(self):
        for protein in self.entries:
            yield self.entries[protein]

    def parsing(self, input, fParsing):
        return fParsing(input)
