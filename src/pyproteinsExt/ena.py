import pyproteins.container.customCollection
import pyproteins.container.Core
from os.path import expanduser
import re
import os 
import io 
from skbio import DNA

enaEntrySet=None

def getENACollection():
    print("Get ENA Collection")
    global enaEntrySet
    enaEntrySet=EntrySet()
    return enaEntrySet    


class EntrySet(pyproteins.container.customCollection.EntrySet):
    def __init__(self, **kwargs):
        home = expanduser("~")
        cachePath = home
      #  if 'collectionPath' in kwargs:
      #      cachePath = kwargs['collectionPath']
      #  if 'pfamCollectionPath' in kwargs:
      #      PfamCache = kwargs['pfamCollectionPath']

        super().__init__(collectionPath=cachePath, constructor=Entry, indexer=strip)

    def serialize(self, ext=''):
        print ("serializing ena collection")
        super().serialize(ext=ext)

class Entry(pyproteins.container.Core.Container):
    def __init__(self, id, baseUrl="https://www.ebi.ac.uk/ena/data/view/", fileName=None, ext="&display.txt", charge_empty=False, rerun=False, url_id=None, **kwargs):
        if not id:
            raise TypeError('identifier is empty')
        if url_id: 
            url = baseUrl + str(url_id) + ext
        else: 
            url = baseUrl + str(id) + ext
        super().__init__(id, url=url, fileName=fileName)
        self.rerun = rerun
        self.type = type
        if charge_empty:
            self.features = []
        else:
            self.embl_parsing_features(self.raw, **kwargs)

            # Search assembly embl set if there is no features in classical url
            if len(self.features) <= 1 and not rerun:
                short_id = id[:6]
                first_letters = id[:2].lower()
                new_base_url = "http://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/"+first_letters+"/"
                ext = ".dat.gz"
                fileName = None
                self.__init__(self.id, new_base_url, fileName, ext, rerun=True, url_id=short_id, **kwargs)
     
    def embl_parsing_features(self, rawData, **kwargs):
        def conserve_feature(feature, type_filter, info_filter):
            if type_filter:
                if feature.type not in type_filter:
                    return False
            if info_filter:
                for k in info_filter:
                    if not feature.info.get(k):
                        return False
                    if not set(feature.info[k]).intersection(set(info_filter[k])):
                        return False
            return True

        keep_sequence = kwargs.get("keep_sequence", False)
        type_filter = kwargs.get("type_filter", [])
        info_filter = kwargs.get("info_filter", {})

        list_features = []
        # print("== PARSING")
        # print(kwargs)
        # print(keep_sequence, type_filter, info_filter)
        feature_header_re = re.compile("^FT\s{3}([^\s]+)\s+([^\s]+)")
        feature_element_re = re.compile("^FT\s+\/(.+)")
        feature_txt_re = re.compile("^FT\s{19}([^\s /].+)")

        if isinstance(rawData, bytes):
            rawData = rawData.decode()
        feature = None
        l_nb = 0
        dic_nb = {}
        nb_source = 0
        for l in rawData.split("\n"):
            l_nb += 1

            if l.startswith("FT"):  # Feature line
                test_header = feature_header_re.match(l)
                if test_header:
                    # Store previous feature if exists
                    if feature:
                        if conserve_feature(feature, type_filter, info_filter):
                            list_features.append(feature)
                    header = True
                    current_type = test_header.group(1)
                    location = test_header.group(2)
                    feature = Feature(current_type, location)   
                    # Set feature source if it's not source itself
                    if current_type == "source":
                        nb_source += 1
                        if self.rerun:
                            feature.name = "Contig"+str(nb_source)
                        else:
                            feature.name = "Genome_element"+str(nb_source)
                        source = feature
                        dic_nb = {}
                    else:
                        if current_type not in dic_nb:
                            dic_nb[current_type] = 0
                        dic_nb[current_type] += 1 
                        feature.source = source
                        feature.name = current_type+str(dic_nb[current_type])+"_"+feature.source.name

                test_element = feature_element_re.match(l)
                if test_element: # Element line. Example : "FT      \pseudo" or "FT      \protein_id=""
                    header = False
                    element = test_element.group(1).split("=")
                    # Complete feature info dictionnary
                    if len(element) == 1:
                        element_key = "other_qualifiers"
                        to_add = element[0].replace('"', '')
                    
                    elif len(element) > 1:
                        element_key = element[0]
                        to_add = element[1].replace('"', '')

                    if not keep_sequence and element_key == "translation":
                        continue

                    if element_key not in feature.info:
                        feature.info[element_key] = []
                    feature.info[element_key].append(to_add)

                test_txt = feature_txt_re.match(l)  # Line with the following text of the element
                if test_txt:
                    if header:
                        feature.location += test_txt.group(1).replace('"','')
                    else:
                        if not keep_sequence and element_key == "translation":
                            continue
                        feature.info[element_key][-1] += test_txt.group(1).replace('"', '')

        if conserve_feature(feature, type_filter, info_filter):
            list_features.append(feature)

        # print("LF",len(list_features))
        self.features = list_features

    def filter(self, fPredicat, **kwargs):
        # print("FILTER")
        # print(self.id)
        # print(self.fileName)
        new_Entry = Entry(self.id, self.id, charge_empty=True, fileName=self.fileName)
        for f in self.features:
            if fPredicat(f, **kwargs):
                new_Entry.features.append(f)
        return new_Entry

def strip(string):
    subString = re.search(".embl$", string)
    if subString:
        return string.split(".")[0]
    return None


class Feature():
    def __init__(self, type, location):
        self.type = type
        self.location = location
        self.source = None
        self.info = {}