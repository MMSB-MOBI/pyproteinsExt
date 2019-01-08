import re
import pyproteinsExt.uniprot
import copy
import math
import pyproteinsExt.matrisome
from queue import Queue
from threading import Thread
from collections import OrderedDict
import pandas as pd
import numpy as np

import copy

# Basic input is a pdFrame, CSV maybe

# Takes an iterable, along with a coherce to uniprot Identifier
# returns and list of tuple with each iterable element and its uniprot Object
#
# Iterable ---->(Iterable, Uniprot(reduce(Iterable))

class pandaUniprotEnricher(object):
    def __init__(self, pdFrame=None, lambdaFn=None, **kwargs):
        self.pd = pdFrame
        if 'n' in kwargs:
            self.pd = pdFrame[:n]
        self.lambdaFn = lambdaFn
        if 'uniAno' in kwargs:
            self.data = kwargs['uniAno']
        else:
            self.data = [ UniAno(lambdaFn(e)) for index, e in self.pd.iterrows() ]
            print("Pooling")
            _threadPool(self.data);

        noUniprotBound = []
        _data = []
        for i,e in enumerate(self.data):
            if not e._uniprotBound:
                noUniprotBound.append(i)
                continue
            _data.append(e)

        if noUniprotBound:
            print ("Following table elements will be dropped due to unbound uniprot object")
            print ( str([self.data[n].uniprotID for n in noUniprotBound ]) )
        #print (noUniprotBound)
        #for n in reversed(noUniprotBound):
        #    print ("Droping table element " + str(n) + " (" + self.data[n].uniprotID + ") due to unbound uniprot object")
        df = self.pd.drop(self.pd.index[noUniprotBound]).reset_index(drop=True)
        self.pd = df
        self.data = _data  

    def __len__(self):
        if len(self.data) != len(self.pd):
            raise ValueError('data length dont match')
        return len(self.data)


# We select out of the padaframe only entries that match annotator
#https://chrisalbon.com/python/data_wrangling/pandas_selecting_rows_on_conditions/

# Create a dataframe
#raw_data = {'first_name': ['Jason', 'Molly', np.nan, np.nan, np.nan],
#        'nationality': ['USA', 'USA', 'France', 'UK', 'UK'],
#        'age': [42, 52, 36, 24, 70]}
#df = pd.DataFrame(raw_data, columns = ['first_name', 'nationality', 'age'])
#df

# Merging dataframes
#https://pandas.pydata.org/pandas-docs/stable/merging.html

    def filter(self, annotator, collapseCol=True, exclude=True):

        _annot, status = annotator.annotateAll(self.data)
        # fix a column order
        colHeadOrdered = [ term for term in _annot ]
        # stringify terms occurences


        #print (_annot)


        for term in _annot:
            _annot[term] = [ ','.join([ str(tOcc) for tOcc in terms ]) for terms in _annot[term] ]
       
        #return _annot
        df = self.pd
        uniAnoList =[]
        # reduce list of annotations and original dataframe
        if exclude:
            annot = {}
            for term in _annot:
                annot[term] = [ e for i,e in enumerate(_annot[term]) if status[i] ]
            print ('Getting rid of ' + str(status.count(False)) + ' entries')
            df = df.drop(df.index[ [ i for i,e in enumerate(status) if not e ] ])
            uniAnoList = [ self.data[i] for i,e in enumerate(status) if e ]
            
            #return df
        else:
            annot = _annot
            uniAnoList = self.data

        if collapseCol:
           
            _df = pd.DataFrame(annot, columns = colHeadOrdered)
            df.reset_index(drop=True, inplace=True)
            _df.reset_index(drop=True, inplace=True)

         
            df = pd.concat([df, _df], axis=1)
         
        if exclude:
            df = df.replace('', np.nan, regex=True).dropna(axis=1,how='all').replace(np.nan, '', regex=True)
        return pandaUniprotEnricher(pdFrame=df, lambdaFn = self.lambdaFn, uniAno=uniAnoList)
# panda integer indexing
# https://stackoverflow.com/questions/16096627/selecting-a-row-of-pandas-series-dataframe-by-integer-index
    def __getitem__(self,i):
        return ( self.data[i], self.pd.iloc[[i]] )

class UniAno(object):
    def __init__(self, _id):        
        self.uniprotID = _id
        self._uniprotBound = None

    #def __str__(self):
    #    return str(self._uniprotBound);

    def _boundUniprot(self):
        if not self._uniprotBound:
            id = self.uniprotID

            if not id:
                return False
            try :
                notIsoform_id = pyproteinsExt.uniprot.strip(id)
                self._uniprotBound = uniprotEntrySet.get(notIsoform_id)#Uniprot.Entry(id)
            except TypeError as msg:
                print ("Following element could not be bound to uniprot entity, reason " + str(msg))
                print (self.uniprotID)
                return False
            except ValueError as msg:
                print ("Following element could not be bound to uniprot entity, reason " + str(msg))
                print (self.uniprotID)
                return False
        return True

def _threadPool(stack):
    num_worker_threads = 10
    def worker():
        while True:
            e = q.get()
            e._boundUniprot()
            q.task_done()

    q = Queue()
    for i in range(num_worker_threads):
        t = Thread(target=worker)
        t.daemon = True
        t.start()

    for e in stack:
        #print("Pouet")
        q.put(e)


    print ('*** Main thread waiting ' + str(q.qsize()))
    q.join()
    print ('*** Done')



