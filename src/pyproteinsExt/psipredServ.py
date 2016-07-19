import paramiko
import re
import uuid
import subprocess
import os
import json
from cStringIO import StringIO
import utils

SOCKET = "/projet/extern/save/gulaunay/tmp/socket/ppsvr"
INPUT_FORMAT='^[\s]*([\d]+)[\s]+([\S]+)[\s]+([\S]+)[\s]+([0-9\.]+)[\s]+([0-9\.]+)[\s]+([0-9\.]+)'





def isValidRecord(line):
    #print line
    m = re.match(INPUT_FORMAT, line)
    if m:
        return True
    return False



class Socket():
    def __init__(self, url="migale.jouy.inra.fr", username="gulaunay"):
        self.client = paramiko.client.SSHClient()
        self.client.load_system_host_keys()
        self.client.connect(hostname=url, username=username)
        self.pool = {}
        self.user = username
        self.host = url
    def push(self, fileList=None, peptidesList=None, _blankShotID=None):
        if _blankShotID:
            name = _blankShotID
        else:
            name = uuid.uuid1().get_hex()
        self.pool[name] = {
            'socket' :  SOCKET + "/" + name
        }
        if not _blankShotID:
            self.client.exec_command("mkdir " + self.pool[name]['socket'])

        if fileList:
            self.pool[name]['inputs'] = self._pushFiles(name, fileList)
        elif peptidesList:
            self.pool[name]['inputs'] = self._pushPeptides(name, peptidesList)
        else:
            raise ValueError("No input to process")

        print 'Chunk content ' + str(self.pool[name]['inputs'])

        cmd = 'submitPsipred.sh ' + self.pool[name]['socket']  + ' ' + ' '.join([self.pool[name]['socket'] + '/' + os.path.basename(f) for f in self.pool[name]['inputs'] ])
        #print cmd + "\n"

        if _blankShotID:
            print "Not executing remote command:\"" + cmd + "\""
            return name


        ans = self.client.exec_command(cmd)
        s =  ans[2].read()
        if s:
            raise ValueError(s)
        return name

    def pull(self, chunkid):
        if chunkid not in self.pool:
            raise ValueError("unknwon id \"" + chunkid + "\"")
        data = []
        for i, e in enumerate(self.pool[chunkid]['inputs']):
            ss2File = self.pool[chunkid]['socket'] + "/qsub_" + str(i + 1) + "/input.ss2"
            cmd = "cat " + ss2File
            #print cmd
            ans = self.client.exec_command(cmd)
            data.append(collection(stream=ans[1]))
        return data

    def _pushFiles(self, name, fileList): # sending fasta files
        cmd = ['scp'] + fileList + [ self.user + '@' + self.host + ':' + self.pool[name]['socket'] ]
        #print cmd
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        return [ os.path.basename(f) for f in fileList ]

    def _pushPeptides(self, name, peptidesList): # dumping peptide object content to remote fasta file
        fList = []
        for i, e in enumerate(peptidesList):
            fname = self.pool[name]['socket'] + "/peptide_" + str(i) + ".fasta"
            cmd = "(echo \"" + e.fasta.replace("\n", "\" && echo \"") + "\") > " + fname # echo format should be moved to Peptide object
            print cmd
            ans = self.client.exec_command(cmd)

            print ans[1].read()
            fList.append(fname)
        return fList

class collection():
    def __init__(self, fileName=None, string=None, stream=None):

        if not fileName and not string and not stream:
            raise ValueError("No input provided")
        bufStream = None

        if fileName:
            bufStream = open (fileName, "r")
        elif string:
            bufStream = StringIO('foo')
        else :
            bufStream = stream

        self.data = [ element(line) for line in bufStream if isValidRecord(line) ]

        if utils.hasMethod(bufStream, "close"):
            bufStream.close()

    @property
    def horiz(self):
        return ''.join([e.state for e in self.data])
    @property
    def aaSeq(self):
        return ''.join([e.aa for e in self.data]).upper()

    @property
    def dict(self):
        return [ e.pos + " " + e.aa + " " + e.state + " " + e.cProb + " " + e.hProb + " " + e.eProb for e in self.data ]

class element():
    def __init__(self, line):
        lBuffer = line.split()
        self.pos = lBuffer[0]
        self.aa = lBuffer[1]
        self.state = lBuffer[2]
        self.cProb = lBuffer[3]
        self.hProb = lBuffer[4]
        self.eProb = lBuffer[5]
    def __repr__(self):
        return "\t".join([ str(val) for val in [self.pos, self.aa, self.state, self.cProb, self.hProb, self.eProb] ])
