import copy
import numpy as np
import json
from collections import OrderedDict as odict

import utils

############################################################

SISPI_DICT = odict([
    ("filter",  None),
    ("program", "des nu"),
    ("seqtot",  None),
    ("seqnum",  None),
    ("expType", "object"),
    ("object",  None),
    ("comment", None),
    #("exptime", None), # NOT NEEDED?
    ("note",    "Added to queue from desnu json file, not obstac"),
    ("seqid",   None),
    ("RA",      None),
    #("propid",  "2017B-0239"), # https://www.noao.edu/noaoprop/abstract.mpl?2017B-0239 (for DECam)
    ("propid",  "2019A-0240"), # https://www.noao.edu/noaoprop/abstract.mpl?2018B-0312 (for DECam)
    ("dec",     None),
    ("exptime", 150),
    ("wait",    "False"),
])

BANDS = ['z', 'r', 'g']
#BANDS = ['g', 'r', 'i'] # Normal
#BANDS = ['r', 'i', 'z'] # Special Moony setting
N_DITHER = 2
N_EXPOSURES = 1

############################################################

def writeJson(outfile, data, **kwargs):
    kwargs.setdefault('indent', 4)
    json.encoder.FLOAT_REPR = lambda o: format(o, '.4f')

    with open(outfile, 'wb') as out:
        # It'd be nice to have a header
        #out.write(header())
        out.write(json.dumps(data, **kwargs))

def readJson(filename, **kwargs):
    with open(filename, 'r') as f:
        return json.loads(f.read(), **kwargs)

############################################################

def makeJson(ra, dec, eventid, optimal_time):
    seqnum = 1
    seqtot = len(BANDS) * N_DITHER * N_EXPOSURES
    #seqid = 'IceCube_%i_%s'%(eventid, utils.datestring(optimal_time))
    seqid = 'IceCube_%i'%(eventid)

    sispi = []
    for band in BANDS:
        for index_dither in np.arange(N_DITHER):
            for index_exposure in np.arange(N_EXPOSURES):
                sispi_dict = copy.deepcopy(SISPI_DICT)

                tag = "DESNU: IceCube event %s: %i of %i"%(seqid, seqnum, seqtot)

                #sispi_dict["RA"]      = ra
                #sispi_dict["dec"]     = dec
                sispi_dict["RA"]      = ra + (index_dither * (1. / 60.) / np.cos(np.radians(dec)))
                sispi_dict["dec"]     = dec + (index_dither * (1. / 60.))
                sispi_dict["filter"]  = band
                sispi_dict["seqnum"]  = seqnum
                sispi_dict["seqtot"]  = seqtot
                sispi_dict["object"]  = tag
                sispi_dict["comment"] = tag
                sispi_dict["seqid"]   = seqid

                sispi.append(sispi_dict)
                seqnum += 1
    return sispi
    
############################################################


#seqid = '80127519'
#ra, dec = 45.8549, 15.7851
#optimal_time = 


"""
"count": 1, 
        "filter": "i", 
        "program": "des gw", 
        "seqtot": 84, 
        "seqnum": 1, 
        "expType": "object", 
        "object": "DESGW: LIGO NS event M249148: 1 of 84, hex ['10.16252' '9.14855' '8.13457' '-8.08903' '-5.0471' '6.10662']", 
        "comment": "DESGW: LIGO NS event M249148: 1 of 84, hex ['10.16252' '9.14855' '8.13457' '-8.08903' '-5.0471' '6.10662']", 
        "exptime": 90, 
        "note": "Added to queue from desgw json file, not obstac", 
        "seqid": "M249148", 
        "RA": -34.173876, 
        "propid": "2015B-0187", 
        "dec": 10.16252, 
        "expTime": 90, 
        "wait": "False"


data = [{'a': 1,
         'b': 2},
        {'c': 3,
         'd': 4}]

write_json('out.json', data)
"""



#FILENAME = optimal time
#M249148-6-UTC-2016-8-17-5-23-00.json
# The time in the file name is the "optimal time" to run this script.
# The  "object" field is parsed by diffimg pipeline.
# The "comment" field is parsed by our jobmanager code.

"""
SISPI_DICT = odict([
    ("object",  None),
    ("seqnum",  None), # 1-indexed
    ("seqtot",  2),
    ("seqid",   None),
    ("expTime", 90),
    ("RA",      None),
    ("dec",     None),
    ("filter",  None),
    ("count",   1),
    ("expType", "object"),
    ("program", "maglites"),
    ("wait",    "False"),
    ("propid",  "2016A-0366"),
    ("comment", ""),
])

SISPI_MAP = odict([ 
    ('expTime','EXPTIME'),
    ('RA','RA'),
    ('dec','DEC'),
    ('filter','FILTER'),
])

def to_sispi(self):
        sispi = []
        objects = self.object
        seqnums = self.seqnum
        seqids = self.seqid
        comments = self.comment
        for i,r in enumerate(self):
            sispi_dict = copy.deepcopy(SISPI_DICT)
            for sispi_key,field_key in SISPI_MAP.items():
                sispi_dict[sispi_key] = r[field_key]
            sispi_dict['object'] = objects[i]
            sispi_dict['seqnum'] = seqnums[i]
            sispi_dict['seqid']  = seqids[i]
            sispi_dict['comment'] = comments[i]
            sispi.append(sispi_dict)
        return sispi

def write(self, filename, **kwargs):
        base,ext = os.path.splitext(filename)
        logging.debug('Writing %s...'%filename)
        if ext in ('.json'):
            data = self.to_sispi()
            fileio.write_json(filename,data,**kwargs)
        elif ext in ('.csv','.txt'):
            data = self.to_recarray()
            fileio.rec2csv(filename,data,**kwargs)
        else:
            msg = "Unrecognized file extension: %s"%ext
            raise IOError(msg)
"""
