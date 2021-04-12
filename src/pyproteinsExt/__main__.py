"""Go ontology tree manipulation tool and microservice

Usage:
    pyproteinsExt service uniprot start [<xmlProteomeFile>] [--silent] [--port=<portNumber>] [--redis --rh=<redis_host> --rp=<redis_port>]
    pyproteinsExt service uniprot redis wipe [--rh=<redis_host> --rp=<redis_port>]

Options:
  -h --help     Show this screen.
  <xmlProteomeFile>  uniprot file location in xml format
  --port=<portNumber>  port for public API [default: 2332]
  --redis  use redis backend as storage engine
  --rp=<redis_port>  redis DB TCP port [default: 6379]
  --rh=<redis_host>  redis DB http adress [default: localhost]
  --silent  verbosity
  
"""
import os


from docopt import docopt
from .services import uniprot as uniprotService

arguments = docopt(__doc__)
print(arguments)

if arguments['wipe']:
  uniprotService.cleanup(rh=arguments['--rh'], rp=arguments['--rp'])
  exit(1)

app = uniprotService.startup(arguments['<xmlProteomeFile>'],\
    redis=arguments['--redis'],\
    rh=arguments['--rh'], rp=arguments['--rp'])

app.run(debug=False, port=arguments['--port'])
