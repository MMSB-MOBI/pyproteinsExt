"""Uniprot ressources microservice

Usage:
  pyproteinsext service uniprot redis start [<xmlProteomeFile>] [--rh=<redis_host> --rp=<redis_port>] [--silent] [--port=<portNumber>]
  pyproteinsext service uniprot redis wipe [--rh=<redis_host> --rp=<redis_port>]
  pyproteinsext service uniprot xml start <xmlProteomeFile> [--silent] [--port=<portNumber>]
    
Options:
  -h --help     Show this screen.
  --port=<portNumber>  port for public API [default: 2332]
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
    redis=arguments['redis'],\
    rh=arguments['--rh'], rp=arguments['--rp'])

app.run(debug=False, port=arguments['--port'])
