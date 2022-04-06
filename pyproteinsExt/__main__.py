"""Uniprot ressources microservice

Usage:
  pyproteinsext service uniprot redis start [<xmlProteomeFile>] [--rh=<redis_host> --rp=<redis_port>] [--silent] [--port=<portNumber>]
  pyproteinsext service uniprot redis wipe [--rh=<redis_host> --rp=<redis_port>]
  pyproteinsext service uniprot xml start <xmlProteomeFile> [--silent] [--port=<portNumber>]
  pyproteinsext service uniprot repl [--host=<http_host> --port=<http_port>]

Options:
  -h --help     Show this screen.
  --port=<http_port>  port for public API [default: 2332]
  --host=<http_host>  hostname for public API [default: localhost]
  --rp=<redis_port>  redis DB TCP port [default: 6379]
  --rh=<redis_host>  redis DB http adress [default: localhost]
  --silent  verbosity
  
"""
import os


from docopt import docopt
from .services import uniprot as uniprotService
from .services.uniprot.client.repl import run_repl
arguments = docopt(__doc__)


if arguments['repl']:
  run_repl(port=arguments['--port'], host=arguments['--host'])


if arguments['wipe']:
  uniprotService.cleanup(rh=arguments['--rh'], rp=arguments['--rp'])
  exit(1)

app = uniprotService.startup(arguments['<xmlProteomeFile>'],\
    redis=arguments['redis'],\
    rh=arguments['--rh'], rp=arguments['--rp'])

app.run(debug=False, port=arguments['--port'])
