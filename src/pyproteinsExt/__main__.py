"""Go ontology tree manipulation tool and microservice
Usage:
    pyproteinsExt service uniprot <xmlProteomeFile> [--silent] [--port]

Options:
  -h --help     Show this screen.
  <xmlProteomeFile> uniprot file location in xml format
  --silent  verbosity
"""
import os


from docopt import docopt
from .services import uniprot as uniprotService

arguments = docopt(__doc__)
print(arguments)

proteomeXML = arguments['<xmlProteomeFile>'] if arguments['<xmlProteomeFile>'] else None
apiPort = arguments['--port']     if arguments['--port'] else 5000

if arguments['service'] and arguments['uniprot']:
    app = uniprotService.startup(proteomeXML)
    app.run(debug=False, port=apiPort)
