from repl_core.application import get_response, Application
from repl_core import run, Response
from repl_core import print_formatted_text, HTML
from ....uniprot import EntrySet
import json

app = Application(port=2332, route="/handshake", auto_connect=True)


@app.viewer("/uniprot/{uniprot_id}",
            "find {uniprot_id:_string}",
            help_msg="Find a protein"
            )
def find(uniprot_id):
    print_formatted_text(f"Searching for protein ID {uniprot_id}")
    resp = get_response()
    color = "ansigreen" if resp.status_code == 200 else "ansired"
    if resp.content:
        print_formatted_text( HTML(f"<{color}>{json.dumps(resp.json(), indent=4)}</{color}>") )

@app.mutator("/uniprot/put",
            "add {uniprot_xml_file:_path}",
            help_msg="Add proteins from xml file"
            )
def add(uniprot_xml_file):
    def process(response:Response):   
        color = "ansigreen" if response.status_code == 200 else "ansired"     
        print_formatted_text(HTML(f"<{color}>{response.text}</{color}>"))
    entrySet = EntrySet(collectionXML=uniprot_xml_file)
    print_formatted_text(f"Loading {len(entrySet)} uniprot objects from {uniprot_xml_file}")
    data_to_post = None
    for e in entrySet:
        data_to_post = e.toJSON()
        break
    
    print_formatted_text(HTML(f"Sending <ansigreen>{data_to_post}</ansigreen>"))
    
    return data_to_post, process



def run_repl():
    run(app)

# We neef a app.bulk which splice and push while tracking advances
