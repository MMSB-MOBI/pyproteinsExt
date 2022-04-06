from repl_core.application import get_response, Application
from repl_core import run, Response
from repl_core import print_formatted_text, HTML
from ....uniprot import EntrySet
import json
import sys

def run_repl(host="locahost", port=2339, auto_connect=True):

    app = Application(port=port, host=host, route="/handshake", auto_connect=auto_connect)

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
    
    @app.viewer("/uniprot/list/{interval}",
                "list {interval:_string}",
                help_msg="list proteins in a given interval, expresses as a slice expression"
                )
    def list(interval=":1000"):
        print_formatted_text(f"Listing proteins in range {interval}")
        resp = get_response()
        color = "ansigreen" if resp.status_code == 200 else "ansired"
        if resp.content:
            print_formatted_text( HTML(f"<{color}>{json.dumps(resp.json(), indent=4)}</{color}>") )

    @app.mutator("/uniprot/put",
                "add {uniprot_xml_file:_path}",
                help_msg="Add one protein from xml file"
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


    def g(msg_chunks):
        ttl = 0
        err_stack = []
        for chunk in msg_chunks:
            d = chunk["results"]
            ttl += len(d["errors"]) + len(d["success"])           
            err_stack += [ e["id"] for e in d["errors"] ]
        
        if err_stack:
            m = f"<ansired>Failed to insert ({len(err_stack)} / {ttl}): {err_stack}</ansired>"
            print_formatted_text(HTML(m))
        else:
            print_formatted_text(HTML(f"<ansigreen>All {ttl} insertions successfull</ansigreen>"))

        return True
    @app.bulk("/uniprot/put_many",
                "add_many {uniprot_xml_file:_path} {start:_number} {stop:_number}",
                help_msg="Add proteins ranked start to stop from xml file",
                prefix="entrySet",
                size=2,
                validator=g
                )
    def add_many(uniprot_xml_file, start, stop):

        entrySet = EntrySet(collectionXML=uniprot_xml_file)
        print_formatted_text(f"Loading {len(entrySet)} uniprot objects from {uniprot_xml_file}")
        data_to_put = []
        summary_to_put = []
        for i,  e in enumerate(entrySet):
            if i >= start:
                data_to_put.append(e.toJSON())
                summary_to_put.append(e.id)
            if i == stop:
                break
        print_formatted_text(HTML(f"Bulk insert of <ansigreen>{summary_to_put}</ansigreen>"))
        
        return data_to_put

   

    @app.viewer("/uniprot/wipe",
                "wipe",
                help_msg="Wipe entiere database"
                )
    def wipe():
        print_formatted_text(f"Deleting entiere database")
        resp = get_response()
        color = "ansigreen" if resp.status_code == 200 else "ansired"
        if resp.content:
            print_formatted_text( HTML(f"<{color}>{json.dumps(resp.json(), indent=4)}</{color}>") )


    run(app)

# We neef a app.bulk which splice and push while tracking advances
