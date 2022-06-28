import pyproteinsext.hmmrContainerFactory as hmmr
import pyproteinsext.tmhmmContainerFactory as tmhmm
import pyproteinsext.fastaContainerFactory as fasta
import pyproteinsext.proteinContainer
from collections import OrderedDict
from ete3 import NCBITaxa
from statistics import mean
from igraph import Graph
import copy


def check_if_same_proteins(dic_container):
    hmmr_proteins = set()
    tmhmm_proteins = set()
    fasta_proteins = set()
    if dic_container['hmmr']:
        hmmr_proteins = set([h.prot for h in dic_container['hmmr'].hmmrEntries])
    if dic_container['tmhmm']:
        tmhmm_proteins = set([e.prot for e in dic_container['tmhmm']])
    if dic_container['fasta']:
        fasta_proteins = set([e.prot for e in dic_container['fasta']])

    list_proteins = [hmmr_proteins, tmhmm_proteins, fasta_proteins]
    if len(list_proteins) <= 1:
        return True
    for i in range(len(list_proteins)):
        for j in range(i+1, len(list_proteins)):
            prot1 = list_proteins[i]
            prot2 = list_proteins[j]
            if len(prot1) != len(prot2):
                return False
            if prot1.difference(prot2):
                return False
    return True


def parse(hmmrOut, tmhmmOut, fastaOut):
    container = TopologyContainer()
    hmmrContainer = hmmr.parse(hmmrOut)
    tmhmmContainer = tmhmm.parse(tmhmmOut)
    fastaContainer = fasta.parse(fastaOut)

    dic_container = {'hmmr': hmmrContainer, 'tmhmm': tmhmmContainer, 'fasta': fastaContainer}
    if not check_if_same_proteins(dic_container):
        raise Exception("not same proteins ", hmmrOut, tmhmmOut, "Check !")
    container.addParsing(TopologyContainer(input=dic_container))
    return container


def _parseBuffer(dic_container):
    hmmrContainer = dic_container['hmmr']
    tmhmmContainer = dic_container['tmhmm']
    fastaContainer = dic_container['fasta']
    dic_obj = {}
    for f in fastaContainer:
        p = f.prot
        hmmr = [h for h in hmmrContainer.hmmrEntries if h.prot == p]
        tmhmm = tmhmmContainer.entries[p]
        fasta = fastaContainer.entries[p]
        obj = Topology(p, hmmr, tmhmm, fasta)
        dic_obj[p] = obj
    return dic_obj


class TopologyContainer(pyproteinsext.proteinContainer.Container):
    def __init__(self, input=None, ete3_tree=None, domain_entries=None):
        super().__init__(_parseBuffer, input)
        self.ete3_tree = ete3_tree
        self.domain_entries = domain_entries

    def filter(self, fPredicat, **kwargs):
        new_container = TopologyContainer()
        for e in self:
            if fPredicat(e, **kwargs):
                new_container.addEntry(e)
        return new_container

    def __getitem__(self, index):
        return list(self.entries.values())[index]

    def get_domain_mfasta(self, domain):
        mfasta = ''
        for e in self:
            for hit in e.hmmr:
                if hit.domain == domain:
                    header = ">" + hit.prot + " " + hit.domain
                    seq = self.get_seq(hit.hit)
                    mfasta += header + "\n" + seq + "\n"
        return mfasta

    def get_seq(self, hit):
        seq = hit.aliStringLetters
        seq = seq.replace("-", "")
        seq = seq.upper()
        return seq

    def proteins_mfasta(self):
        mfasta = ''
        for e in self:
            mfasta += ">" + e.fasta.header + "\n" + e.fasta.seq + "\n"
        return mfasta

    def complete_hmmr(self, hmmscan_out):
        container = hmmr.parse(hmmscan_out)
        new_proteins = set([h.prot for h in container.hmmrEntries])
        self_proteins = set(self.entries.keys())
        if len(new_proteins) != len(self_proteins):
            raise Exception("full Pfam proteins and original proteins are not the same")
        else:
            if new_proteins.difference(self_proteins):
                raise Exception("full Pfam proteins and original proteins are not the same")
        for e in self:
            hits = [h for h in container.hmmrEntries if h.prot == e.prot]
            e.hmmr += hits

    def filter_hit(self, fPredicat, **kwargs):
        new_container = TopologyContainer()
        for e in self:
            add = False
            hit_to_add = []
            for hit in e.hmmr:
                if fPredicat(hit, **kwargs):
                    add = True
                    hit_to_add.append(hit)
            if add:
                new_e = Topology(e.prot, hit_to_add, e.tmhmm, e.fasta, e.taxo,
                                 e.uniprot_entry,
                                 e.annotated_domains_fragments,
                                 e.Nter_UR_fragment,
                                 e.Cter_UR_fragment, e.helix_fragments,
                                 e.loop_fragments)
                new_container.addEntry(new_e)
        return new_container

    def filter_last_helix(self, distance=90):
        new_container = TopologyContainer()
        c = 0
        for e in self:
            new_helix_fragments = copy.deepcopy(e.helix_fragments)
            new_loop_fragments = copy.deepcopy(e.loop_fragments)
            if len(e.helix_fragments) == 7:
                c += 1
                new_helix_fragments = new_helix_fragments[:-1]
                new_loop_fragments = new_loop_fragments[:-1]
            elif e.helix_fragments[-1]["start"]-e.helix_fragments[-2]["end"] > distance:
                c += 1
                new_helix_fragments = new_helix_fragments[:-1]
                new_loop_fragments = new_loop_fragments[:-1]
            new_e = Topology(e.prot, e.hmmr, e.tmhmm, e.fasta, e.taxo,
                             e.uniprot_entry,
                             e.annotated_domains_fragments,
                             e.Nter_UR_fragment, e.Cter_UR_fragment,
                             new_helix_fragments, new_loop_fragments,)
            new_container.addEntry(new_e)
        print(c, "proteins have been filtered")
        return new_container

    def compute_overlapped_domains(self, overlap_accept_size):
        self.reinitialize_overlapped_domains()
        for e in self:
            for i in range(len(e.hmmr)):
                for j in range(i+1, len(e.hmmr)):
                    hit1 = e.hmmr[i]
                    hit2 = e.hmmr[j]
                    if hit1.is_overlapping(hit2, overlap_accept_size):
                        hit1.overlapped_hits.append(hit2)
                        hit2.overlapped_hits.append(hit1)

    def reinitialize_overlapped_domains(self):
        for h in [h for e in self for h in e.hmmr]:
            h.reinitialize_overlapped_hits()

    def create_ete3_tree(self):
        ncbi = NCBITaxa()
        taxids = set([e.taxo.taxid for e in self])
        if None in taxids:
            raise Exception("Entries doesn't have taxids")
        tree = ncbi.get_topology(list(taxids))

        # Complete Tree object with list of domains and proteins for each node
        node_list = []
        for n in tree.traverse('postorder'):  # Browse tree in postorder, starts from leaf and ascend to root
            n.sameDomainNode = set()
            node_list.append(n)
            n.domains = set([h.domain for e in self for h in e.hmmr if e.taxo.taxid == n.name])
            n.proteins = set([e.prot for e in self if e.taxo.taxid == n.name])
            if n.get_descendants():
                for child in n.children:
                    n.domains.update(child.domains)
                    n.proteins.update(child.proteins)

        # Complete Tree object with list of nodes with same domains for each node
        c = 0
        for i in range(len(node_list)):
            c += 1
            for j in range(i+1, len(node_list)):
                n1 = node_list[i]
                n2 = node_list[j]
                if len(n1.domains) == len(n2.domains):
                    if not n1.domains.difference(n2.domains):
                        n1.sameDomainNode.add(n2)
                        n2.sameDomainNode.add(n1)
        self.ete3_tree = tree

    def create_domain_entries(self):
        def initialize_domain_entries():
            domain_entries = {}
            domains = set([h.domain for e in self for h in e.hmmr])
            for d in domains:
                domainObj = Domain(d, set(), set(), set())
                domain_entries[d] = domainObj
            return domain_entries
        domain_entries = initialize_domain_entries()
        for e in self:
            for h in e.hmmr:
                domain_entries[h.domain].hits.add(h)
                domain_entries[h.domain].proteins.add(e.prot)
                domain_entries[h.domain].taxo.add(e.taxo)
        domain_entries = OrderedDict(sorted(domain_entries.items(),
                                     key=lambda kv: len(kv[1].proteins),
                                     reverse=True))
        self.domain_entries = domain_entries

    def compute_upper_node_and_distance(self, core_domains=[]):
        ncbi = NCBITaxa()
        if not self.domain_entries:
            raise Exception("Compute domain_entries first.")
        if not self.ete3_tree:
            raise Exception("Compute ete3_tree first.")

        for d in self.domain_entries.values():
            if d.name not in core_domains:
                distances = []
                if len(d.taxo) == 1:
                    taxo = list(d.taxo)[0]
                    d.upper_node = self.ete3_tree.search_nodes(name=taxo.taxid)[0]
                    d.mean_distance = 0
                else:
                    list_taxids = list(set([t.taxid for t in d.taxo]))
                    domain_tree = ncbi.get_topology(list_taxids)
                    traverse_generator = domain_tree.traverse()
                    d.upper_node = next(traverse_generator)
                    for i in range(len(list_taxids)):
                        for j in range(i+1, len(list_taxids)):
                            dist = self.ete3_tree.get_distance(list_taxids[i], list_taxids[j])
                            distances.append(dist)
                    d.mean_distance = mean(distances)

    def create_domain_graph(self, core_domains):

        def get_vertex_size(nb_domain, interval, min_size):
            if nb_domain == 1:
                return min_size
            else:
                size = min_size+interval*(nb_domain-1)
                return size

        def get_edge_size(nb_occ, interval):
            return nb_occ*interval

        dic_edges = {}
        dic_nb_domain = {}
        all_domains = set()
        for p in self:
            domains = set([h.domain for h in p.hmmr])
            for cd in core_domains:
                domains.discard(cd)
            related_domains = domains.copy()
            for d in domains:
                if d not in core_domains:
                    if d not in dic_nb_domain:
                        dic_nb_domain[d] = 0
                    dic_nb_domain[d] += 1
                related_domains.remove(d)
                all_domains.add(d)
                for d2 in related_domains:
                    edge = tuple(sorted((d, d2)))
                    if edge not in dic_edges:
                        dic_edges[edge] = 0
                    dic_edges[edge] += 1

        list_edges = []
        list_weight = []
        for e in dic_edges:
            list_edges.append(e)
            list_weight.append(dic_edges[e])

        g = Graph()
        g.add_vertices(len(all_domains))
        g.vs["name"] = list(all_domains)
        g.vs["label"] = g.vs["name"]
        g.add_edges(list_edges)
        for vertex in g.vs:
            vertex["weight"] = dic_nb_domain[vertex["name"]]
            vertex["size"] = get_vertex_size(vertex["weight"], 2, 5)
        for e in g.es:
            source = [v for v in g.vs if v.index == e.source][0]
            target = [v for v in g.vs if v.index == e.target][0]
            edge_tuple = tuple(sorted((source["name"], target["name"])))
            nb_occ = dic_edges[edge_tuple]
            min_domains = min(source["weight"], target["weight"])
            e["weight"] = nb_occ/min_domains
            e["label"] = round(nb_occ/min_domains, 2)
            e['width'] = e["weight"]*5

        return g

    def separate_seq_into_fragments(self):
        for e in self:
            e.set_annotated_domains_fragments()
            e.set_Nter_UR_fragment()
            e.set_Cter_UR_fragment()
            e.set_helix_fragments()
            e.set_loop_fragments()

    def add_neighborhood_clusters(self, cluster_tsv):
        cluster_nb = 0
        browse_representative = set()
        f = open(cluster_tsv, "r")
        for l in f:
            l_split = l.rstrip().split("\t")
            representative = l_split[0].split(" ")[0].strip('"')
            seq = l_split[1].split(" ")[0].strip('"')
            prot = seq.split("+")[0]
            seq_name = seq.split("+")[1]
            if representative not in browse_representative:
                cluster_nb += 1
                browse_representative.add(representative)

            if not hasattr(self.entries[prot], "neighborhood_clusters"):
                self.entries[prot].neighborhood_clusters = {}

            self.entries[prot].neighborhood_clusters[seq_name] = cluster_nb
        f.close()

        for e in self:
            if not hasattr(e, "neighborhood_clusters"):
                e.neighborhood_clusters = {}

    def neighborhood_common_matrix(self):
        def get_common_size(list1, list2):
            common_size = 0
            for n in list1:
                if n in list2:
                    common_size += 1
            return common_size

        matrix = {}
        for i in range(100):
            for j in range(i+1, 100):
                neighbors_cluster_list_1 = self[i].get_neighborhood_clusters_number()
                neighbors_cluster_list_2 = self[j].get_neighborhood_clusters_number()
                if not self[i].prot in matrix: 
                    matrix[self[i].prot] = {}
                matrix[self[i].prot][self[j].prot] = get_common_size(neighbors_cluster_list_1, neighbors_cluster_list_2)

        return matrix            
                

class Topology():
    def __init__(self, prot, hmmr, tmhmm, fasta, taxo=None, uniprot_entry=None,
                 annotated_domains_fragments=None, Nter_UR_fragment=None,
                 Cter_UR_fragment=None, helix_fragments=None,
                 loop_fragments=None):
        self.prot = prot
        self.hmmr = hmmr
        self.tmhmm = tmhmm
        self.fasta = fasta
        self.taxo = taxo
        self.uniprot_entry = uniprot_entry
        self.annotated_domains_fragments = annotated_domains_fragments
        self.Nter_UR_fragment = Nter_UR_fragment
        self.Cter_UR_fragment = Cter_UR_fragment
        self.helix_fragments = helix_fragments
        self.loop_fragments = loop_fragments

    def get_taxo(self, function_get_taxid):
        ncbi = NCBITaxa()
        taxname = None
        taxrank = None
        taxid = function_get_taxid(self)
        taxname_dic = ncbi.get_taxid_translator([taxid])
        if taxname_dic:
            taxname = taxname_dic[int(taxid)]
            taxrank_dic = ncbi.get_rank([taxid])
            if taxrank_dic:
                taxrank = taxrank_dic[int(taxid)]
        self.taxo = Taxo(taxid, taxname, taxrank)

    def set_annotated_domains_fragments(self):
        annotated_domains_fragments = []
        for hit in self.hmmr:
            dic_core_domains = {}
            if hit.domain in dic_core_domains:
                raise Exception("Several hits for domain. Handle this part")
            dic_core_domains = {'name': hit.domain, 'seq': hit.sequence,
                                'start': hit.start, 'end': hit.end}
            annotated_domains_fragments.append(dic_core_domains)
        annotated_domains_fragments.sort(key=lambda r: r["start"])
        self.annotated_domains_fragments = annotated_domains_fragments

    def set_Nter_UR_fragment(self):
        dic = {'name': 'N-ter_UR', 'start': 1}
        Nter = self.tmhmm.fragments[0]
        if Nter.end > self.annotated_domains_fragments[0]["start"]:
            dic["end"] = self.annotated_domains_fragments[0]["start"]
        else:
            dic["end"] = Nter.end
        dic["seq"] = self.fasta.get_subsequence(dic["start"], dic["end"])
        self.Nter_UR_fragment = dic

    def set_Cter_UR_fragment(self):
        dic = {'name': "C-ter_UR"}
        Cter = self.tmhmm.fragments[-1]
        dic["end"] = Cter.end
        if Cter.start < self.annotated_domains_fragments[-1]["end"]:
            dic["start"] = self.annotated_domains_fragments[-1]["end"]
        else:
            dic["start"] = Cter.start
        dic["seq"] = self.fasta.get_subsequence(dic["start"], dic["end"])
        self.Cter_UR_fragment = dic

    def set_helix_fragments(self):
        list_fragments = []
        helixes = [f for f in self.tmhmm.fragments if f.cellular_location == "TMhelix"]
        helix_number = 1
        for h in helixes:
            dic = {'name': "TMhelix_" + str(helix_number), 'start': h.start,
                   'end': h.end}
            helix_number += 1
            dic["seq"] = self.fasta.get_subsequence(dic["start"], dic["end"])
            list_fragments.append(dic)
        self.helix_fragments = list_fragments

    def set_loop_fragments(self):
        list_fragments = []
        loops = [f for f in self.tmhmm.fragments[1:-1] if f.cellular_location != "TMhelix"]
        count_inside = 0
        count_outside = 0
        for l in loops:
            dic = {}
            if l.cellular_location == "inside":
                count_inside += 1
                dic["name"] = "inside_loop_" + str(count_inside)
                dic["start"] = l.start
                dic["end"] = l.end
            elif l.cellular_location == "outside":
                count_outside += 1
                dic["name"] = "outside_loop_"+str(count_outside)
                dic["start"] = l.start
                dic["end"] = l.end
            dic["seq"] = self.fasta.get_subsequence(dic["start"], dic["end"])
            list_fragments.append(dic)
        self.loop_fragments = list_fragments

    def set_uniprot_xref(self, uColl):
        p_id = self.prot.split("|")[1]
        try:
            uniprot_entry = uColl.get(p_id)
            self.uniprot_xref = uniprot_entry.xref
        except ValueError as ve:
            ve = str(ve)
            if ve != "Error, empty xmlHandler":
                raise Exception()
            self.uniprot_xref = None

    def get_genome_features(self, enaColl, **kwargs):
        # print("GET GENOME FEATURES")
        # print(kwargs)
        if not self.uniprot_xref:
            raise Exception("Get uniprot xref first.")
        if not self.uniprot_xref.get("EMBL"):
            raise Exception("No EMBL xref. Handle this.")

        # For now, take first id
        ena_id = list(self.uniprot_xref["EMBL"].keys())[0]
        if len(ena_id) != 8: 
            ena_id = ena_id[:6]
        return enaColl.get(ena_id, force_reading_cache=True, **kwargs)

    def set_all_genome_features(self, enaColl, **kwargs):
        self.genome_ena_entry = self.get_genome_features(enaColl, type="all_genome", **kwargs)
        self.genome_ena_id = list(self.uniprot_xref["EMBL"].keys())[0]

    def set_neighborhood(self, number_neighbors, enaColl):
        def filter_type(feature, **kwargs):
            type = kwargs.get("type")
            if not type: 
                raise Exception("Give type argument to filter_type()")
            if feature.type == type:
                return True
            return False

        def filter_no_pseudogene(feature):
            if feature.info.get("other_qualifiers") and "pseudo" in feature.info["other_qualifiers"]:
                return False
            return True

        # Check if required attributes exists
        if not hasattr(self, "uniprot_xref") or not self.uniprot_xref:
            raise Exception("Search uniprot xref first. set_uniprot_xref() method.")
        if not hasattr(self, "genome_ena_entry") or not self.genome_ena_entry:
            raise Exception("Search all genome ena entry first. set_all_genome_features() method")

        # Filter CDS
        # print("filter CDS")
        cds = self.genome_ena_entry.filter(filter_type, type="CDS").features

        # Find protein
        ena_protein_ref = self.uniprot_xref["EMBL"][self.genome_ena_id]
        protein_feature = [f for f in cds if ena_protein_ref in f.info.get("protein_id", [])]
        if not protein_feature:
            raise Exception(self.prot, "protein not found in genome")
        if len(protein_feature) > 1:
            raise Exception(self.prot, "More than 1 protein with id has been found. Handle this")
        protein_feature = protein_feature[0]
        protein_index = cds.index(protein_feature)

        # Get neighborhood
        start = protein_index - number_neighbors
        if start < 0 : 
            start = 0
        end = protein_index + number_neighbors + 1
        if end > len(cds) : 
            end = len(cds)    
        neighborhood = cds[start:end]

        # Correct neighborhood if not in same_contig and delete anchor protein
        neighborhood = [n for n in neighborhood if n.source == protein_feature.source]
        #print(neighborhood)
        neighborhood.remove(protein_feature)

        # Create filters for relaunch genome feature getting
        filter_info = {}
        for n in neighborhood:
            if not n.info.get("locus_tag"):
                # Test with protein_id
                if not n.info.get("protein_id"):
                    raise Exception("No locus_tag, no_protein_id for", self.prot, "Handle this.")
                if "protein_id" not in filter_info:    
                    filter_info["protein_id"] = []
                filter_info["protein_id"] += n.info["protein_id"]
            else:
                if "locus_tag" not in filter_info:
                    filter_info["locus_tag"] = []
                filter_info["locus_tag"] += n.info["locus_tag"]
        # Get genome features again, just keep neighbors and this time keep their sequences.
        # print("NEIGHBORHOOD FEATURES")
        self.set_neighborhood_features(enaColl, info_filter=filter_info, keep_sequence=True, type_filter=["CDS"])

        # Filter pseudogenes
        self.neighborhood_ena_entry = self.neighborhood_ena_entry.filter(filter_no_pseudogene)

    def set_neighborhood_features(self, enaColl, **kwargs):
        self.neighborhood_ena_entry = self.get_genome_features(enaColl, type="neighbors", **kwargs)

    def delete_all_genome_features(self):
        self.genome_ena_entry = None

    @property
    def neighbors_sequences_fasta(self):
        fasta_str = ''
        if not hasattr(self, "neighborhood_ena_entry"):
            raise Exception("Compute neighborhood features first")
        for f in self.neighborhood_ena_entry.features:
            if not f.info.get("translation", None):
                raise Exception("Don't have sequence")
            fasta_str += ">" + self.prot + "+" + f.name + "+" + ";".join(f.info.get("protein_id",'')) + " " + ";".join(f.info.get("product", "")) + "\n"
            fasta_str += f.info.get("translation", [])[0] + "\n"
        return fasta_str

    def get_neighborhood_clusters_number(self):
        if not hasattr(self, "neighborhood_clusters"):
            raise Exception("Compute neighborhood clusters first. Method add_neighborhood_clusters() from TopologyContainer")    

        return list(self.neighborhood_clusters.values())


class Domain():
    def __init__(self, name, hits, proteins, taxo):
        self.name = name
        self.hits = hits
        self.proteins = proteins
        self.taxo = taxo
        self.upper_node = None
        self.mean_distance = None


class Taxo():
    def __init__(self, taxid, taxname, taxrank):
        self.taxid = taxid
        self.taxname = taxname
        self.taxrank = taxrank
