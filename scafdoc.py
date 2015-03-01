import MySQLdb
from igraph import *
from igraph import _igraph
import Queue

global n
max_round = 4

db = MySQLdb.connect(host="carnot.ncats.nih.gov", # your host, usually localhost
                     user="chembl_18", # your username
                      passwd="chembl_18", # your password
                      db="chembl_18") # name o
cur = db.cursor()

def node_exists(g,name):
    try:
        g.vs.find(Name=name)
        return True
    except:
        return False
    
def pg(g):
    color_dict = {"doc": "red", "fragment": "black"}
    g.vs["color"] = [color_dict[klass] for klass in g.vs["klass"]]
    g.vs['label'] = g.vs['Name']
    layout = g.layout("fr")
    
    plot(g, layout = layout, vertex_size=5)

def get_docs_for_fragid(fid, min_year = 1997, max_year = 2015, cur = None):
    sql = """select distinct pubmed_id, title, abstract, year
  from ncgc_fragment_class x, ncgc_fragment_instances a, activities b, docs c
where x.class_id = %d
and x.class_id = a.class_id
and a.molregno = b.molregno
and c.doc_id = b.doc_id
and year >=  %d and year <= %d
and pubmed_id is not null""" % (fid, min_year, max_year)

    cur.execute(sql)
    return cur.fetchall()

def get_frag_for_docid(pmid, max_inst = 300, min_inst = 20, min_complex = 300,  cur = None):
    sql = """select distinct x.class_id, x.smiles, acount, complexity, symmetry, instances
  from ncgc_fragment_class x, ncgc_fragment_instances a, activities b, docs c
where c.pubmed_id = %d
and c.doc_id = b.doc_id
and a.molregno = b.molregno
and x.class_id = a.class_id and instances < %d and instances > %d and complexity > %d""" % (pmid, max_inst, min_inst, min_complex)
    cur.execute(sql)
    return cur.fetchall()

def handle_doc_nodes(g, node, new_nodes):
    assert node['klass'] == 'fragment'
    node_name = node['Name']
    for pmid, title, abstract, year in new_nodes:
        try:
            eid = g.get_eid( g.vs.find(Name=node_name).index, g.vs.find(Name=pmid), directed=False, error=True )
        except (ValueError,_igraph.InternalError) as e:
            if not node_exists(g, pmid): g.add_vertex(Name=pmid, pubmed_id = pmid, title=title, year=year, klass='doc')
            g.add_edge(g.vs.find(Name=pmid), g.vs.find(Name=node_name))
    return (g, [x[0] for x in new_nodes])

def handle_fragment_nodes(g, node, new_nodes):
    assert node['klass'] == 'doc'
    node_name = node['Name']
    for class_id, smi, acount, complexity, sym, n in new_nodes:
        try:
            eid = g.get_eid( g.vs.find(Name=node_name).index, g.vs.find(Name=class_id), directed=False, error=True  )
        except (ValueError,_igraph.InternalError) as e:
            if not node_exists(g, class_id): g.add_vertex(Name=class_id, fid=class_id, smi=smi, acount=acount, complexity=complexity, sym=sym, klass='fragment')
            g.add_edge(g.vs.find(Name=node_name), g.vs.find(Name=class_id))
    return (g, [x[0] for x in new_nodes])
        
def extend_graph(g, node, cur):
    node_name = node['Name']
    node_type = node['klass']
    if node_type == 'doc':
        new_nodes = get_frag_for_docid(node_name, cur=cur)
        new_nodes = list(set(new_nodes))
        g, newids = handle_fragment_nodes(g, node, new_nodes)
    else:
        new_nodes = get_docs_for_fragid(node_name, cur=cur)
        new_nodes = list(set(new_nodes))
        
        g, newids = handle_doc_nodes(g, node, new_nodes)
    return g, newids
        
    
if __name__ == '__main__':

    g = Graph()
    
    seed_doc = 6687479
    ##seed_doc = 20871596 ## JQ1 paper
    ##seed_doc = 12477366
    ##seed_doc = 21568322 ## BET family bromodomains
    g.add_vertex(Name=seed_doc, pubmed_id = seed_doc, title = None, klass='doc')
    g, newids = extend_graph(g, g.vs.find(Name=seed_doc), cur=cur)

    q = newids[:]

    n = 0
    while n < max_round:
        nq = []
        qlen = len(q)
        print 'Round %d, queue size is %d, graph [V=%d] [E=%d]' % (n, qlen, g.vcount(), g.ecount())

        i = 1
        for node_name in q:
            node = g.vs.find(Name=node_name)
            g, newids = extend_graph(g, node, cur=cur)
            sys.stdout.write('\r [%d/%d] Extending from %d [%s] and got %d new nodes    ' % (i, qlen, node_name, node['klass'], len(newids)))
            sys.stdout.flush()
            nq.extend(newids)
            i += 1
        print
        n += 1
        q = nq[:]

    print g.vcount(), g.ecount()
    if g.vcount() > 1: g.save('%d-%d.gml' % (max_round,seed_doc),format='gml')
