""" Perform analysis for HW2 assignment. """
from os.path import join
import sys
import pygraphviz as pgv
from collections import defaultdict, Counter


def read_fasta(read_fn):
    f = open(read_fn, 'r')
    first_line = True
    all_reads = []
    count = 0
    for line in f:
        count += 1
        if count % 1000 == 0:
            print(count, " reads done.")
        if first_line:
            first_line = False
            continue  # We skip the first line, since it
            # only contains the name of the chromosome the reads
            # came from.
        line = line.strip()
        all_reads.append(line)  # forward direction
        all_reads.append(''.join(reversed(line)))  # backward direction
    return all_reads


def read_fastq(read_fn):
    f = open(read_fn, 'r')
    all_reads = []
    next_line = False
    count = 0
    for line in f:
        if count % 1000 == 0:
            print(count, " reads done.")
        line = line.strip()
        # only append if its a true read
        if next_line and line[0] != "@":
            count += 1
            all_reads.append(line)  # forward direction
            all_reads.append(''.join(reversed(line)))  # backward direction
        # true read always follows the line that starts with '@' character
        if line[0] == "@":
            next_line = True
        else:
            next_line = False
    return all_reads


def read_assembly_reads(read_fn):
    if read_fn[-1] == "q":
        print("importing fastq")
        reads = read_fastq(read_fn)
    else:
        print("importing fasta")
        reads = read_fasta(read_fn)
    return reads


def simple_de_bruijn(sequence_reads, k):
    """
    Creates A simple DeBruijn Graph with nodes that correspond to k-mers of size k.
    :param sequence_reads: A list of reads from the genome
    :param k: The length of the k-mers that are used as nodes of the DeBruijn graph
    :return: A DeBruijn graph where the keys are k-mers and the values are the set
                of k-mers that
    """
    read_count = 0
    de_bruijn_counter = defaultdict(Counter)
    # loop through each read
    for read in sequence_reads:
        kmers = [read[i: i + k] for i in range(len(read) - k)]
        # loop through each kmer in this read's spectrum
        for i in range(len(kmers) - 1):
            pvs_kmer = kmers[i]
            next_kmer = kmers[i + 1]      
            de_bruijn_counter[pvs_kmer].update([next_kmer])  # connect the kmers

    # This line removes the empty nodes from the DeBruijn graph
    de_bruijn_graph = {key: de_bruijn_counter[key] for key in de_bruijn_counter if de_bruijn_counter[key]}
    return de_bruijn_graph


def build_edges(db_graph):
    """ Build a dict of dict of Counters where the counter is the edge and the coverage between the two nodes. For example {A: {B: [AB_edge, AB_coverage]}}. """
    db_edges = defaultdict()
    for A in db_graph:
        i = 0
        for B in db_graph[A]:
            AB_edge = merge_strings(A, B)
            AB_coverage = db_graph[A][B]
            if i == 0:
                db_edges[A] = {B: [AB_edge, AB_coverage]}
            else:
                db_edges[A].update({B: [AB_edge, AB_coverage]})
            i += 1
    return db_edges


def condense_graph(db_graph, db_edges):
    """ Condenses de Bruijn graph by collapsing nodes that have one in-edge and one out-edge. """
    found_match = False
    repeat = True
    while repeat:
        #print("size of db graph:", len(db_graph))
        repeat = False
        for X in db_graph:
            found_match = False
            for Y in db_graph[X]:
                [XY_edge, XY_cov] = db_edges[X][Y]
                if Y not in db_graph:
                    pass
                elif (len(db_graph[Y]) == 1):
                    for Z in db_graph[Y]:
                        [YZ_edge, YZ_cov] = db_edges[Y][Z]
                        collapse_node(db_graph, db_edges, X, Y, Z)
                        found_match = True
                        repeat = True # loop through entire graph again since we changed its size
                        break
            if found_match:
                break


def collapse_node(db_graph, db_edges, X, Y, Z):
    """ Given a linear region of the graph where X --XY--> Y --YZ--> Z we collapse the Y node to get X --XYZ--> Z. """
    # find strings and coverages for original edges
    [XY_edge, XY_cov] = db_edges[X][Y]
    [YZ_edge, YZ_cov] = db_edges[Y][Z]

    XYZ_edge = merge_strings(XY_edge, YZ_edge)
    XYZ_cov = (XY_cov*len(XY_edge) + YZ_cov*len(YZ_edge)) / (len(XY_edge) + len(YZ_edge))
    
    # find edge leading out of X that points to Y
    new_children = Counter()
    for child in db_graph[X]:
        if child != Y:
            new_children.update(Counter({child: db_graph[X][child]}))  # add children that aren't Y
    
    # add Z to new_children
    new_children.update(Counter({Z: XYZ_cov}))
    
    # have X point to new_children
    db_graph[X] = new_children
    
    # delete all traces of Y from the graph
    del db_edges[X][Y]
    del db_edges[Y]
    del db_graph[Y]
    
    # create db_edges[X][Z]
    db_edges[X].update({Z: [XYZ_edge, XYZ_cov]})


def merge_strings(A, B):
    """ Merges strings A and B into AB assuming that the right end of A matches the left end of B. """
    max_overlap = min(len(A), len(B))
    
    for i in range(max_overlap, 0, -1):
        A_end = A[-i:]  # get last i characters of A
        B_beg = B[:i]  # get first i characters of B
        if A_end == B_beg:
            A_beg = A[:(len(A)-i)]  # get first len(A)-i characters of A
            B_end = B[-(len(B)-i):]  # get last len(B)-i characters of B
            AB = A_beg + A_end + B_end  # collapse the region of overlap
            break
    return AB


def plot_db_graph(db_edges, output_file, min_cov=0, min_len=0):
    """ Plots the De Bruijn Graph into a dot file using pygraphviz. Can specify miminal coverage and edge length for graph simplification. """
    A = pgv.AGraph()
    for X in db_edges:
        for Y in db_edges[X]:
            [XY_edge, XY_cov] = db_edges[X][Y]
            if (len(XY_edge) >= min_len and XY_cov >= min_cov):
                A.add_edge(X, Y, label=("cov = " + str(round(XY_cov, 2)) + ", len = " + str(len(XY_edge))))
    A.node_attr.update(label=0, fontsize=0)
    A.write(output_file)

    
def plot_db_tip_removal(db_edges, output_file, tip_cov=0, tip_len=0):
    """ Plots de Bruijn graph into a dot file using pygraphvis. Any tip that falls below the specified tip_cov or tip_len will not be included in the graph. """
    A = pgv.AGraph()
    for X in db_edges:
        if len(db_edges[X]) > 1:  # use if there's a fork in the graph the edge must pass tip thresholds
            for Y in db_edges[X]:
                [XY_edge, XY_cov] = db_edges[X][Y]
                if (len(XY_edge) >= tip_len and XY_cov >= tip_cov):
                        A.add_edge(X, Y, label=("cov = " + str(round(XY_cov, 2)) + ", len = " + str(len(XY_edge))))
        else:  # plot edges normally if there isn't a fork in the graph
            for Y in db_edges[X]:
                [XY_edge, XY_cov] = db_edges[X][Y]
                A.add_edge(X, Y, label=("cov = " + str(round(XY_cov, 2)) + ", len = " + str(len(XY_edge))))
    A.node_attr.update(label=0, fontsize=0)
    A.write(output_file)


def fasta_edges(db_edges, output_fn, min_cov=0, min_len=0):
    """ Writes a fasta format file of all the edges in the graph that meet the specified coverage & length criteria. """
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + reads_fn + '\n')
        output_file.write('>EDGES,COVERAGE\n')
        for X in db_edges:
            for Y in db_edges[X]:
                [XY_edge, XY_cov] = db_edges[X][Y]
                if (len(XY_edge) >= min_len and XY_cov >= min_cov):
                    output_file.write(str(XY_edge) + ',' + str(round(XY_cov, 2)) + '\n')
    return


def fasta_edges_tip_removal(db_edges, output_fn, tip_cov=0, tip_len=0):
    """ Same as fasta_edges() except this function removes tips from the graph. """
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + reads_fn + '\n')
        output_file.write('>EDGES,COVERAGE\n')
        for X in db_edges:
                if len(db_edges[X]) > 1:  # check tip requirements when there's a branch in the graph
                    for Y in db_edges[X]:
                        [XY_edge, XY_cov] = db_edges[X][Y]
                        if (len(XY_edge) >= tip_len and XY_cov >= tip_cov):
                            output_file.write(str(XY_edge) + ',' + str(round(XY_cov, 2)) + '\n')
                else:  # just plot the edge when there's no branch
                    for Y in db_edges[X]:
                        [XY_edge, XY_cov] = db_edges[X][Y]
                        output_file.write(str(XY_edge) + ',' + str(round(XY_cov, 2)) + '\n')
    return


if __name__ == "__main__":
    # perform analysis for s_6.first1000 file
    reads_fn = "s_6.first1000.fastq"
    reads = read_assembly_reads(reads_fn)
    print("making db_graph")
    db_graph = simple_de_bruijn(reads, 55)
    print("done with db_graph... making db_edges")
    db_edges = build_edges(db_graph)
    print("made db_edges... condensing graph")
    condense_graph(db_graph, db_edges)
    print("condensed graph... making plots")
    plot_db_graph(db_edges, "s6_normal_db.dot")  # use "dot -Tpng s6_normal_db.dot > s6_normal_db.png" to convert to png
    fasta_edges(db_edges, "s6_normal.edges.fasta")
    plot_db_tip_removal(db_edges, "s6_tip_removal.dot", tip_cov=10, tip_len=100)  # use "dot -Tpng s6_tip_removal.dot > s6_tip_removal.png" to convert to png
    fasta_edges_tip_removal(db_edges, "s6_tip_removal.edges.fasta", tip_cov=10, tip_len=100)
    plot_db_graph(db_edges, "s6_high_quality.dot", min_cov=10, min_len=100)  # use "dot -Tpng s6_high_quality.dot > s6_high_quality.png" to convert to png
    fasta_edges(db_edges, "s6_high_quality.edges.fasta", min_cov=10, min_len=100)

