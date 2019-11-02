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
    de_bruijn_counter = defaultdict(Counter)
    # You may also want to check the in-degree and out-degree of each node
    # to help you find the beginnning and end of the sequence.
    for read in sequence_reads:
        # Cut the read into k-mers
        kmers = [read[i: i + k] for i in range(len(read) - k)]
        for i in range(len(kmers) - 1):
            pvs_kmer = kmers[i]
            next_kmer = kmers[i + 1]
            de_bruijn_counter[pvs_kmer].update([next_kmer])

    # This line removes the nodes from the DeBruijn Graph that we have not seen enough.
    #de_bruijn_counter = {key: {val for val in de_bruijn_counter[key] if de_bruijn_counter[key][val] > 1}
    #                   for key in de_bruijn_counter}

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
        #print("db_edges[A]:", db_edges[A])
    return db_edges


def condense_edges(db_graph, db_edges):
    """ Condenses de Bruijn graph by collapsing nodes that have one in-edge and one out-edge. """
    found_match = False
    repeat = True
    while repeat:
        print("starting over")
        repeat = False
        for X in db_graph:
            print("X:", X)
            found_match = False
            for Y in db_graph[X]:
                print("Y:", Y)
                print("db_graph[X][Y]", db_graph[X][Y])
                [XY_edge, XY_cov] = db_edges[X][Y]
                if Y not in db_graph:
                    pass
                elif (len(db_graph[Y]) == 1):
                    for Z in db_graph[Y]:
                        print("Z:", Z)
                        [YZ_edge, YZ_cov] = db_edges[Y][Z]
                        collapse_node(db_graph, db_edges, X, Y, Z)
                        print("collapsed node")
                        found_match = True
                        repeat = True # loop through entire graph again since we changed its size
                        break
            if found_match:
                break
        #if found_match:
        #    break
        #break


def collapse_node(db_graph, db_edges, X, Y, Z):
    """ Given a linear region of the graph where X --XY--> Y --YZ--> Z we collapse the Y node to get X --XYZ--> Z. """
    print("db_edges[X]", db_edges[X])
    print("db_graph[X]", db_graph[X])
    
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
    
    print("db_edges[X]", db_edges[X])
    print("db_graph[X]", db_graph[X])


def condense_db_graph(db_graph):
    """ Condenses de Bruijn graph by collapsing nodes that have one in-edge and one out-edge """
    found_match = False
    repeat = True
    while repeat:
        repeat = False
        print("enterring for loop")
        for parent in db_graph:
            found_match = False
            X = parent
            print("X:", X)
            if (len(db_graph[parent]) == 1):  # if only one child node is present
                for child in db_graph[parent]:  # get the name & coverage of the child node
                    Y = child
                    print("Y:", Y)
                    XY_cov = db_graph[parent][child]
                    if child not in db_graph:  # don't search for grandchildren keys if there aren't any grandchildren
                        pass
                    elif (len(db_graph[child]) == 1):  # if the child has only one child node (grandchild)
                        for grandchild in db_graph[child]:
                            Z = grandchild
                            print("Z:", Z)
                            YZ_cov = db_graph[child][grandchild]
                            collapse_node_og(db_graph, X, Y, Z, XY_cov, YZ_cov)  # collapse the linear region
                            print("collapsed node")
                            found_match = True
                            repeat = True # loop through entire graph again since we changed its size
                    #if found_match:
                    #    break
            if found_match:
                break

    return


def collapse_node_og(db_graph, X, Y, Z, XY_cov, YZ_cov):
    """ Given a linear region of the graph where X --> Y --> Z we collapse the Y node to get XY --> YZ. """
    # new XY and YZ strings
    XY = merge_strings(X, Y)
    YZ = merge_strings(Y, Z)
    #print("XY:", XY)

    # use weighted avg to find new coverage
    new_cov = (XY_cov*len(XY) + YZ_cov*len(YZ)) / (len(XY) + len(YZ))

    # change all values of X in the dict to be XY
    for key in db_graph:
        if X in db_graph[key]:
            #print("key:", key)
            #print("db_graph[key]:", db_graph[key])
            for entry in db_graph[key]:
                if entry == X:
                    db_graph[key][XY] = db_graph[key][X]  # add new XY with same counter value
                    #print("db_graph[key][XY]:", db_graph[key][XY])
                    del db_graph[key][X]  # delete pointer of key (X's parent) to X
            #print("db_graph[key]:", db_graph[key])
            
    # store the children of node Z
    if Z in db_graph:
        temp = db_graph[Z]
    else:
        temp = Counter({})  # set to empty counter if Z has no children (end of graph)
    
    # delete nodes X, Y, and Z
    if X in db_graph:
        del db_graph[X]
    if Y in db_graph:
        del db_graph[Y]
    if Z in db_graph:
        del db_graph[Z]
    
    # create node XY that points to YZ with coverage of new_cov
    db_graph[XY] = Counter({YZ: new_cov})
    #print("db_graph[XY]:", db_graph[XY])
    
    # have YZ node point to temp
    db_graph[YZ] = temp
    #print("db_graph[YZ]:", db_graph[YZ])
    return


def merge_strings(A, B):
    """ Merges strings A and B into AB assuming that the right end of A matches the left end of B. """
    max_overlap = min(len(A), len(B))
    
    for i in range(max_overlap, 0, -1):
        A_end = A[-i:]  # get last i characters of A
        B_beg = B[:i]  # get first i characters of B
        #print("A_end:", A_end)
        #print("B_beg:", B_beg)
        if A_end == B_beg:
            A_beg = A[:(len(A)-i)]  # get first len(A)-i characters of A
            B_end = B[-(len(B)-i):]  # get last len(B)-i characters of B
            #print("A_beg:", A_beg)
            #print("B_end:", B_end)
            AB = A_beg + A_end + B_end  # collapse the region of overlap
            break
    return AB


def plot_db_graph(db_graph):
    """ Plots the De Bruijn Graph into a dot file using pygraphviz. """
    A = pgv.AGraph()
    for key in db_graph:
        for cntr in db_graph[key]:
            A.add_edge(key, cntr, label=("cov = " + str(round(db_graph[key][cntr], 3)) + ", len = " + str(len(merge_strings(key, cntr)))))
    A.node_attr.update(label=0, fontsize=0)
    A.write("test.dot")  # use "dot -Tpng test.dot > test.png" to convert to png
    #A.draw('test.png')


if __name__ == "__main__":

    #A = "ABCDE"
    #B = "CDEF"
    #print(merge_strings(A,B))
    #reads_fn = "MG1655-K12.first1K.fasta"
    reads_fn = "s_6.first1000.fastq"
    reads = read_assembly_reads(reads_fn)
    db_graph = simple_de_bruijn(reads, 55)
    db_edges = build_edges(db_graph)
    #test = "CCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGACGCGT"
    #print("test case:", db_edges[test])
    condense_edges(db_graph, db_edges)
    #condense_db_graph(db_graph)
    plot_db_graph(db_graph)

    #output_fn = "fastq_reads.txt"
    #with open(output_fn, 'w') as output_file:
    #    output_file.write('>' + reads_fn + '\n')
    #    output_file.write('>READS\n')
    #    output_file.write('\n'.join(reads))
