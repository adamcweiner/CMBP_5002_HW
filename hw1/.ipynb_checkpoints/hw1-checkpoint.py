import numpy as np


def superstring(s, t):
    """ Given: two strings s and t. 
        Return: shortest common supersequence of s and t. """
    s = "$" + s
    t = "$" + t
    s_list = list(s)
    t_list = list(t)
    grid = np.zeros((len(s), len(t)))
    grid[0, :] = np.arange(len(t))
    grid[:, 0] = np.arange(len(s))

    for i in range(1, len(s)):
        for j in range(1, len(t)):
            if s[i] == t[j]:
                grid[i, j] = 1 + grid[i-1, j-1]  # add one & move diag if they match
            else:
                grid[i, j] = 1 + min(grid[i, j-1], grid[i-1, j])  # add one and move either up or down if they don't match
    
    # now backtrack through the grid and print the shortest common supersequence
    i = len(s) - 1
    j = len(t) - 1
    out = ""  # output string
    while i > 0 and j > 0:
        if s[i] == t[j]:  # if chars match
            out = s[i] + out  # put string in beginning of out
            i -= 1  # decrease both i and j
            j -= 1
        elif grid[i, j-1] > grid[i-1, j]:  # case where it's better to go left than up
            out = s[i] + out
            i -= 1  # only decrease i
        else:  # case where it's better to go up than left
            out = t[j] + out
            j -= 1
    
    if i > 0:  # if you are to the right of top left corner
        out = s[1:i+1] + out  # print the first i letters of s
    
    if j > 0:  # if you are underneath the top left corner
        out = t[1:j+j] + out  # print the first j letters of t

    return out


def overlap_alignment(s, t):
    """ Find the optimal local alignment of a suffix of s with a prefix of t. Return the max alignment score and the alignment that achieves this score. """
    # store lengths as constants
    m = len(s)
    n = len(t)

    mu = -2  # penalty for mismatches
    sigma = -2  # penalty (linear) for indels
    tau = 1  # reward for matches
    
    # holds the number of substitutions between versions of the two strings
    edit_mat = np.zeros((m+1, n+1), dtype=int)
    backtrack = edit_mat.copy()
    
    # give indel penalty for first row but leave first col as 0
    for j in range(1, n+1):
        edit_mat[0, j] = sigma*j
        backtrack[0, j] = 2  # 2 for going right

    # loop through and find number of edit distances
    for i in range(1, m+1):
        for j in range(1, n+1):
            if s[i-1] == t[j-1]:
                edit_mat[i, j] = edit_mat[i-1, j-1] + tau
                backtrack[i, j] = 4  # 4 for a match
            # if the two characters don't match at the current position then take the min value of a SNP, insertion, or deletion and add 1
            else:
                edit_mat[i, j] = max(edit_mat[i-1, j] + sigma, edit_mat[i, j-1] + sigma, edit_mat[i-1, j-1] + mu)
                backtrack[i, j] = np.argmax([edit_mat[i-1, j] + sigma, edit_mat[i, j-1] + sigma, edit_mat[i-1, j-1] + mu]) + 1

    i = m
    j = np.argmax(edit_mat[m, :])  # start in column with highest score in bottom row
    alignment_score = edit_mat[i, j]
    
    s_aligned = ""
    t_aligned = ""
    while i > 0 or j > 0:
        if backtrack[i, j] == 0: # break out of while-loop if we've reached the left edge
            break
        elif backtrack[i, j] == 3 or backtrack[i, j] == 4:
            s_aligned = s[i-1] + s_aligned # add current letters to beginning of aligned strings
            t_aligned = t[j-1] + t_aligned
            i -= 1 # move diagonally up and left
            j -= 1
        elif backtrack[i, j] == 1: # going up
            s_aligned = s[i-1] + s_aligned # add current letter from s
            t_aligned = "-" + t_aligned # add dash to represent gap in t
            i -= 1 # only move up
        elif backtrack[i, j] == 2:  # going left
            s_aligned = "-" + s_aligned  # add dash to represent gap in s
            t_aligned = t[j-1] + t_aligned  # add current letter from t
            j -= 1  # only move left
        else:
            print("invalid entry found in backtrack")

    return s_aligned, t_aligned, alignment_score
    

def load_matrix(matrix_filename):
    """ Imports the blosum62.txt file as a dict/matrix. """
    with open(matrix_filename) as matrix_file:
        matrix = matrix_file.read()
    lines = matrix.strip().split('\n')

    header = lines.pop(0)
    columns = header.split()
    matrix = {}

    for row in lines:
        entries = row.split()
        row_name = entries.pop(0)
        matrix[row_name] = {}

        if len(entries) != len(columns):
            raise Exception('Improper entry number in row')
        for column_name in columns:
            matrix[row_name][column_name] = entries.pop(0)

    return matrix


def get_blos_score(matrix, a, b):
    """ Given a BLOSUM62 matrix and two characters (a & b), find the score of this relationship in the matrix. """
    a = a.upper()
    b = b.upper()

    if a not in matrix or b not in matrix[a]:
        raise InvalidPairException('[%s, %s]' % (a, b))
    return int(matrix[a][b])


def protein_alignment(s, t):
    """ Finds BLOSUM62 alignment between two protein sequences s and t. This alignment uses an affine (non-linear) gap penalty. """
    gap_ext = -1
    gap_open = -11
    blos = load_matrix("blosum62.txt")
    m = len(s)
    n = len(t)

    # holds the number of substitutions between versions of the two strings
    edit_mat = np.zeros((m+1, n+1), dtype=int)
    backtrack = edit_mat.copy()
    
    # fill left column
    for i in range(1, m+1):
        backtrack[i, 0] = 1  # move vertically throughout this column
        # open gap at i=1
        if i == 1:
            edit_mat[i, 0] = edit_mat[i-1, 0] + gap_open
        # extend gap for all other rows
        else:
            edit_mat[i, 0] = edit_mat[i-1, 0] + gap_ext
    
    # fill top row
    for j in range(1, n+1):
        backtrack[0, j] = 2  # move horizontally throughout this row
        # open gap at j=1
        if j == 1:
            edit_mat[0, j] = edit_mat[0, j-1] + gap_open
        # extend gap for all other columns
        else:
            edit_mat[0, j] = edit_mat[0, j-1] + gap_ext

    for i in range(1, m+1):
        for j in range(1, n+1):
            # use gap_ext if you moved down in the previous position
            if backtrack[i-1, j] == 1:
                down_pen = gap_ext
            # use gap_open if you aren't continuing and existing downward gap
            else:
                down_pen = gap_open
            
            # repeat logic for gaps that are moving right
            if backtrack[i, j-1] == 2:
                right_pen = gap_ext
            else:
                right_pen = gap_open

            mu = get_blos_score(blos, s[i-1], t[j-1])  # find blosum62 score for moving diagonally

            edit_mat[i, j] = max(edit_mat[i-1, j] + down_pen, edit_mat[i, j-1] + right_pen, edit_mat[i-1, j-1] + mu)
            backtrack[i, j] = np.argmax([edit_mat[i-1, j] + down_pen, edit_mat[i, j-1] + right_pen, edit_mat[i-1, j-1] + mu]) + 1

    
    s_aligned = ""
    t_aligned = ""
    i, j = m, n  # start at bottom-right corner
    alignment_score = edit_mat[i, j]  # score in this corner is the maximum alignment score

    while i > 0 and j > 0:
        if backtrack[i, j] == 0: # break out of while-loop if we've reached the top-left corner
            break
        elif backtrack[i, j] == 3:  # move diagonally
            s_aligned = s[i-1] + s_aligned # add current letters to beginning of aligned strings
            t_aligned = t[j-1] + t_aligned
            i -= 1 # move diagonally up and left
            j -= 1
        elif backtrack[i, j] == 1: # going up
            s_aligned = s[i-1] + s_aligned # add current letter from s
            t_aligned = "-" + t_aligned # add dash to represent gap in t
            i -= 1 # only move up
        elif backtrack[i, j] == 2:  # going left
            s_aligned = "-" + s_aligned  # add dash to represent gap in s
            t_aligned = t[j-1] + t_aligned  # add current letter from t
            j -= 1  # only move left
        else:
            print("invalid entry found in backtrack")

    return s_aligned, t_aligned, alignment_score
    

def import_fasta(file_name):
    """ Loads fasta file into a list where first element is s and second element is t according to problem statements for problems 2 and 3. """
    file = open(file_name)
    out = []
    for num, line in enumerate(file, 1):
        if line[0] != '>':
            out.append(line.rstrip('\n'))
    assert(len(out) == 2)
    return out


def main():
    """ Perform main routine. """
    # problem 1
    print("starting problem 1")
    str1 = "ACGTC"
    str2 = "ATAT"
    prob1 = superstring(str1, str2)
    print("shortest common supersequence:", prob1)
    
    # problem 2
    print("starting problem 2")
    [s, t] = import_fasta("prob2_ex.fasta")
    s_suffix, t_prefix, score = overlap_alignment(s, t)
    print(s_suffix)
    print(t_prefix)
    print(score)

    # problem 3
    print("starting problem 3")
    [s, t] = import_fasta("prob3_ex.fasta")
    s_aligned, t_aligned, score = protein_alignment(s, t)
    print(s_aligned)
    print(t_aligned)
    print("score:", score)


if __name__ == "__main__":
    main()
