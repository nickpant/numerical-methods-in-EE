# -*- coding: utf-8 -*-
"""
@author: Nick Pant
Affiliation: McGill University
Class: ECSE 543 - Numerical Methods in Electrical Engineering

Assignment 1, Problem 2 part a

"""

import os

# create an Nx2N matrix
def construct_matrix(N):

    return [[0]*2*N for x in range(N)]

# returns the node number for index (i,j), ex:
# 0 1 2 GROUND
# 3 4 5 6
# 7 8 9 10
def index_to_node(i, j, J):


    if (i == 0 and j < J-1):

        return i * J + j

    elif (i == 0 and j == J-1):

        return 'GROUND'

    else:

        return i * J + j - 1



# return the branch number than connects nodes N1 and N2
# A is the incidence matrix
def nodes_to_branch(N1, N2, A):

    for i in range(len(A[0])):

        if (A[N1][i] != 0 and A[N2][i] != 0):

            return i

    return 'NaN' # if no branch exists, return NaN


# solve the matrix equation L y = b
def solve(A_passed, b_passed):

    if isinstance(b_passed, list) == False:

        return b_passed / A_passed[0][0]

    n = len(A_passed)

    # perform the elimination step
    for j in range(n):

        if (A_passed[j][j] <= 0):

            return -1

        A_passed[j][j] = (A_passed[j][j])**(0.5) # square root of A[j][j]

        b_passed[j] = b_passed[j]/A_passed[j][j]

        for i in range(j):
            A_passed[i][j] = 0

        for i in range(j+1, n):

            A_passed[i][j] = A_passed[i][j] / A_passed[j][j]

            b_passed[i] = b_passed[i] - A_passed[i][j] * b_passed[j]

            for k in range(j+1, i+1):

                A_passed[i][k] = A_passed[i][k] - A_passed[i][j] * A_passed[k][j]

    return solve_x(A_passed, b_passed)

# solve the matrix equation L^T x = y
def solve_x(L_passed, y_passed):

    if isinstance(y_passed, list) == False:
        return y_passed/L_passed[0][0]

    n = len(L_passed)

    for j in range(n-1, -1, -1):

        y_passed[j] = y_passed[j]/L_passed[j][j]

        for i in range(j):

            y_passed[i] = y_passed[i] - L_passed[j][i] * y_passed[j]

    return y_passed

# print any m x n matrix as long as m != 1
def print_matrix(M):

    for row in range(0, len(M)):

        temp = ''

        for col in range(len(M[row])):

            temp += str(round(M[row][col], 4)) + '\t'

        print(temp)

# returns the i-th column of Matrix M
def column(M, i):

    return [row[i] for row in M]

# transpose a matrix M to M^T
def transpose(M):

    if len(M) == 1:
        return M

    # if m != n, EXIT -> non-square matrices cannot be transposed
    if len(M) != len(M[0]):
        return -1

    transposed = [[0]*len(M) for i in range(len(M))]

    for i in range(len(M)):

        for j in range(len(M)):

            transposed[i][j] = M[j][i]

    return transposed

# M1 + M2, m = 1 implies add, m = -1 implies substract
def add(M1, M2, m):

    if isinstance(M1, list):

        if isinstance(M1[0], list):

            summation = [[0]*len(M1[i]) for i in range(len(M1))]

            for row in range(len(M1)):

                for col in range(len(M1[row])):

                    summation[row][col] = M1[row][col] + m * M2[row][col]

            return summation

        else:

            summation = [0] * len(M1)

            for row in range(len(M1)):

                summation[row] = M1[row] + m * M2[row]

            return summation
    else:

        return 'NAV' # if neither M1 nor M2 are matrices/vectors, return not a vector


# dot a vector v with another vector v
def dot_vv(v1, v2):

    prod = 0

    for row in range(0, len(v1)):

        prod = prod + v1[row]*v2[row]

    return prod


# dot a matrix M with a vector v
def dot_Mv(A_passed, x_passed):

    if (isinstance(A_passed[0], list) == False or len(A_passed) == 1):

        return dot_vv(A_passed[0], x_passed)

    prod = [0] * len(A_passed)

    for row in range(len(A_passed)):

        for col in range(len(A_passed[row])):

            prod[row] = prod[row] + A_passed[row][col] * x_passed[col]

    return prod


# multiply two matrices together by dot product
def dot_MM(M1, M2):

    m = len(M1)
    n = len(M1[0])
    p = len(M2)
    q = len(M2[0])

    # if the inner dimensions don't match, matrix multiplication not possible
    if n != p:

        return -1

    # create an mxq matrix for the product
    prod = [[0]*q for row in range(m)]

    for row in range(m):

        for col in range(q):

            prod[row][col] = dot_vv(M1[row][:], column([row[:] for row in M2], col))

    return prod


def dot_MMT(M1, M2):

    m = len(M1)
    n = len(M1[0])
    p = len(M2[0])
    q = len(M2)

    if n != p:

        return -1

    prod = [[0]*q for row in range(m)]

    for row in range(m):

        for col in range(q):

            prod[row][col] = dot_vv(M1[row][:], M2[col][:])

    return prod


# read the circuit file
def read_circuit(file):

    # read from fine a list of Jk, Rk, Ek and a reduced incidence matrix

    with open(file) as f:

        line = f.readline().replace('\n', '').split('\t')

        # read the incidence matrix A

        if line[0] == 'A':
            # line[2] contains row length and line[4] col length
            nrows = int(line[1])
            ncols = int(line[2])

            A = [[0]*ncols for r in range(nrows)]

            J = [0]*len(A[0])

            E = [0]*len(A[0])

            R = [0]*len(A[0])

        for i in range(nrows):

            row = f.readline().replace('\n', '').split('\t')

            for j in range(ncols):

                A[i][j] = int(row[j])

        f.readline()

        # read the R, E, J vectors

        for i in range(len(A[0])):

            row = f.readline().replace('\n', '').split('\t')

            R[i] = int(row[0])

            E[i] = int(row[1])

            J[i] = int(row[2])


        # convert the resistance vector R into an admittance matrix Y

        Y = [[0]*len(R) for y in range(len(R))]

        for i in range(len(Y)):

            Y[i][i] = 1/R[i]

            for j in range(len(Y[i])):

                if i != j:

                    Y[i][j] = 0


        return A, J, E, Y # note that J and E should be vertical vectors
        # but in practice, they are represented as horizontal vectors

def solve_full(A, J, E, Y):

    AYAT = dot_MMT(dot_MM(A,Y), A)

    YE = dot_Mv(Y, E)

    J_YE = add(J, YE, -1)

    b = dot_Mv(A, J_YE)

    # now we solve AYAT v = b

    x = solve(AYAT, b)

    return x


# apply voltage to the bottom right node from the top left node, which is grounded
def apply_voltage(A, J, E, Y, temp_node):

    new_A = [[0] * (len(A[0])+1) for row in range(len(A))]

    new_Y = [[0] * (len(Y[0])+1) for row in range(len(Y[0]) + 1)]

    # set the new A matrix with values of A

    for i in range(len(A)):

        for j in range(len(A[0])):

            new_A[i][j] = A[i][j]

    # new branch connects to temp node with an incoming current

    new_A[temp_node][-1] = -1

    # set the new Y matrix with values of Y

    for i in range(len(Y)):

        new_Y[i][i] = Y[i][i]

    # the new branch has a resistance of 1 kOhm

    new_Y[-1][-1] = 100

    # set the new E and J matrix with their respective values

    new_E = [0] * (len(E) + 1)

    new_J = [0] * (len(J) + 1)

    for i in range(len(E)):

        new_E[i] = E[i]

        new_J[i] = J[i]

    new_E[-1] = 10

    return new_A, new_J, new_E, new_Y


def main():

    # change directory
    os.chdir('C:/Users/aagni/Google Drive/McGill/U3/ECSE 543/Assignment 1/1.2')

    A, J, E, Y = read_circuit('1.2_test.txt') # read the circuit from file

    N = 2

    N_ROWS = N

    N_COLS = 2*N

    temp_node = index_to_node(N_ROWS - 1, 0, N_COLS) # bottom right node

    new_A, new_J, new_E, new_Y = apply_voltage(A, J, E, Y, temp_node)

    v = solve_full(new_A, new_J, new_E, new_Y)

    adj_node_right = index_to_node(N_ROWS - 1, 1, N_COLS) # adjacent to right

    adj_node_up = index_to_node(N_ROWS - 2, 0, N_COLS) # adajacent to top

    branch_right = nodes_to_branch(temp_node, adj_node_right, new_A)

    branch_up = nodes_to_branch(temp_node, adj_node_up, new_A)

    print(v)

    i_1 = (v[temp_node] - v[adj_node_up]) * new_Y[branch_up][branch_up]

    i_2 = (v[temp_node] - v[adj_node_right]) * new_Y[branch_right][branch_right]

    R_eq = v[temp_node] / (i_1 + i_2)

    print(str(round(R_eq, 5)) + ' ohms')


if __name__ == '__main__':
    main()