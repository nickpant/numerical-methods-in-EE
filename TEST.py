# -*- coding: utf-8 -*-
"""
Created on Sun Oct 15 17:19:33 2017

@author: aagni
"""

import numpy as np

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

# print any m x n matrix as long as m != 1
def print_matrix(M):

    for row in range(0, len(M)):

        temp = ''

        for col in range(len(M[row])):

            temp += str(round(M[row][col], 4)) + '\t'

        print(temp)

def main():

    os.chdir('C:/Users/aagni/Google Drive/McGill/U3/ECSE 543/Assignment 1/1.1')

    A, J, E, Y = read_circuit('1.1_test5.txt')

    print_matrix(A)

    print()

    print(J)

    print()

    print(E)

    print()

    print_matrix(Y)

    print()

    AYAT = np.dot(np.dot(np.matrix(A), np.matrix(Y)), np.transpose(np.matrix(A)))

    L = np.linalg.cholesky(AYAT)


    print(L)

    '''

    print()


   #s print_matrix(AYAT)

    print(np.vstack(J) - np.dot(np.matrix(Y), np.vstack(E)))

    print()

    b = np.dot(np.matrix(A), np.vstack(np.array(J)) - np.dot(np.matrix(Y),np.vstack(np.array(E))))

    print()
    print(b)

    print()

    print(np.linalg.solve(AYAT,b))

'''

if __name__ == '__main__':

    main()