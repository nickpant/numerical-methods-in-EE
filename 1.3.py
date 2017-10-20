# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 17:36:06 2017

@author: Nick Pant
Affiliation: McGill University
Class: ECSE 543 - Numerical Methods in Electrical Engineering

Assignment 1, Problem 3

"""

import os


# return an index on the vector corresponding to its (i, j) position
def index_to_node(i, j, J):

        return i * J + j

# fill the x matrix with boundary condition

def boundary_condition(N, M, h, w, V):

    b = [0]*(N * M)

    for i in range(0, h):

        for j in range(M - 1, M - 1 - w - 1, -1):

            ind = index_to_node(i, j, M)

            b[ind] = V

    return b

# check if the residuals are all below 10^-5
def residuals_satisfied(x, N, M, h, w):

    Rmax = 0

    for i in range(1, N-1):

        for j in range(1, M-1):

            if not ((j > N - 1 - w) and (i < h)):

                ind = index_to_node(i, j, M)

                R = -4 * x[ind] + x[ind - 1] + x[ind + 1] + x[ind - M] + x[ind + M]

                print(ind)

                if (abs(R) > abs(Rmax)):

                    Rmax = R

                if (abs(R) > 10**-5):

                    print(R, ' at (', i, j, ')')

                    return False

    return True

# run the SOR process
def SOR(x, x_next, N, M, h, w, omega):

    for i in range(0, N):

        for j in range(0, M):

            ind = index_to_node(i, j, M)

            if ((j == M - 1) and (i > h-1 and i < N - 1)): # Neumann condition for right edge

                x_next[ind] = (4/3) * x_next[ind - 1] - (1/3) * x_next[ind - 2]

                x_next[ind] = (1 - omega) * x[ind] + omega * x_next[ind] # over relax with w

                #print ('At (' + str(i) + ',' + str(j) + ') we have ' + str(round(x_next[ind], 3)) + ' which is Neumann')

            elif ((i == 0) and (j > 0 and j < M - w)): # Neumann for top edge

                x_next[ind] = (4/3) * x[ind + M] - (1/3) * x[ind + 2*M]

                x_next[ind] = (1 - omega) * x[ind] + omega * x_next[ind] # ower relax with w
                #print ('At (' + str(i) + ',' + str(j) + ') we have ' + str(round(x_next[ind], 3)) + ' which is Neumann')

            elif (j == 0): # Dirchlet edge

                x_next[ind] = x[ind]

                #print ('At (' + str(i) + ',' + str(j) + ') we have ' +  str(round(x_next[ind], 3)) + ' which is Dirichlet')

            elif (i == N - 1): # Dirichlet edge

                x_next[ind] = x[ind]

                #print ('At (' + str(i) + ',' + str(j) + ') we have ' +  str(round(x_next[ind], 3)) + ' which is Dirichlet')


            elif ((i <= h - 1) and (j >= M - w)): # Dirichlet edge

                x_next[ind] = x[ind]

                #print ('At (' + str(i) + ',' + str(j) + ') we have ' +  str(round(x_next[ind], 3)) + ' which is Dirichlet')


            else: # free nodes

                x_next[ind] = 0.25 * (x_next[ind - 1] + x_next[ind - M] + x[ind + 1] + x[ind + M])

                x_next[ind] = (1 - omega) * x[ind] + omega * x_next[ind] # over relax with w

                #print ('At (' + str(i) + ',' + str(j) + ') we have ' +  str(round(x_next[ind], 3)) + ' which is free')


# run the Jacobi process
def Jacobi(x, x_next, N, M, h, w):

    for i in range(0, N):

        for j in range(0, M):

            ind = index_to_node(i, j, M)

            if ((j == M - 1) and (i > h-1 and i < N - 1)): # Neumann condition for right edge

                x_next[ind] = (4/3) * x[ind - 1] - (1/3) * x[ind - 2]

                print ('At (' + str(i) + ',' + str(j) + ') we have ' + str(round(x_next[ind], 3)) + ' which is Neumann')

            elif ((i == 0) and (j > 0 and j < M - w)): # Neumann for top edge

                x_next[ind] = (4/3) * x[ind + M] - (1/3) * x[ind + 2*M]

                print ('At (' + str(i) + ',' + str(j) + ') we have ' + str(round(x_next[ind], 3)) + ' which is Neumann')

            elif (j == 0): # Dirchlet edge

                x_next[ind] = x[ind]

                print ('At (' + str(i) + ',' + str(j) + ') we have ' +  str(round(x_next[ind], 3)) + ' which is Dirichlet')

            elif (i == N - 1): # Dirichlet edge

                x_next[ind] = x[ind]

                print ('At (' + str(i) + ',' + str(j) + ') we have ' +  str(round(x_next[ind], 3)) + ' which is Dirichlet')


            elif ((i <= h - 1) and (j >= M - w)): # Dirichlet edge

                x_next[ind] = x[ind]

                print ('At (' + str(i) + ',' + str(j) + ') we have ' +  str(round(x_next[ind], 3)) + ' which is Dirichlet')


            else: # free nodes

                x_next[ind] = 0.25 * (x[ind - 1] + x[ind - M] + x[ind + 1] + x[ind + M])

                print ('At (' + str(i) + ',' + str(j) + ') we have ' +  str(round(x_next[ind], 3)) + ' which is free')



# perform successive over relaxation
def run(x, x_next, N, M, h, w, omega):

    k = 0

    while True:

        #SOR(x, x_next, N, M, h, w, omega)

        Jacobi(x, x_next, N, M, h, w)

        print('The potential at node i = 1, j = 1 is', x_next[index_to_node(1,1, M)])

        if (residuals_satisfied(x_next, N, M, h, w)):

            break

        x = x_next[:]

        print('iteration: ', k)

        k += 1

    print('The number of iterations is: ', k)

# write to file
def potential_at(a, b, h_x, h_y, N, M, V):

    i = int(a/h_x)

    j = int(b/h_y)

    ind = index_to_node(N - 1 - j, i, M)

    return V[ind]

def main():

    os.chdir('C:/Users/aagni/Google Drive/McGill/U3/ECSE 543/Assignment 1/1.3')

    s_out = 0.2 # sive length of the region in the problem

    s = s_out / 2 # side length of region that we will solve

    h_c_og = 0.04 # height of conductor

    h_c = h_c_og / 2 # height of conductor in the region where we will solve

    w_c_og = 0.08 # width of conductor

    w_c = w_c_og / 2 # width of conductor in the region where we will solve

    spacing = 0.01

    h_x = spacing # space between successive mesh nodes in the x direction

    h_y = spacing # space between successive mesh nodes in the y direction

    N = int(s/h_x)

    M = int(s/h_y)

    h = int(h_c/h_y)

    w = int(w_c/h_x)

    print('N = ' + str(N))
    print('M = ' + str(M))
    print('h = ' + str(h))
    print('w = ' + str(w))

    V = 15

    b = boundary_condition(N, M, h, w, V)

    x = b

    x_next = [0] * (N*M)

    omega = 1.35

    run(x, x_next, N, M, h, w, omega)

    print('The potential at (0.06, 0.04) is V = ' + str(potential_at(0.06, 0.04, h_x, h_y, N, M, x_next)) + ' V')

if __name__ == '__main__':

    main()

