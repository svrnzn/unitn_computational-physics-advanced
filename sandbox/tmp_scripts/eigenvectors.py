import sys
import numpy
import matrix

ITERATIONS = int(sys.argv[1])
eigenvectors = numpy.zeros( (matrix.N, matrix.N) )
for i in range(ITERATIONS):
    eigenvalues, eigenvectors = matrix.radial_r(eigenvectors[:, 0])
    for i in range(1, matrix.N):
        eigenvectors[i] = eigenvectors[i] / (i * matrix.Dx)
    eigenvectors[0] = eigenvectors[1]#DIRT
    for i in range(matrix.N):
        eigenvectors[:, i] = eigenvectors[:,i] / (numpy.sum(eigenvectors[:,i] * matrix.Dx))

print('eigenvalues', end = '\t')
for i in range(matrix.N):
    print('"{}"'.format(eigenvalues[i]), end = '\t')
print()
for i in range(matrix.N):
    print(i * matrix.Dx, end = '\t')
    for j in range(matrix.N):#i want them to be printed in columns
        print(eigenvectors[i][j], end = '\t')
    print()
