import sys
import numpy
import matrix

ITERATIONS = int(sys.argv[1])
eigenvectors = numpy.zeros( (matrix.N, matrix.N) )
for i in range(ITERATIONS - 1):
    eigenvalues, eigenvectors = matrix.radial_r(eigenvectors[:, 0])
    for i in range(1, matrix.N):
        eigenvectors[i] = eigenvectors[i] / (i * matrix.Dx)
    eigenvectors[0] = eigenvectors[1]#DIRT
    for i in range(matrix.N):
        eigenvectors[:, i] = eigenvectors[:,i] / (numpy.sum(eigenvectors[:,i] * matrix.Dx))

for i in range(matrix.N):
    print(i * matrix.Dx, matrix.effective_potential((i * matrix.Dx), eigenvectors[:, 0]))
