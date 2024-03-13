import numpy as np

# phân rã ma trận A=LL.T
def decomposition(matrix):
	a = np.array(matrix, dtype = complex)
	n,_ = np.shape(a)
	L = np.zeros((n,n), complex)
	for j in range(n):
		for i in range(j,n):
			if i == j:
				sumk = 0
				for k in range(j):
					sumk += L[j,k]**2
				L[i,j] = np.sqrt(a[i,j]-sumk)
			else:
				sumk = 0
				for k in range(j):
					sumk += L[i,k]*L[j,k]
				L[i,j] = (a[i,j]-sumk)/L[j,j]
	print('ma tran tam giac duoi:\n', L, '\n')
	return L

# Det = 0 ??
def DetZeros(matrix):
	n,_ = np.shape(matrix)
	prod = prod2 = 1
	for i in range(n):
		prod *= matrix[i,i]
		prod2 *= (matrix[i,i]**2)
	if prod == 0:
		print('Det(A) = 0')
		return True
	else:
		print('Det(A)=', prod2.real)
		return False

# ma trận đối xứng ??
def IsSymmetry(matrix):
	n,_ = np.shape(matrix)
	for i in range(n):
		for j in range(i+1,n):
			if matrix[i,j] != matrix[j,i]:
				print('A ko doi xung')
				return False
	print('A doi xung')
	return True

# Giải hệ Ax = B (A là ma trận tam giác trên/dưới)
def Solve(matrixA, matrixB, L = True):
	y = np.zeros_like(matrixB, complex)
	x = np.zeros_like(matrixB, complex)
	n,_ = np.shape(matrixB)
	for i in range(n):
		sumk_y = sumk_x = 0
		j = n-1-i
		for k in range(i):
			sumk_y += matrixA[i,k] * y[k,0]
			sumk_x += matrixA[j,n-1-k] * x[n-1-k,0]
		y[i,0] = (matrixB[i,0]-sumk_y) / matrixA[i,i]
		x[j,0] = (matrixB[j,0]-sumk_x) / matrixA[j,j]
	if L == True:
		return y
	else:
		return x.real

def Cholesky(matrixA,matrixB):
	doixung = IsSymmetry(matrixA)
	if doixung == True:
		print("Ma tran thoa man")
	if doixung == False:
		matrixB = matrixA.T @ matrixB
		matrixA = matrixA.T @ matrixA
		print('ma tran sau khi bien doi ve doi xung: A=', matrixA, 'B=', matrixB)
	L = decomposition(matrixA)
	det0 = DetZeros(L)
	if det0 == True:
		print('Ko giai dc, det = 0')
	print('Nghiem cua he pt la:')
	y = Solve(L,matrixB,True)
	x = Solve(L.T,y,False)
	return x

A = np.array([[1,3,-2,0,-2],
			[3,4,-5,1,-3],
			[-2,-5,3,-2,2],
			[0,1,-2,5,3],
			[-2,-3,2,3,4]])
B = np.array([[0.5],
			[5.4],
			[5],
			[7.5],
			[3.3]])

A1 = np.array([[-1,2,1],
			[2,1,-3],
			[1,-3,-2]])
B1 = np.array([[1],
			[2],
			[3]])
#print(np.linalg.solve(A,B))
print(Cholesky(A1,B1))