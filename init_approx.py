import numpy as np
import matplotlib.pyplot as plt


def f(k, b, alpha):
	return (k/alpha)*(-b*k**(alpha) + (alpha+1))

def F(k, alpha, b = 1, d=8):
	for i in range(d):
		k = f(k, b, alpha)
	return k

def find_bounds(d, alpha):
	delta = 1e-5
	#Upper Bound
	# l = 1; r = 2*3**0.5 - 1; err = 7e-3
	l = 1; r = 2*3**0.5 - 1; err = 1e-4
	# l = 1; r = 3**0.5; err = 7e-3
	while r-l >= delta:
		mid = (l+r)/2
		val = abs(F(mid, alpha, d=d) - 1)
		# print(l, r, mid, F(mid, alpha, d=d), val)
		if val <= err:
			l = mid 
		else:
			r = mid - delta
	# print(l)
	k1 = l

	# print('-'*20)
	#Lower Bound
	# l = 0; r = 1
	l = 2*delta-1; r = 1
	while r-l >= delta:
		mid = (l+r)/2
		val = abs(F(mid, alpha, d=d) - 1)
		# print(l, r, mid, F(mid, alpha, d=d), val)
		if val <= err:
			r = mid
		else:
			l = mid + delta
	# print(r)
	k2 = r

	return k1, k2

if __name__ == "__main__":
	print(find_bounds(8, 2))