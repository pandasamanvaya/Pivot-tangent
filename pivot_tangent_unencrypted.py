import numpy as np
from scipy import optimize
from matplotlib import pyplot as plt
from init_approx import find_bounds

np.random.seed(0)

a = 1e-4; b = 1e3; d=9
# a = 1e-3; b=750; d=8

(k1, k2) = find_bounds(d, 2)
def f(x):
  return (35*(x) - 35*(x)**3 + 21*(x)**5 - 5*(x)**7)/(2**4)

def g(x):
	return (4589*(x) - 16577*(x)**3 + 25614*(x)**5 - 12860*(x)**7)/(2**10)

def sign(x, d_g, d_f):
	for i in range(d_g):
		x = g(x)
		# print(max(x))
	for i in range(d_f):
		x = f(x)
	return x

def comp(x, p):
	return (sign((x - p)/b, d_g, d_f)+1)/2

def getIndex(ele):
	l = 0; r = N-1
	# 	(ele)
	while l < r:
		mid = (l+r)//2
		# print(mid, x[mid], ele, l, r)
		if abs(x[mid]-ele) < 1e-5:
			return mid
		elif x[mid] < ele:
			l = mid+1
		else:
			r = mid-1
	return l

def iterInd(i):
	if i == 1:
		return '1st'
	elif i == 2:
		return '2nd'
	else:
		return str(i)+'th'

def eq1(k1, k2, x1):
	coeff = [-4*k2**2, 9*x1*k1**2, -6*(x1**2)*k1**2, (k1**2)*(x1**3)]
	# print(coeff)
	return coeff

def eq2(k1, k2, x1):
	coeff = [k1**2, -6*x1*k1**2, 9*(k1**2)*x1**2, -4*(k2**2)*x1**3]
	# print(coeff)
	return coeff

c1 = eq1(k1, k2, b)
def p(x):
	return c1[0]*x**3 + c1[1]*x**2 + c1[2]*x + c1[3]

def p_prime(x):
	return 3*c1[0]*x**2 + 2*c1[1]*x + c1[2]

def q(x):
	return c2[0]*x**3 + c2[1]*x**2 + c2[2]*x + c2[3]

def q_prime(x):
	return 3*c2[0]*x**2 + 2*c2[1]*x + c2[2]

def line1(x, x1):
	# print(-0.5*k1*x1**-1.5, 1.5*k1*x1**-0.5)
	return (-0.5*k1*x1**-1.5)*x + 1.5*k1*x1**-0.5

def line(x, x1):
	m = (k2*b**-0.5 - k2*x1**-0.5)/(b-x1)
	return m*(x-b) + k2*b**-0.5

sol = optimize.root_scalar(p, fprime=p_prime, method='newton', x0=400)
# print(sol)
ans1 = sol.root
# print("First tangent point =", ans1)
# print(k2*b**-0.5, line1(b, ans))

c2 = eq2(k1, k2, ans1)
point = optimize.root_scalar(q, fprime=q_prime, method='newton', x0=1)
pivot = point.root
# print("Pivot point =", pivot)
c1 = eq1(k1, k2, pivot)
x0 = a
while p_prime(x0) < 0:
	x0 += 0.001
	sol = optimize.root_scalar(p, fprime=p_prime, method='newton', x0=x0)
ans2 = sol.root
# print("Second tangent point =", ans2)
# print(line1(pivot, ans2), k2*pivot**-0.5)

N = 2**13; d_g = 7; d_f = 2
points1 = np.linspace(1, b, N//2)
points2 = np.linspace(a, 1, N//2)
print("k1 =", k1, "k2 =", k2, "x1 =", ans1, "x2 =", ans2, "pivot =", pivot)

x = np.append(points2, points1)
# x = np.linspace(a, b, N)
y = x**-0.5

fact = comp(x, pivot)
y_pred = (1-fact)*line1(x, ans2) + fact*line1(x, ans1)

print("Initial error =", np.mean((y-y_pred)**2))
for i in range(d):
	y_pred = (y_pred/2)*(-x*y_pred**2 + 3)
	print("Error after %s iteration = %f and max error = %f"%(iterInd(i+1), np.mean(abs(y-y_pred)), np.max(abs(y-y_pred))))
print(y_pred[0])

# ind = getIndex(pivot)

# plt.plot(x, y, label='Original')
# plt.plot(x, k1*x**-0.5, label='Upper bound')
# plt.plot(x, k2*x**-0.5, label='Lower bound')
# plt.plot(x, y_pred, label='Predicted')
# plt.plot(x[:ind], line1(x[:ind], ans2), label='L1')
# plt.plot(x[ind:], line1(x[ind:], ans1), label='L2')
# plt.xlim(-0.5, 8)
# plt.ylim(-0.5, 8)
# plt.legend()
# plt.savefig('inv_sqrt.png', dpi=500)
# plt.show()
