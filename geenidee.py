import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression


#independent variables
thl = np.array([0.2, 0.22, 0.24, 0.26, 0.28, 0.3])
qt = np.array([1e-4, 11e-5, 12e-5, 13e-5, 14e-5, 15e-5])

x1 = thl.reshape((-1,1))
x2 = qt.reshape((-1,1))
#dependent variables
lwp_thl = np.array([0.062, 1.143, 2.094, 1.047, 0.245, 8.937])
rwp_thl = np.array([5.1e-7, 3.58e-5, 5e-5, 7.75e-6, 1.66e-6, 2.53e-4])
ql_thl = np.array([0.29, 1.03, 1.77, 1.56, 0.67, 3.72])
cf_thl = np.array([0.01, 0.16, 0.36, 0.1, 0.07, 1.07])

lwp_qt = np.array([0.062, 2.095, 2.849, 5.705, 6.158, 10.225])
rwp_qt = np.array([5.1e-7, 5.78e-5, 1.1e-4, 3.84e-4, 6.76e-4, 2.73e-3])
ql_qt = np.array([0.29, 1.98, 2.25, 2.99, 4.42, 6.79])
cf_qt = np.array([0.01, 0.08, 0.14, 0.56, 0.33, 0.63])

x = x2
y = rwp_qt

model = LinearRegression()
model.fit(x,y)

r_sq = model.score(x,y)
print("coefficient of determination:", r_sq)	#R^2
print("intercept:", model.intercept_)   	#b0
print("slope:", model.coef_)			#b1


y_pred = model.intercept_ + model.coef_ * x 		#model.predict(x)
print('predicted response:', y_pred, sep='\n')

#Plots
plt.scatter(x, y, color="black")
plt.plot(x, y_pred, linewidth=2, label='R$^2$ = {:.2f}'.format(r_sq))
plt.xlabel('qt [kg kg-1]')
plt.ylabel('RWP [g m-2]')
plt.legend()

plt.show()

