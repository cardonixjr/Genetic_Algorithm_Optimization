import sys
import numpy as np

try:
	import matplotlib.pyplot as plt 		# Library to plot
except:
	print("matplotlib.pyplot library not found, unable to plot your function.")
	print("Try to download it:")
	print("https://matplotlib.org/3.1.1/users/installing.html")

	sys.exit()

try:
	from mpl_toolkits.mplot3d import Axes3D
except:
	print("mpl_toolkits library not found, unable to plot your function.")
	print("Try to download it:")
	print("https://stackoverflow.com/questions/37661119/python-mpl-toolkits-installation-issue")

def my_plot(f, min_x, max_x, min_y, max_y, markers = []):
	subdivisions = 50.0

	x_values = np.arange(min_x,max_x, abs((max_x - min_x))/subdivisions)
	y_values = np.arange(min_y,max_y, abs((max_y - min_y))/subdivisions)
	x_values, y_values= np.meshgrid(x_values, y_values)
	z_values = f(x_values,y_values)



	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	ax.plot_surface(x_values,y_values,z_values)

	for marker in markers:
		ax.scatter(marker[0], marker[1], marker[2], s=30, c="red")

	plt.show()

if __name__ == "__main__":
	f = lambda x,y: x**2 + y**2

	my_plot(f,-15,15,-3,3, [[1,1,f(1,1)], [-5,-1, f(-5,-1	)], [10,0,f(10,0)]])