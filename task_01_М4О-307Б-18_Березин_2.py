import matplotlib.pyplot as plt
import numpy as np
import math
import os

#Папка для txt файла  результатами
try:
 os.mkdir('results')
except OSError:
 pass

complete_file = os.path.join('results', 'task_01_307Б_Березин.txt')
print("Текущая деректория:", os.getcwd())

#Заданная функция
def f(x):
    return -20*np.exp(-0.2*np.sqrt(0.5*x**2))-np.exp(0.5*(np.cos(2*math.pi*x)+1))+np.e+20

for x in list(range(-5, 6, 1)):
    y=-20*np.exp(-0.2*np.sqrt(0.5*x**2))-np.exp(0.5*(np.cos(2*math.pi*x)+1))+np.e+20
    print(x,y)

x= np.linspace(-5,5,300)
y=f(x)

plt.plot(x, y)
plt.xlabel("x")
plt.ylabel("y")
plt.grid()
plt.show()

