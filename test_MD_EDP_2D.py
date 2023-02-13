
import scipy.integrate as itg
import numpy as np
import matplotlib.pyplot as plt


#paramètres utilisés
K = 50
N = 50
L = 50
pas = 0.1

t_in = 0.3
t_out = 2.4
t_open = 130
t_close = 150
v_gate = 0.13
v_stim = 1
t_stim = 0.1

def MS_2D(pas = 0.1, K = 50, N = 50, L = 50, Tmax = 500, t_in = 0.03, t_out = 0.6, t_open = 12, t_close = 3, v_gate = 0.1, t_stim = 0.1):

  M = np.array([[0.0 for i in range(N**2)] for j in range(N**2)])
  for i in range(1,N-1):
    for j in range(1, N-1):
      M[N*i +j][N*i+j] = -4
      M[N*(i-1)+ j][N*i + j] = 1
      M[N*(i+1)+j][N*i+j] = 1
      M[N*i+ j-1][N*i + j] = 1
      M[N*i + j+1][N*i +j] = 1
    M[N*i][N*i] = -3
    M[N*(i+1)][N*i] = 1
    M[N*i + N-1][N*i +N-1] = -3
    M[N*(i-1) + N-1][N*i +N-1] = 1
  for j in range(N):
    M[j][j] = -3
    M[j+1][j] = 1
    M[j+N*(N-1)][j+N*(N-1)] = -3
    M[j +N*(N-1) -1][j + N*(N-1)] = 1
  M[0][0] = -2
  M[N**2 - 1][N**2 -1] = -2
  M[N*(N-1)][N*(N-1)] = -2
  M[N-1][N-1] = -2
  M = (N/L)**2 * M



  def dh_dt(vh, t):
    v = vh[:N**2]
    h = vh[N**2:]
    dh = [0. for i in range(N**2)]
    for k in range(N**2):
      if v[k] < v_gate:
        dh[k] = (1-h[k]) /t_open
      else : 
        dh[k] =  -h[k]/t_close
    return np.array(dh)

  def J_stim(t):
    I = np.zeros(N**2)
    if t<t_stim :
      for k in range(int(0.1*N)):
        for j in range((int(0.1*N))):
          I[N*N//2+N//2] = 1
      return v_stim * I
    else : 
      return I
  
  np.vectorize(J_stim)

  def dv_dt(vh, t):
    v = vh[:N**2]
    h = vh[N**2:]
    dv = K*np.dot(v, M) + J_stim(t)+ h/t_in * (np.ones(N**2)-v)*(v**2) - v/t_out 

    return np.array(dv)

  v0 = np.array([0 for i in range(N**2)])
  h0 = np.array([1 for i in range(N**2)])
  vh0 = np.concatenate((v0, h0), axis = None)  
  
  def F(t, Y):
   a = dv_dt(Y, t)
   b = dh_dt(Y, t)
   return  np.concatenate((a, b), axis = None)
  
  t = np.arange(0., Tmax, pas)
  X = np.linspace(0, L, N)
  l = itg.odeint(F, vh0, t, tfirst = True)

  return l

L = 30
N = 30
Tmax = 200
K = 0.1
pas = 0.05
t_stim = 0.001*Tmax


u = MS_2D(L = L, N = N, Tmax = Tmax, K = K, t_stim = t_stim, pas = pas)
t = np.arange(0, Tmax, pas)
x = np.arange(0, L, L/N)
print(u.shape)

#plt.plot(x, u[300, :N])
#plt.xlabel("x")
#plt.ylabel("u(x)")
#plt.show()

#plt.plot(t, u[:, N//2])
#plt.xlabel("t")
#plt.ylabel("u(t)")
#plt.show()


from IPython.display import HTML
import matplotlib.animation  as animation
from mpl_toolkits.mplot3d import Axes3D 

plt.rcParams["figure.figsize"] = [7.00, 3.50]
plt.rcParams["figure.autolayout"] = True

fig3 = plt.figure()
ax3 = plt.axes(projection='3d') 
ax3.set_zlim(0,1.2)

x = np.arange(0, L, L/N)
y = np.arange(0, L, L/N)
x, y = np.meshgrid(x, y)
zarray = []
for t in range(0, int(Tmax/pas), 10):
  z = np.array([[0. for k in range(N)] for j in range(N)])
  for i in range(N):
      for j in range(N):
          z[i][j] = u[t,  N*i + j]
  zarray.append(z)

def surface_plot(t):
  plot[0].remove()
  plot[0] = ax3.plot_surface(x, y, zarray[t], cmap = 'viridis')

plot = [ax3.plot_surface(x, y, zarray[0], cmap = 'viridis')]
anim3 = animation.FuncAnimation(fig3, surface_plot,frames = 400, interval=100, repeat = True)
plt.rc('animation', html='jshtml')
plt.show()
#writergif = animation.PillowWriter(fps=1)
#fn = 'plot_surface_animation_funcanimation'
#anim3.save(fn+'.gif',writer=writergif)


