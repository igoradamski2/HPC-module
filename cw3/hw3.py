"""M3C 2017 Homework 3 solution
Igor Adamski, CID: 01069152
"""
import numpy as np
import matplotlib.pyplot as plt
import timeit, time
import scipy.spatial as spatial
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.axes3d as p3
import math
import matplotlib.animation as animation
from m1 import hw3 #assumes hw3_template.f90 has been compiled with f2py to generate m1.so

def rwnet(H,n,dim=2,debug=False):
    """
    Simulate random network model
    Input variables
    H: Height at which nodes are initially introduced
    n: Total number of nodes
    dim: dimension of node coordinates
    Output variables
    X: Final node coordinates
    output: a tuple containing any other information you would
    like the function to return, may be left empty
    """
    output=()
    dstar = 1
    X = np.zeros((n,dim))

    for i in range(1,n):
        XN = np.squeeze(X[:i,:]) # current nodes
        xx = np.zeros(dim) # location of new node
        xx[-1] = H

        deltaX = np.max(np.abs(XN-np.outer(np.ones(i),xx)),axis=1)
        dist = np.min(deltaX) #"distance" to nearest node

        #Compute trajectory of node, stopping when appropriate
        while xx[-1]>0 and dist>dstar:
            R = np.random.randint(0,2,dim)
            R[0:-1] = R[0:-1]*2-1
            R[-1] = -R[-1]
            xx = xx + R
            deltaX = np.max(np.abs(XN-np.outer(np.ones(i),xx)),axis=1)
            dist = np.min(deltaX)

        X[i,:] = xx

        if debug:
            if np.mod(i,100)==0:
                print("i,xx,yy=",i,xx)

    return X,output

def performance1():
    """Analyze performance of python and serial fortran codes
    Modify input/output as needed

    On the first figure we can see the average time (average of 5 runs)
    it takes python to simulate a network with H=140 and varying
    numbers of points and dimensions. We can see clearly a pattern that
    it takes less time for python to simulate the network of dimension 2
    than it takes to simulate a network of dimension 3, which is kind of
    a obvious observation. The second observation is that the greater the
    number of points the greater the time it takes to run. We can also see
    that there is a decline in time moving from 1000 nodes to 14000
    nodes. This is an effect of the inate randomness of the generation algorithm
    because the times should be the same, as clearly above 1000 nodes
    the network height is big enough that all subsequent nodes spawned
    at the position (0,H) dont do the walk at all because they are
    in the neighbourhood of other nodes immediately after spawning.

    On the second figure we can see the same relationship but this
    time the networks have been simulated by fortran. We can see that
    the trends depending on number of points and dimension are the
    same: the greater the number of points and dimension the longer it
    takes to simulate. It is not surprising because the more nodes we want to
    add to the network the more calculations we need to execute. Also in 3d its
    a lot slower because the random walk has one more dimension so we need to
    make one more calculation at each step. Between 1000 and 1400 number of points
    we see no speedup in the 2d case because of the reason described above -
    the network is high enough that new nodes dont need to make the random walk.

    We can also see that the effect of dimension is greater the more points we are
    simulating. For small number of points, the times to generate the 2d and 3d
    network are comparable, because the effect of calculating in one more dimension
    is not visible at such a small scale.

    In the second figure we for n=10 we cant see the times for the fortran code, that
    is because they are too small to be displayed on the graph.

    We can also see that the simulation done in fortran is much faster than
    simulation done in python. Much faster in this context means that fortran is
    faster by an order of magnitude approximately 20 (in some cases for 500 nodes
    fortran is 30 times faster). This is not surprising, as we know that
    fortran is a compiled language and therefore has its for loops maximally optimiezed.
    On the other hand python is a interpreted language and the loops in python are not
    optimized as much as in fortran.
    """
    #initialize number of points and dimensions
    num_points = [10, 100, 500, 1000, 1400]
    dim = [2, 3]

    #keep some variables constant so we can analyze the time
    H = 140

    #initialize the vector in which we'll keep the python and fortran times
    times_py = np.zeros([len(num_points), len(dim)])
    times_f90 = np.zeros([len(num_points), len(dim)])

    #loop
    i = 0
    j = 0
    for d in dim:
        i = 0
        for n in num_points:
            print('we just checked n=',num_points[i])
            def rw_py():
                return rwnet(H, n, d)
            def rw_f90():
                return hw3.rwnet(H, d, n, 0)
            times_py[i, j] = timeit.timeit(rw_py, number = 5)/float(5)
            times_f90[i, j] = timeit.timeit(rw_f90, number = 5)/float(5)
            i += 1
        j += 1

    x = list(range(1, len(num_points) + 1))
    #fig1 shows time varying for number of points and dimensionality
    fig1 = plt.figure(figsize=(10,7))
    ax = fig1.add_subplot(111)
    ax.bar(np.array(x)-0.1, times_py[:,0], width=0.2, color='b', align='center', label='dim 2')
    ax.bar(np.array(x)+0.1, times_py[:,1], width=0.2, color='r', align='center', label='dim 3')
    ax.set_xticks(x)
    ax.set_xticklabels(num_points)
    ax.set_xlabel('Number of points in the network')
    ax.set_ylabel('Mean time to generate network(s)')
    ax.set_title('Average times (mean of 5 runs) to generate networks in Python for different dimensions \n and number of points, keeping H=140. \n Igor Adamski, performance1()')
    for i, v in enumerate(times_py[:,0]):
        if v*100 < 1:
            continue
        ax.text(i+0.79, v+0.01, "{0:.2f}".format(v), color='black')
    for i, v in enumerate(times_py[:,1]):
        if v*100 < 1:
            continue
        ax.text(i+0.99, v+0.01, "{0:.2f}".format(v), color='black')
    ax.legend()

    plt.show()

    #fig2 shows the same but for gfortran
    fig2 = plt.figure(figsize=(10,7))
    ax1 = fig2.add_subplot(111)
    ax1.bar(np.array(x)-0.1, times_f90[:,0], width=0.2, color='b', align='center', label='dim 2')
    ax1.bar(np.array(x)+0.1, times_f90[:,1], width=0.2, color='r', align='center', label='dim 3')
    ax1.set_xticks(x)
    ax1.set_xticklabels(num_points)
    ax1.set_xlabel('Number of points in the network')
    ax1.set_ylabel('Mean time to generate network (s)')
    ax1.set_title('Average times (mean of 5 runs) to generate networks in Fortran for different dimensions \n and number of points, keeping H=140. \n Igor Adamski, performance1()')
    for i, v in enumerate(times_f90[:,0]):
        if v*100 < 1:
            continue
        ax1.text(i+0.79, v+0.01, "{0:.2f}".format(v), color='black')
    for i, v in enumerate(times_f90[:,1]):
        if v*100 < 1:
            continue
        ax1.text(i+0.99, v+0.01, "{0:.2f}".format(v), color='black')
    ax1.legend()

    plt.show()

    return None



def performance2():
    """Analyze performance of parallel m-simulation fortran code
    Modify input/output as needed

    On the created figure we can see a very self explanatory trend,
    that using twice as many cores speeds up the generation of m networks
    twice. We can see that the trend is visible straight from m=2.
    However for m=2 the times are very little either way (take less than
    1 second) and we can attribute a speedup of perfectly twice to randomness.
    We can see that for m=10 the speedup is almost twice, and for subsequent
    m's the speedup is either exactly twice or a bit more (also due to randomness).

    The trend is therefore that using twice as many cores speeds up the
    calculation twice, and the trend is more visible for greater m's.
    """
    #initialize variables
    h = 600
    dim = 2
    n = 100
    m = [2, 10, 50, 100, 200]

    times_1 = np.zeros([len(m)])
    times_2 = np.zeros([len(m)])

    i = 0
    for simulations in m:
        print('im calculating for m=', m[i])
        def one_thread():
            return hw3.simulate_m(h,dim,n,simulations,1)

        def two_threads():
            return hw3.simulate_m(h,dim,n,simulations,2)

        times_1[i] = timeit.timeit(one_thread, number=3)/float(3)
        times_2[i] = timeit.timeit(two_threads, number=3)/float(3)

        i += 1

    fig1 = plt.figure(figsize=(7,6))
    ax = fig1.add_subplot(111)
    x = list(range(1,len(m)+1))
    ax.bar(np.array(x)-0.1, times_1, width=0.2, color='b', align='center', label='1 core')
    ax.bar(np.array(x)+0.1, times_2, width=0.2, color='r', align='center', label='2 cores')
    ax.set_xticks(x)
    ax.set_xticklabels(m)
    ax.set_xlabel('Number of simulations (m)')
    ax.set_ylabel('Time (mean of 3 runs in seconds)')
    ax.set_title('Comparison of simulating m networks using 1 and 2 cores. \n Igor Adamski, performance2().')
    for i, v in enumerate(times_1):
        ax.text(i+0.79, v+0.01, "{0:.1f}".format(v), color='black')
    for i, v in enumerate(times_2):
        ax.text(i+0.99, v+0.01, "{0:.1f}".format(v), color='black')
    ax.legend()
    plt.show()


    return None

def analyze1():
    """Analyze dependence of network height on number of nodes
    Modify input/output as needed

    In the figure I am plotting the mean height fo the network
    from running 10 simulations for different numbers of points
    in the network. We can see that the mean height for the
    dimension 2 graphs is a lot bigger than for dimension 3 graphs
    which is not surprising as dropping a point in 3d we have an even
    smaller probability to come closer to a previous node because
    we are moving in 3 directions instead of 2. We can see on the
    graphs that as the number of points increases the mean network
    height increases (which is not surprising), until it reaches a certain point
    and the height stays constant. This is because as the number of points get
    bigger, our network will reach the height H before all the points are added to
    the network and thus we will have a situation when the node will be seeded
    at (0,H) and will immediately stay there because it will connect to other nodes
    even before it starts its random walk. We can see that for small n the 2d network
    height grows rather linearly and fast, and then the height grows stalls when its
    nearing H. Once we reach a sufficent n we should see the mean heights to be at H.
    We can see that the network height for the 3d case grows almost linearly with the
    number of points added to the network. This also happens for 2d, with the exception
    that the height stalls at some particular number of nodes, beyond which it cannot
    increase its height.

    The trend we can see on the graph, summarizing, is that as n (number of points) increases
    so do the mean heights of the networks, and they increase rather linearly until they reach a
    point when they have a constant (or roughly constant) mean height close to H.
    """
    #we will analyze the mean height from 50 simulations

    #initialize variables
    n = [250, 500, 1000, 2000, 3000, 4000, 5000]
    mean_h1 = np.zeros([len(n)])
    mean_h2 = np.zeros([len(n)])
    i = 0
    for numbers in n:
        print('im calculating n=', numbers)
        mean_h1[i] = hw3.simulate_m(200,2,numbers,10,2)[3]
        mean_h2[i] = hw3.simulate_m(100,3,numbers,10,2)[3]
        i += 1

    fig1 = plt.figure(figsize=(10,8))
    ax = fig1.add_subplot(111)

    x = list(range(1,len(n)+1))
    ax.bar(np.array(x)-0.1, mean_h1, width=0.2, color='b', align='center', label='Dimension 2, H=200')
    ax.bar(np.array(x)+0.1, mean_h2, width=0.2, color='r', align='center', label='Dimension 3, H=100')
    ax.set_xticks(x)
    ax.set_xticklabels(n)
    ax.set_xlabel('Number of network points')
    ax.set_ylabel('Mean network height (10 runs)')
    ax.set_title('Mean network heights (10 runs) for different numbers of points in the network. \n Igor Adamski, analyze1().')
    ax.legend()
    for i, v in enumerate(mean_h1):
        ax.text(i+0.79, v+0.01, "{0:.1f}".format(v), color='black')
    for i, v in enumerate(mean_h2):
        ax.text(i+0.99, v+0.01, "{0:.1f}".format(v), color='black')
    plt.show()

    return None


def analyze2():
    """ Estimate dimensions of simulated networks with d=2,3
        Output: estimated dimensions

        This function tries to fit a curve of the form given in the question
        to the data C(e). If we think about it, both deltas (for 2d and 3d) should
        converge to 1 as N goes to infinity. This is because as N goes to infinity
        we will have that most of the nodes will be concentrated just at one point
        of seed ((0,H) for 2d and (0,0,H) for 3d). Therefore since most of the
        network nodes will be in one place at the seed node, the fraction N(e)/N^2
        will tend to be linear in epsilon. For 2d the N=10000 that is being used in thus
        code is sufficient to result in delta close to 1 but for 3d 10000 is not enough.
        I hypothesise that if we did N=100000 then both deltas would be close to 1
        but such N would take too much time to run.

        Warning: This function takes quite a lot of time to run (round 15 minutes),
        but n is quite big and this is expected. If you want to make it more managable
        change n=5000 in the first line below.
    """
    #set n
    n=5000

    #run the 2d network generation and extract points
    positions = hw3.rwnet(4,2,n,0)
    positions = positions.T.reshape(-1,2)
    point_tree = spatial.cKDTree(positions)
    height = max(positions[:,1])

    #run the 3d network and extarct points
    positions_d3 = hw3.rwnet(3,3,n,0)
    positions_d3 = positions_d3.T.reshape(-1,3)
    point_tree_d3 = spatial.cKDTree(positions_d3)
    height_d3 = max(positions_d3[:,2])

    #initialize epsilons and empty arrays
    eps = np.linspace(2, height/2, 100)
    eps_d3 = np.linspace(2, height_d3/2, 100)
    C_eps_array = np.zeros([100])
    C_eps_array_d3 = np.zeros([100])

    #begin the loop over different epsilons
    index = 0
    for epsilon, epsilon_d3 in zip(eps, eps_d3):
        N_eps = 0
        N_eps_d3 = 0
        print("We are at iteration", index, "of", len(eps))
        for points, points_d3 in zip(positions, positions_d3):
            #get indeces of points connected to the current point
            indeces = point_tree.query_ball_point(points, epsilon)
            indeces_d3 = point_tree_d3.query_ball_point(points_d3, epsilon_d3)

            #now update the N_eps and N_eps_d3 (len(indeces) are never less than 1)
            N_eps += len(indeces) - 1
            N_eps_d3 += len(indeces_d3) - 1

        #no extract the true N_eps because we counted twice
        if N_eps % 2 == 1:
            N_eps = (N_eps + 1)/2
        else:
            N_eps = N_eps/2

        if N_eps_d3 % 2 == 1:
            N_eps_d3 = (N_eps_d3 + 1)/2
        else:
            N_eps_d3 = N_eps_d3/2

        #update array for 2d case
        C_eps = N_eps/(n**2)
        C_eps_array[index] = C_eps

        #update array for 3d case
        C_eps_d3 = N_eps_d3/(n**2)
        C_eps_array_d3[index] = C_eps_d3

        index += 1

    #define function to which we will fit curve
    def func(eps, delta, scale, shift):
        return scale*(eps**delta)+shift

    #fit curve for 2d case
    popt, pcov = curve_fit(func, eps, C_eps_array, maxfev = 10000)

    #fit curve for 3d case
    popt3, pcov3 = curve_fit(func, eps_d3, C_eps_array_d3, maxfev = 10000)

    #extract coefficients
    delta, scale, shift = popt
    delta2, scale2, shift2 = popt3

    fig = plt.figure(figsize=(10,7))
    ax = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax.plot(eps, C_eps_array, label='$C(\epsilon)$ from 2d data')
    ax.plot(eps, func(eps, delta, scale, shift), label='Fitted $C(\epsilon)$ ~ $\epsilon^\delta$')
    ax.set_xlabel('$\epsilon$')
    ax.set_ylabel('$C(\epsilon)$')
    ax.set_title('Plot of $C(\epsilon)$ vs. $\epsilon$ for the 2D network \n using N=7000. $\delta$ \u2248 {:.2f}. Igor Adamski, analyze2()'.format(delta))
    ax2.plot(eps_d3, C_eps_array_d3, label='$C(\epsilon)$ from 3d data')
    ax2.plot(eps_d3, func(eps_d3, delta2, scale2, shift2), label='Fitted $C(\epsilon)$ ~ $\epsilon^\delta$')
    ax2.set_xlabel('$\epsilon$')
    ax2.set_ylabel('$C(\epsilon)$')
    ax2.set_title('Plot of $C(\epsilon)$ vs. $\epsilon$ for the 3D network \n using N=7000. $\delta$ \u2248 {:.2f}. Igor Adamski, analyze2()'.format(delta2))
    ax.legend()
    ax2.legend()
    plt.show()

    return (delta, delta2)


def visualize():
    """Generate an animation illustrating the generation of a network
    """
    #first simulate the network and we will be adding
    #the nodes to the animation one by one
    positions = hw3.rwnet(200,3,500,0)
    positions = positions.T.reshape(-1,3)
    point_tree = spatial.cKDTree(positions)
    x = positions[:,0]
    y = positions[:,1]
    z = positions[:,2]

    #this bit of code is my failed attempt to add
    #edges between nodes, with more work it should work
    #but for now dont uncomment it
    '''
    def get_edges():
        index = 0
        dic = {}
        for points in positions:
            indeces = point_tree.query_ball_point(points, 1)
            indeces.remove(index)
            dic[index] = indeces
            index += 1
        return dic

    dic = get_edges()

    X_s = []
    Y_s = []
    Z_s = []
    Z_e = []
    X_e = []
    Y_e = []

    def add_edges(dic, j):
        #get the connections

        if len(dic[j]) != 0:
            for a in range(len(dic[j])):
                if dic[j][a] > j:
                    continue
                else:
                    X_s.append(x[j])
                    Y_s.append(y[j])
                    Z_s.append(z[j])
                    X_e.append(x[dic[j][a]])
                    Y_e.append(y[dic[j][a]])
                    Z_e.append(z[dic[j][a]])
                    print("we are connecting ", j, "node to ", dic[j][a], "node")
        return X_s,Y_s,Z_s,X_e,Y_e,Z_e


    for i in range(500):
        if len(dic[i]) != 0:
            for points in dic[i]:
                X_s.append(x[i])
                Y_s.append(y[i])
                Z_s.append(z[i])
                X_e.append(x[points])
                Y_e.append(y[points])
                Z_e.append(z[points])
    '''

    def update_figure(i,x,y,z):
        #gather x and y data
        xy_data = np.c_[x[:i], y[:i]].T
        #simulate the axis rotation as sine
        g=np.linspace(0,4*math.pi,500)
        g=45*np.sin(g)+45
        #rotate axis every iteration
        ax.view_init(30, azim = g[i])
        #add new point to data
        line.set_data(xy_data)
        line.set_3d_properties(z[:i])
        #make the newest point have blue color
        xy_curr = np.c_[x[i],y[i]].T
        line_temp.set_data(xy_curr)
        line_temp.set_3d_properties(z[i])
        return line, line_temp

    fig = plt.figure(figsize=(8,6))
    ax = p3.Axes3D(fig)
    ax.set_xlim3d([min(x)-1, max(x)+1])
    ax.set_ylim3d([min(y)-1, max(y)+1])
    ax.set_zlim3d([0, max(z)])
    ax.set_xlabel('x-axis')
    ax.set_ylabel('y-axis')
    ax.set_zlabel('z-axis')
    ax.set_title('Animation of hw3.rwnet(200,3,500,0) network generation. Igor Adamski, visualize().')
    line = ax.plot(x[:1], y[:1], z[:1], 'o', color='red')[0]
    line_temp = ax.plot([],[],[], 'o', color='blue')[0]

    line_ani = animation.FuncAnimation(fig, update_figure, fargs=(x,y,z), frames=500,interval=10, blit=False, repeat=True)
    #uncomment to save animation (takes more time)
    #line_ani.save('hw3movie.mp4', writer="ffmpeg", bitrate=9000)
    plt.show()

    return line_ani



if __name__ == '__main__':
    #The code here should call the performance and analyze routines above and generate the
    #figures that you are submitting with your codes
    performance1()
    performance2()
    analyze1()
    analyze2()
    visualize()
