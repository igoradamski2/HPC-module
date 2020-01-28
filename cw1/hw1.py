"""Igor Adamski
CID: 01069152
"""
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import scipy.spatial as spatial
from collections import Counter


def rw2d(Nt,M,a=0,b=0):
    """Input variables
    Nt: number of time steps
    M: Number of simulations
    a,b: bias parameters
    """
    assert a >= -1, "ValueError: a must be greater or equal to 0"
    assert b <= 1, "ValueError: b must be less or equal to 1"
    #simulate the steps
    Xsteps = np.random.randint(0,2,[M,Nt+1])
    Ysteps = np.random.randint(0,2,[M,Nt+1])
    #change the values to proper steps
    Xsteps[Xsteps == 1] = -1
    Xsteps[Xsteps == 0] = 1+a
    Ysteps[Ysteps == 1] = -1
    Ysteps[Ysteps == 0] = 1-b
    #convention is we start at 0
    Xsteps[:,0] = 0.
    Ysteps[:,0] = 0.
    #calculate on-the-go position from steps
    X = np.cumsum(Xsteps, axis = 1)
    Y = np.cumsum(Ysteps, axis = 1)
    #calculate all needed statistics
    Xm = np.mean(X, axis = 0)
    Ym = np.mean(Y, axis = 0)
    Xv = np.var(X, axis = 0)
    Yv = np.var(Y, axis = 0)
    Xv = Xv + Xm ** 2
    Yv = Yv + Ym ** 2
    XY = np.mean(np.multiply(X,Y), axis = 0)

    return Xm,Ym,Xv,Yv,XY


def rwnet1(H,Hf,a=0,display=False):
    """Input variables
    H: Height at which new nodes are initially introduced
    Hf: Final network height
    a: horizontal bias parameter, should be -1, 0, or 1
    display: figure displaying the network is created when true
    Output variables
    X,Y: Final node coordinates
    positions: a matrix giving the coordinates of the positions of the nodes
    """

    assert a in [-1,0,1], "ValueError: horizontal bias parameter a should be -1, 0, or 1"

    #initialize the graph, positions vector, distance and add the initial node
    G = nx.Graph()
    G.add_node(1, pos=(0,0))
    positions = np.zeros((1,2))
    index = 1
    dstar = np.sqrt(1+(1+a)**2)
    while positions[-1][1] < Hf:
        index += 1
        #place a point at (0,H) before the beginning of walk
        temp_pos = np.array([0,H])
        #make a tree of points in positions vector
        point_tree = spatial.cKDTree(positions)
        #check the y-coordinate of the highest point in the tree
        ver_dist = max(positions[:,1])
        #vectorize steps (using 2* because on average that will be enough)
        steps_x, steps_y, *b = rw2d(int(2*(H-ver_dist)), 1, a=a, b=1)
        steps_y += H
        #if any of the steps in y-direction are greater than the ball around the highest point
        if np.any(steps_y > (ver_dist + dstar)):
            #index of the furthest value
            idx = np.where(steps_y > (ver_dist + dstar))[0][-1]
            temp_pos = np.array([steps_x[idx], steps_y[idx]])
        else:
            temp_pos = temp_pos

        while True:
            #check which nodes in the network are within dstar of current node
            indeces = point_tree.query_ball_point(temp_pos, np.sqrt(1+(1+a)**2))
            if len(indeces) != 0 or temp_pos[1] <= 0:
                break
            #if none nodes are within dstar continue walking
            x_move = np.random.choice([1+a, -1])
            y_move = np.random.choice([0, -1])
            temp_pos[0] += x_move
            temp_pos[1] += y_move

        G.add_node(index, pos=(temp_pos[0], temp_pos[1]))
        positions = np.vstack((positions, temp_pos))
        #add edges between connected nodes
        for idx in indeces:
            G.add_edge(index, int(idx)+1)

    #this part is to make the figure look better
    if len(positions) >= 200:
        with_labels = False
        node_size = 40
    else:
        with_labels = True
        node_size = 100

    #display the figure of the network if True
    if display is True:
        plt.figure(figsize=(12,8))
        nx.draw_networkx(G, pos = nx.get_node_attributes(G, 'pos'), with_labels=with_labels, node_size=node_size, width = 1.5, font_size = 7)
        plt.title('Network for a={} created by rwnet1(H={},Hf={},a={}). Igor Adamski'.format(a, H, Hf, a))
        plt.axis('on')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()

    #extract final output from vector of final positions
    X = positions[:,0]
    Y = positions[:,1]

    return X,Y,positions

def rwnet2(L,H,Hf,a=0,display=False):
    """Input variables
    L: Walls are placed at X = +/- L
    H: Height at which new nodes are initially introduced
    Hf: Final network height
    a: horizontal bias parameter, should be -1, 0, or 1
    display: figure displaying the network is created when true
    Output variables
    X,Y: Final node coordinates
    positions: a matrix giving the coordinates of the positions of the nodes
    """

    assert a in [-1,0,1], "ValueError: horizontal bias parameter a should be -1, 0, or 1"

    G = nx.Graph()
    G.add_node(1, pos=(0,0))
    positions = np.zeros((1,2))
    index = 1
    dstar = np.sqrt(1+(1+a)**2)
    while positions[-1][1] < Hf:
        index += 1
        temp_pos = np.array([0, H])
        point_tree = spatial.cKDTree(positions)
        ver_dist = max(positions[:,1])
        #ensure we dont do unnecessary steps beyond wall
        num_steps = min(int(2*(H-ver_dist)), L)
        steps_x, steps_y, *b = rw2d(num_steps, 1, a=a, b=1)
        steps_y += H

        if np.any(steps_y > (ver_dist + dstar)):
            idx1 = np.where(steps_y > (ver_dist + dstar))[0][-1]
            #index the first time when x gets to L
            if np.any(abs(steps_x) >= L):
                idx2 = np.where(abs(steps_x) >= L)[0][-1]
            else:
                idx2 = idx1
            idx = min(idx1, idx2)
            temp_pos = np.array([steps_x[idx], steps_y[idx]])
        else:
            temp_pos = temp_pos

        while True:
            indeces = point_tree.query_ball_point(temp_pos, np.sqrt(1+(1+a)**2))
            #get the sign of the x-coordinate of our point
            sgn = np.sign(temp_pos[0])
            if len(indeces) != 0 or temp_pos[1] <= 0:
                break
            x_move = np.random.choice([1+a, -1])
            y_move = np.random.choice([0, -1])
            temp_pos[0] += x_move
            #check if we are not attempting to cross the wall
            if abs(temp_pos[0]) > L:
                #place node at wall if above happens
                temp_pos[0] = sgn * L
            temp_pos[1] += y_move

        G.add_node(index, pos=(temp_pos[0], temp_pos[1]))
        positions = np.vstack((positions, temp_pos))
        for idx in indeces:
            G.add_edge(index, int(idx)+1)

    if len(positions) >= 500:
        with_labels = False
        node_size = 40
    else:
        with_labels = True
        node_size = 100

    if display is True:
        plt.figure(figsize=(12,8))
        nx.draw_networkx(G, pos = nx.get_node_attributes(G, 'pos'), with_labels=with_labels, node_size=node_size, width = 1.5, font_size = 7)
        plt.title('Network for a = {}, L = {} created by rwnet2(L={},H={},Hf={},a={}). Igor Adamski'.format(a,L,L,H,Hf,a))
        plt.axis('on')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()

    X = positions[:,0]
    Y = positions[:,1]

    return X,Y,positions

def cumulative_rate(array):
    """Input variables:
    array: A 1-dimensional array
    Output variables:
    A 2 dimensional array of values of the array and indexes at which they were
    first achieved in the array
    """
    #for the purpose of question 3 we start at the second value of array
    curr_max = array[1]
    cumulative_rate = np.array([curr_max, 1])
    for i in range(2,len(array)):
        if array[i] > curr_max:
            curr_max = array[i]
            cumulative_rate = np.vstack((cumulative_rate, [curr_max, i]))
        else:
            continue

    return cumulative_rate


def analyze(H = 200, Hf = 150, display1 = True, display2 = True):
    """ Figure 1: represents a histogram of rates of height increase for different
        values of parameters a and L.
        The rate is calculated as ∆height/(∆ in no. of nodes)=Hf/no.of nodes.
        We can see that for a=1 and L=30 the rate is the highest, which makes sense because
        the network will build very fast as the nodes are biased to the right and when they reach
        the wall they quickly fall down making almost a vertical strip of nodes. The second fastest
        is a=1 L=150 and it also makes sense as since we start at (0,200) and the expected step in
        x direction is +1, then we are quite sure to reach the wall, where the nodes accumulate
        similarly to L=30. However, of course, they have to travel a lot further and can spread from each other
        on the way thats why the rate is significantly lower than with L=30.
        We can see that for a=0 the rates are very similar which is not surprising because there is no
        bias towards any x-direction and the nodes just fall down, accumulating and creating a tower.
        This is supposed to happen as fast with the walls and without them so its not surprising.
        A surprising trend we see that for a no-wall scenario, the rate is higher for a=1 than for a=0,
        suggesting that for a=0 the tower is buildt which is wide and for a=1 the tower accumulates
        height faster making a kind of arm shape reaching from the far x-directions to (0,Hf).

        Figure 2: represents a plot of the biggest height at different times of
        the network generation for fixed a=0 and different values of L.
        It is calculated as simply checking the networks height at different 'times' of the
        network generation so at different number of nodes present in the network.
        We can see that the curve increases steadily and very similarly for a=0. This confirms
        what we have found in the histogram that the rate of increase is roughly the same
        for a=0. We can see that for different L's the network reaches final height for slightly
        different 'times', but that could be just the effect of randomness as there is no
        significant difference in the rate at which the height increases.

        Figure 3: represents a plot of the biggest height at different times of
        the network generation for fixed a=1 and different values of L.
        It is calculated in the same way as above, so we check the biggest height in the network
        at different number of nodes present in the network.
        We can see a pattern unfold here. For L=30 the curve steeps very fast towards Hf,
        a lot faster than for the other L's. This was also the case in the histogram. For the other
        L's we see that for L=150 it reaches Hf faster than with L=inf, which was also a trend
        shown by the histogram. What is an interesting trend that was missed by the 'general' approach
        of the histogram is that the rate of increase of the height of the network for L=30 and L=150
        is roughly the same up until we have added around 100 nodes. After that point, the height of the
        network increases a lot slower for L=150 and the slope doesnt change for L=30.

        The fact that for a=1 the rates of height increase are generally higher than for a=0
        could be also explained by the fact that dstar for a=1 is greater than dstar for a=0
        hence we have a greater range of reaching a node, hence the network will build quicker.

    """
    #create all the random networks we are going to analyze
    X0, Y0, positions_a0 = rwnet1(H,Hf, a = 0)
    X1, Y1, positions_a1 = rwnet1(H,Hf, a = 1)
    X2, Y2, positions_a0_l30 = rwnet2(30, H, Hf, a = 0)
    X3, Y3, positions_a1_l30 = rwnet2(30, H, Hf, a = 1)
    X4, Y4, positions_a0_l150 = rwnet2(150, H, Hf, a = 0)
    X5, Y5, positions_a1_l150 = rwnet2(150, H, Hf, a = 1)

    if display1 is True:
        rate_1 = Hf/float(len(X0))
        rate_2 = Hf/float(len(X1))
        rate_3 = Hf/float(len(X2))
        rate_4 = Hf/float(len(X3))
        rate_5 = Hf/float(len(X4))
        rate_6 = Hf/float(len(X5))

        rate = [rate_3, rate_5, rate_1, rate_4, rate_6, rate_2]

        plt.figure(figsize=(7,6))
        x = np.arange(6)
        plt.bar(x, height = rate)
        plt.xticks(x, ['a=0 \n L=30','a=0 \n L=150','a=0 \n L=inf','a=1 \n L=30','a=1 \n L=150','a=1 \n L=inf'])
        plt.title("Figure1 : Histogram showing the rate of height increase \n with respect to the node growth \n (generated by analyze(). Igor Adamski)")
        plt.xlabel("Different combiantions of parameters a and L")
        plt.ylabel(r"$\frac{Hf}{no. of nodes}$" + "   ")
        #plt.show(block=False)

    if display2 is True:

        rate_a0 = cumulative_rate(Y0)
        rate_a0_l30 = cumulative_rate(Y2)
        rate_a0_l150 = cumulative_rate(Y4)

        rate_a1 = cumulative_rate(Y1)
        rate_a1_l30 = cumulative_rate(Y3)
        rate_a1_l150 = cumulative_rate(Y5)

        plt.figure(figsize=(7,6))
        plt.plot(rate_a0[:,1], rate_a0[:,0], 'r', label="a=0, L=inf")
        #plt.hold(True)
        plt.plot(rate_a0_l30[:,1], rate_a0_l30[:,0], 'b', label="a=0, L=30")
        plt.plot(rate_a0_l150[:,1], rate_a0_l150[:,0], 'g', label="a=0, L=150")
        plt.legend()
        plt.xlabel("Number of nodes present in the network")
        plt.ylabel("Height of the graph")
        plt.title("Figure 2 : Plot showing the increase of network height \n with respect to the increase of the number of nodes \n (generated by analyze(). Igor Adamski)")
        #plt.show(block=False)

        plt.figure(figsize=(7,6))
        plt.plot(rate_a1[:,1], rate_a1[:,0], 'r', label="a=1, L=inf")
        #plt.hold(True)
        plt.plot(rate_a1_l30[:,1], rate_a1_l30[:,0], 'b', label="a=1, L=30")
        plt.plot(rate_a1_l150[:,1], rate_a1_l150[:,0], 'g', label="a=1, L=150")
        plt.legend()
        plt.xlabel("Number of nodes present in the network")
        plt.ylabel("Height of the graph")
        plt.title("Figure 3 : Plot showing the increase of network height \n with respect to the increase of the number of nodes \n (generated by analyze(). Igor Adamski)")
        plt.show()

    return None

def network(X,Y,dstar,display=False,degree=False):
    """ Input variables
    X,Y: Numpy arrays containing network node coordinates
    dstar: Links are placed between nodes within a distance, d<=dstar of each other
        and dstar2 = dstar*dstar
    display: Draw graph when true
    degree: Compute, display and return degree distribution for graph when true
    Output variables:
    G: NetworkX graph corresponding to X,Y,dstar
    D: degree distribution, only returned when degree is true

    Degree histogram generation:

    The histogram is generated in the following way. A function degree from
    the networkx package calculates the dictionary with node numbers as keys
    and their respective degrees as values. Then using Counter from collections
    the frequency of each degree is calculated, creating essentially a distribution
    of degrees. Then to make the histogram all the missing degrees (values of D)
    are added to D so that we can create a nice looking histogram and see the
    degree distribution.

    Weakness of network model:

    (I assume we are talking about the networks generated by rwnet1 and rwnet2)
    In this model for the network, nodes are dropped from above and are let to
    randomly walk downwards until they reach another node or nodes within a
    certain distance. Therefore there is relatively small possibility for the
    nodes to cluster next to each other. If by some chance nodes were to cluster into
    some shape, then a falling node would not attach to all the nodes in the shape
    but would only attach to the nodes on one side of the boundary of that shape.
    It is the fact that nodes stop walking until they reach a proximity of a SINGLE
    node that prevents the nodes from clustering and making many connections and hence
    making the maximum degree constrained.
    Modification: We could modify the model in the following way to avoid the constraint
    explained above. A very easy modification that could have been done is to walk until
    our node is connected to a fixed number of nodes. For example, we could walk until our node is
    in proximity of at least 3-4 other nodes and only then connect it. Of course at first
    we would not get many connections and all the nodes will be placed at y=0 value and we wouldnt get far.
    Therefore, we could not only place a initial node at (0,0) but place a whole lattice of points
    for example place nodes at (i,j) for i,j in range(4). Then we would have assurance that the network would
    reach crazy big degrees. To change that we would just need to seed the initial nodes at these positions
    and drop a node like we did before. Also we would need to change the constraint in line 85
    from len(indeces) != 0 to len(indeces) >= 3 which would ensure that we stop the walk when a node is in
    proximity of 3 nodes.

    Another idea is to just allow steps in non-discrete manner. Doing it discretely we constrain that one node
    can be connected to 8 nodes at once, because of the grid they are moving on, the grid being in steps of one.
    If we scaled our steps to be a bit more continuous (so instead of step +1 do step +0.5). This would generally
    have the same effect as increasing dstar and increasing Hf.
    """
    #initialize graph
    G = nx.Graph()
    mat = np.array([X,Y])
    positions = np.transpose(mat)
    #add all nodes at their respective positions
    for i in range(1,len(X)+1):
        G.add_node(i, pos = (positions[i-1][0], positions[i-1][1]))

    #create a tree of points
    point_tree = spatial.cKDTree(positions)
    index = 0
    for points in positions:
        index += 1
        indeces = point_tree.query_ball_point(points, dstar)
        #if node has no connections then remove it
        if len(indeces) == 1:
            G.remove_node(index)
        else:
            for idx in indeces:
                G.add_edge(index, idx+1)

    if display is True:
        plt.figure(figsize=(12,8))
        nx.draw_networkx(G, pos = nx.get_node_attributes(G, 'pos'), with_labels = False, node_size =50, width = 1.5, font_size = 6)
        plt.axis('on')
        plt.title('Network generated by network(X=X, Y=Y, dstar={}). Igor Adamski'.format(dstar))
        plt.xlabel("x")
        plt.ylabel("y")
        plt.show()

    D = None
    if degree is True:
        deg = G.degree()
        #D will give the frequency of occurence of a degree in D, hence the distribution
        D = Counter(deg.values())
        plt.figure(figsize=(7,6))
        x = list(range(1, max(list(deg.values()))+1))
        #add all missing in-between values to D
        for i in x:
            if i not in list(D.keys()):
                D[i] = 0
        height = []
        #sorting the D.values() in order of ascending keys
        for keys in sorted(D.keys()):
            height.append(D[keys])

        plt.bar(x, height = height)
        plt.xlabel("Degree")
        plt.ylabel("Frequency of that degree in the graph")
        plt.title('The degree distribution for a graph \n generated by network(X=X, Y=Y, dstar={}). Igor Adamski'.format(dstar))
        plt.show()



    return G, D




if __name__ == '__main__':
    analyze()
