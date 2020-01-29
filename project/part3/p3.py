"""Final project, part 3
Igor Adamski, CID 01069152
"""
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from f_p3 import part3 #assumes rwnetmod.f90 and p31.f90 have been compiled with f2py into module, m3.xxx.so
from f_p3 import rwnetmod
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation


def gen_net(c, dt):

    G = nx.Graph()
    part3.c = c
    energies = part3.wave(100,2,100,dt,200)
    nodes = part3.xnew.copy()
    edges = part3.e.copy()
    amplitudes = part3.u.copy()

    num = np.arange(len(nodes[0,:])) + 1
    for i in num:
        G.add_node(i, pos=(nodes[0,i-1], nodes[1,i-1]))

    G.add_edges_from(edges)
    #nx.draw_networkx(G, pos = nx.get_node_attributes(G, 'pos'),with_labels=False)

    #compute degree and find max degree node
    deg = G.degree()
    max_deg = max(deg.values())
    #initialize empty node list. we will add nodes with highest degre to it
    node_list = []
    for node in deg.keys():
        if deg[node] == max_deg:
            node_list.append(node)

    #now check which one of these has the shortest path
    length = []
    for node in node_list:
        try:
            p = nx.shortest_path(G, source=num[-1], target=node)
            length.append(len(p))
        except:
            node_list.remove(node)

    connected_to_n = []
    if len(length) == 0:
        for n in num:
            #n is the current node we are checking
            if nx.has_path(G, source=num[-1], target=n) is True:
                connected_to_n.append(n)

        max_deg = 0
        connected_to_n.remove(num[-1])
        for n in connected_to_n:
            if max_deg <= deg[n]:
                curr_n = n
                max_deg = deg[n]

        closest_node = curr_n

    else:
        closest_node = node_list[length.index(min(length))]


    path = nx.shortest_path(G, source=num[-1], target=closest_node)
    #compute things for the third plot
    times = []
    for nds in path:
        for i in range(len(amplitudes[0,:])):
            if amplitudes[nds-1,i] > 0:
                times.append(i)
                break

    return amplitudes, nodes, closest_node, path, times


def wave_analyze():
    """ analyze effectiveness of global b-d method


    ANALYSIS:

    First, a disclaimer. Running this function will produce a different network each time
    hence it is not very easy to produce a analysis basing on plots that generally come out
    of this function. Instead I will focus my analysis on the general trends that one can
    see when running this function multiple times (at least 10 times). Running this function
    can sometimes produce weird error that 'there is no shortest path' but try again and it will
    definitely work (its to do with the try operation). This function produces three figures:
    a animation of the wave propagation on a network when c=1 and two figures showing plots which I
    will describe below, one for c=1 and the other for c=0.5.

    On the animation we can eyeball that the amplitude is greatest for the nth node, so at the one
    from which we are starting the wave generation, which is not very surprisng as at the initial condition
    its the only node with non zero amplitude.

    On the plot figures, the first plot is a plot of the amplitude in time, for the nth node,
    the highest degree node (or when the highest degree node is not connected to the nth node
    the highest degree node becomes the node with the highest degree connected to n) and for a node
    halfway through the path from node n and the highest degree node. For some networks where the highest
    degree node lies close to the nth node, we can almost see a perfect sync of the amplitudes which is
    understandable. For these scenarios (when the highest degree node is around 1-2 nodes away from the
    nth node) we have that the wave decays a lot quicker. You can see that on one of the plots that I attached
    named 'wave_decay.png' where the amplitudes on both nodes steadily go to zero as time increases. This doesnt
    happen when the highest node is not that close to the nth node, it happens for a lot bigger time that I am not
    plotting (im only plotting for 100 timesteps with dt=1).
    We generally have that the amplitudes are highest for the nth node and the further we get
    from the nth node the smaller the amplitudes become.

    On the second plot in the plot figure we can see the time it takes for the wave to reach a node, which is
    on the path from nth node to the highest degree node. An interesting observation is that everytime a wave reaches
    a node it also reaches the next node. So the wave spreads through nodes two at a time and that is connected to
    the way these networks tend to build, that the degree of a node is most of the time 1 or 2. And therefore as
    the wave reaches one node it will automatically affect the other.

    When we vary c, the things get a bit complicated. When c is greater than 1, the wave generation blows up.
    I have not plotted this, but basically what happens is that the highest degree node gets extremely huge values on
    its amplitude. This is an effect of the fact that Aij when the highest degree node is i, is going to contain
    lots of ones, therefore making the sum bigger and also we multiply the sum by c^2 later, making the sum very big.
    For the case that is plotted we have c=0.5, and we generally we can see that the speed of the node is
    twice less than normally. So if we call this function enough times we will see that the timestep until the highest
    degree node has amplitude 0 (so wave still hasnt reached it) is generally twice that of the timestep in a wave with c=1.

    """

    #case when c=1
    amplitudes,nodes,closest_node,path,times = gen_net(1, 1)

    x = np.arange(len(path)) + 1

    fig2 = plt.figure(figsize=(12,8))
    fig2.suptitle('Plots of amplitudes and times of wave reaching a node for c=1. \n wave_analyze(), Igor Adamski.')
    ax2 = fig2.add_subplot(121)
    ax2.plot(amplitudes[-1,:], label='Last node (n)')
    ax2.plot(amplitudes[closest_node-1,:], label='Highest degree node')
    ax2.plot(amplitudes[path[int(len(path)/2)]-1,:], label='Node in between')
    ax2.set_title('Amplitudes in time, for chosen nodes \n in the path from n to the max deg. node')
    ax2.legend()

    ax3 = fig2.add_subplot(122)
    ax3.bar(np.array(x), times, width=0.2, color='b', align='center')
    ax3.set_xticks(x)
    ax3.set_xlabel('Node number on the path from n to max deg. node')
    ax3.set_ylabel('Time for the wave to reach a node')
    ax3.set_title('First time that a wave reaches a node, shown for nodes on the path from \n node n to max deg node.')

    fig = plt.figure(figsize=(10,8))
    fig.suptitle('Animation showing the wave propagation on a 2d network. \n 3d coorcinate represents the wave magnitude. wave_analyze(), Igor Adamski')
    ax = p3.Axes3D(fig)
    ax.set_xlim3d([min(nodes[0,:]), max(nodes[0,:])])
    ax.set_ylim3d([min(nodes[1,:]), max(nodes[1,:])])
    ax.set_zlim3d([-2, 2])
    line1 = ax.plot(nodes[0,:], nodes[1,:], amplitudes[:,0], 'o', color='red')[0]
    line_temp1 = ax.plot([],[],[], 'o', color='blue')[0]

    def animate(i):
        line1.set_data(nodes[:,1:len(nodes[0,:])-1])
        line1.set_3d_properties(amplitudes[1:len(amplitudes[:,1])-1,i])
        line_temp1.set_data(nodes[:,-1])
        line_temp1.set_3d_properties(amplitudes[-1,i])
        return line1

    line_ani1 = animation.FuncAnimation(fig, animate, frames=len(amplitudes[1,:]),interval=60, blit=False, repeat=True)



    #case when c=0.5
    amplitudes,nodes,closest_node,path,times = gen_net(0.5, 1)

    x = np.arange(len(path)) + 1

    fig3 = plt.figure(figsize=(12,8))
    fig3.suptitle('Plots of amplitudes and times of wave reaching a node for c=0.5. \n wave_analyze(), Igor Adamski.')
    ax2 = fig3.add_subplot(121)
    ax2.plot(amplitudes[-1,:], label='Last node (n)')
    ax2.plot(amplitudes[closest_node-1,:], label='Highest degree node')
    ax2.plot(amplitudes[path[int(len(path)/2)]-1,:], label='Node in between')
    ax2.set_title('Amplitudes in time, for chosen nodes \n in the path from n to the max deg. node')
    ax2.legend()

    ax3 = fig3.add_subplot(122)
    ax3.bar(np.array(x), times, width=0.2, color='b', align='center')
    ax3.set_xticks(x)
    ax3.set_xlabel('Node number on the path from n to max deg. node')
    ax3.set_ylabel('Time for the wave to reach a node')
    ax3.set_title('First time that a wave reaches a node, shown for nodes on the path from \n node n to max deg node.')


    plt.show()


    return line_ani1



if __name__=='__main__':
    #Modify input/output as needed to generate final figures with call(s) below
    wave_analyze()
