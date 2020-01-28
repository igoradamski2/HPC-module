"""Final project, part 1
Igor Adamski, CID 01069152
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from m2 import part2 #assumes p2.f90 has been compiled with f2py into module, m2.xxx.so
import timeit



def sim_animation(m,nt,s0,L,d,a):
    x,y,psi = part2.simulate(m,nt,s0,L,d,a)
    fig = plt.figure(figsize=(8,6))
    fig.suptitle('Animation of the pedestrian model with parameters m={},nt={},s0={},L={},d={} \n and alpha={}. sim_animation(), Igor Adamski'.format(m,nt,s0,L,d,a))
    ax = plt.axes(xlim=(-25, 25), ylim=(-25, 25))
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    line, = ax.plot([], [],'o')
    def animate(i):
        data_x = x[:,i-1]
        data_y = y[:,i-1]
        line.set_data(data_x, data_y)
        return line,
    anim = animation.FuncAnimation(fig, animate,
                                  frames=nt, interval=20, blit=False, repeat=True)


    plt.show()

    return anim


def analyze1():
    """Analyze behavior of Psi(t) for no-noise pedestrian model
        Add input/output variables as needed

    ANALYSIS:

    This function generates two figures, first one is a plot of psi(t) for varying parameters
    m and d and the second one is an animation of the network. Lets first look at the plot of
    psi(t). For the case when m=10, we can see that there is a big shoot in the value of psi in the
    first time. The shoot is bigger the lower the value of d is. This is because the pedestrians
    start walking on a grid each close to each other and each with a different random velocity.
    The fact that the lower the d the higher the shoot should be an effect of randomness in this particular
    situation because the first velocities (so psi(2)) do not depend on the parameter d. We can see
    that as t increases, psi(t) stabilizes and stabilizes wuicker for higher values of d. This is expected
    as when d is high then more pedestrains start being in 'close proximity' and more pedestrains
    start walking in the same direction (the model works so that if two people get close enough,
    at the distance d, then they start walking in the same direction). Threfore the higher the d, the
    faster the velocities stabilize and become the same for most of the nodes, and hence psi converges.
    We can see that for d=5 and all the m's, psi(t) converges very quickly to 1. psi(t) = 1 means that
    there is no variation in the directions of velocities in the pedestian model and that everyone is walking
    in the same direction. So theoretically with enough walking we should eventually see that every psi
    converges to 1. For m=10 we see that the smaller d's do not converge to 1, they converge to some value
    below 1, stay there, but then as time progress the psi increases in some steps. Fo the case m=20,
    we now have more points so statistically speaking its more likely for points to bump to each other
    on a LxL grid. This is visible on the graph where again for d=5 psi converges almost instantly to 1, but now
    also psi for d=3 converges very quickly to 1 and stays there. d=1 and d=0.5 still do not converge but
    converge to a value a bit lower than 1 and stay there. If we increased the timesteps enough we would
    see both of them climb to 1. For the case m=50 we have a similar situation but both d=1 and d=0.5 converges
    closer to 1 than for the m=20 case. We also have to make the observation that as m increases, the shoot
    up in the first timestep for psi is lower. This is because as we have more nodes, the directions of
    their initial velocities, which are random, should kind of cancel out when the number of m increases.
    We also can see that for smaller m and smaller d, the psi curves are more rough, and that is the effect of
    step-like increases because the m is small.
    Summarizing, the higher the m, the more synchronized the nodes become with time, and the same applied to d.
    """
    em = [10,20,50]
    di = [0.5, 1, 2, 5]

    psi = []
    i = 0
    j = 0
    for m in em:
        j = 0
        for d in di:
            *b,dummy = part2.simulate(m,300,0.1,25,d,0)
            psi.append(dummy)
            j +=1
        i += 1

    fig2 = plt.figure(figsize=(13,7))
    fig2.suptitle('Figures showing $\psi(t)$ for varying parameters m and d. \n analyze1(), Igor Adamski')
    ax2 = fig2.add_subplot(131)
    ax3 = fig2.add_subplot(132)
    ax4 = fig2.add_subplot(133)
    ax2.plot(psi[0], label='d=0.5,')
    ax2.plot(psi[1], label='d=1')
    ax2.plot(psi[2], label='d=2')
    ax2.plot(psi[3], label='d=5')
    ax2.set_xlabel('t')
    ax2.set_ylabel('$\psi(t)$')
    ax2.set_title('m=10')

    ax3.plot(psi[4], label='d=0.5')
    ax3.plot(psi[5], label='d=1')
    ax3.plot(psi[6], label='d=2')
    ax3.plot(psi[7], label='d=5')
    ax3.set_xlabel('t')
    ax3.set_ylabel('$\psi(t)$')
    ax3.set_title('m=20')

    ax4.plot(psi[8], label='d=0.5')
    ax4.plot(psi[9], label='d=1')
    ax4.plot(psi[10], label='d=2')
    ax4.plot(psi[11], label='d=5')
    ax4.set_xlabel('t')
    ax4.set_ylabel('$\psi(t)$')
    ax4.set_title('m=50')

    ax3.legend()
    ax2.legend()
    ax4.legend()

    anim = sim_animation(10,500,0.1,25,1,0)

    plt.show()

    return anim

def analyze2():
    """Analyze influence of noise on Psi
        Add input/output variables as needed

    ANALYSIS:

    This function like the function before outputs three figures, two animations and plots.
    Lets start with the plots. On the figure we can see three plots that correspond to different
    values of m. On each plot we can see psi(t) for a network with nt=400 timesteps, s0=0.1, L=25
    and d=2, with varying parameter of noise. We can see on all the plots a common theme which is
    that the higher the alpha the lower psi oscillates. In other words we can see that for all the noises
    psi oscillates, but the higher that noise, the lower the value psi oscillates about. This is a suspected
    behavior because as noise gets close to 1, the angle of the velocity is mostly dominated by the random value
    dependent on the noise. So for high noise we get a situation in which none of the nodes actually moves in
    any direction, but all of them oscillate in one place. The lower the alpha the less this effect is visible.
    For example for alpha=0.25 we can see that psi oscillates (m=10) but oscillates around a value close to 1.
    So we have a situation that the nodes move in one direction constantly but just oscillate around as well,
    but the overall direction is kind of constant. We can see that on the first animation, when the nodes oscillate
    highly around themselves but maintain an overall direction of travel as a group. We can also see on the second
    animation, for noise alpha=1, we have that the psi oscillates around a value close to 0 (so generally we shouldnt
    see any aggregate direction of travel) and the nodes seem like moving nowehere as a group, just oscillate individually.

    We can also see a trend as we increase m. As we increase the number of walkers the oscillations of psi are a lot higher.
    They also look like they are a bit periodic after a certain amount of t has passed. The psi for m=30 is very unstable
    and oscillates between high values of psi and very low values of psi. This means that the more walkers we introduce
    the more chaotic the behaviour of their aggregate direction of travel will be. This is because as there are more nodes,
    and remember that d=2 for all the plots, then the random direction of the neighbours affects greately the direction of
    velocity of each node. So since we have more nodes the main direction is more affected and thats why we can see such oscillations
    in psi.
    """

    a = [0, 0.25, 0.5, 0.75, 1]
    psi1 = []
    psi2 = []
    psi3 = []

    for alpha in a:
        *b,dummy = part2.simulate(10,400,0.1,25,2,alpha)
        psi1.append(dummy)

    for alpha in a:
        *b,dummy = part2.simulate(20,400,0.1,25,2,alpha)
        psi2.append(dummy)

    for alpha in a:
        *b,dummy = part2.simulate(30,400,0.1,25,2,alpha)
        psi3.append(dummy)

    fig = plt.figure(figsize=(12,8))
    fig.suptitle(r'Figures showing $\psi(t)$ for varying parameters m and $\alpha$. \n analyze2(), Igor Adamski')
    ax = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)
    ax.plot(psi1[0], label=r'$\alpha=0$')
    ax.plot(psi1[1], label=r'$\alpha=0.25$')
    ax.plot(psi1[2], label=r'$\alpha=0.5$')
    ax.plot(psi1[3], label=r'$\alpha=0.75$')
    ax.plot(psi1[4], label=r'$\alpha=1$')
    ax.set_ylabel('$\psi(t)$')
    ax.set_xticklabels([])
    ax.set_title('m = 10')

    ax2.plot(psi2[0], label=r'$\alpha=0$')
    ax2.plot(psi2[1], label=r'$\alpha=0.25$')
    ax2.plot(psi2[2], label=r'$\alpha=0.5$')
    ax2.plot(psi2[3], label=r'$\alpha=0.75$')
    ax2.plot(psi2[4], label=r'$\alpha=1$')
    ax2.set_ylabel('$\psi(t)$')
    ax2.set_xticklabels([])
    ax2.set_title('m = 20')

    ax3.plot(psi3[0], label=r'$\alpha=0$')
    ax3.plot(psi3[1], label=r'$\alpha=0.25$')
    ax3.plot(psi3[2], label=r'$\alpha=0.5$')
    ax3.plot(psi3[3], label=r'$\alpha=0.75$')
    ax3.plot(psi3[4], label=r'$\alpha=1$')
    ax3.set_xlabel('t')
    ax3.set_ylabel('$\psi(t)$')
    ax3.set_title('m = 30')

    ax.legend()
    ax2.legend()
    ax3.legend()

    anim1 = sim_animation(20,400,0.1,25,2,0.25)
    anim2 = sim_animation(20,400,0.1,25,2,1)
    plt.show()


    return anim1,anim2

def performance():
    """Analyze performance of simulate_omp
    Add input/output variables as needed

    ANALYSIS:
    Techically the question doesnt want us to analyze this.
    But the figure is quite self explanatory. Fortran with OMP is
    2 times faster (cant check for more cores as i only have 2) for m big.
    For m small, there is little or no speedup, as we waste more time on communication
    than we gain on parelallizing the loops.


    Q4:

    To estimate the number of floating point operations (flops) we first count the number
    of total operations we do (+,-,*,/,sqrt,^). In my code I first use 2 flops to build the
    vector diff (two subtractions), then I use 3 flops to bulid the vector vector_norm (two
    multiplications and one addition) and then in the end in an if statement i'm taking a square root
    of the vector_norm so thats an additional plot. Also all the variables i mentioned above
    are arrays of size m**2 so we multiply the number of flops by this number. Also all that calculation
    is inside a do loop which goes through all the m**2 walkers - hence the total number of flops for the distance
    calculations in this code is 6*(m**4). Therefore we can say that distance calculations take up most of
    the operations in the simulate code. This explains why the OMP speedup was twice, because OMP paralellizes
    the do loop over m**2 walkers so that we actually have m**2/numprocs operations.

    CELL-LISTS

    To improve the distance calculations it would be better to use the cell-lists algorithm. This algorithm takes
    the whole space and divides it into a grid of squares of sides d. Then when we loop over all the walkers,
    we would take a walker in a cell and search only for walkers within a distance d of him in the cell and the
    neighbouring cells. This is a lot better approach than the one we used in this code. In my code i basically calculate
    the distance to any walker, which is a massive waste because its obvious that some walkers will be too far. With the
    cell lists algorithm we should expect a big speedup of our code. The cost of using such a method is obviously
    that it would be much more involved to code, especially in fortran.
    """

    em = [10, 20, 30, 40, 50]
    time_f90 = np.zeros([len(em)])
    time_omp = np.zeros([len(em)])
    i=0
    for m in em:
        def sim_f90():
            part2.simulate(m,200,0.1,25,0.5,0)
        def sim_omp():
            part2.simulate_omp(m,200,0.1,25,0.5,0,2)
        time_f90[i] = timeit.timeit(sim_f90, number=3)/float(3)
        time_omp[i] = timeit.timeit(sim_omp, number=3)/float(3)
        i += 1

    x = np.arange(len(em))
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(111)
    ax.set_xticks(x)
    ax.set_xticklabels(em)
    ax.bar(np.array(x)-0.1, time_f90, width=0.2, color='b', align='center', label='Unparalellized fortran')
    ax.bar(np.array(x)+0.1, time_omp, width=0.2, color='r', align='center', label='Paralellized with OMP')
    for i, v in enumerate(time_f90):
        ax.text(i+0.79-1, v+0.01, "{0:.1f}".format(v), color='black')
    for i, v in enumerate(time_omp):
        ax.text(i+0.99-1, v+0.01, "{0:.1f}".format(v), color='black')
    ax.set_ylabel('Time')
    ax.set_title('Time (mean of 3 runs) to simulate walkers for nt=200,s0=0.1,L=25,d=0.5 and alpha=0 \n performance(), Igor Adamski')

    plt.show()



if __name__=='__main__':
    #Modify input/output as needed to generate final figures with call(s) below
    analyze1()
    analyze2()
    performance()
