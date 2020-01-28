"""Assumes cost.f90 has been compiled with f2py to generate the module,
hw2mod.so (filename may also be of form hw2mod.xxx.so where xxx is system-dependent text) using
IGOR ADAMSKI: 01069152
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from hw2mod import cost
from hw2mod import hw2
import scipy.optimize as scopt
import time


def visualize(Nx,Ny):
    """Display cost function with and without noise on an Ny x Nx grid
    """
    #create grid
    x_1 = np.linspace(-10, 10, Nx)
    x_2 = np.linspace(-10, 10, Ny)
    points = np.array(np.meshgrid(x_1,x_2)).T.reshape(-1,2)
    cost_z = []
    cost_noise = []
    cost.c_noise = False
    #append cost to vector for each point
    for i in range(len(points)):
        cost_z.append(cost.costj(points[i]))

    #plot
    fig1 = plt.figure(figsize=(10,8))
    fig1.suptitle('Surface and contour plots of the cost function without noise. \n Called by visualize(1000,1000), Igor Adamski')
    ax = fig1.add_subplot(211, projection='3d')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('cost')
    ax2 = fig1.add_subplot(212)
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    X, Y = np.meshgrid(x_1,x_2)
    Z = np.array([cost_z]).reshape(X.shape)
    ax.plot_surface(X,Y,Z)
    ax2.contour(X,Y,Z)
    plt.show()

    #noise initialization
    cost.c_noise = True
    cost.c_noise_amp = 1.
    for i in range(len(points)):
        cost_noise.append(cost.costj(points[i]))
    #plot
    fig2 = plt.figure(figsize=(10,8))
    fig2.suptitle('Surface and contour plots of the cost function with added noise 1.\n Called by visualize(1000,1000), Igor Adamski')
    ax = fig2.add_subplot(211, projection='3d')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('cost')
    ax2 = fig2.add_subplot(212)
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    Z_new = np.array([cost_noise]).reshape(X.shape)
    ax.plot_surface(X,Y,Z_new)
    ax2.contour(X,Y,Z_new)
    plt.show()


def newton_test(xg, display=False):
    """ Use Newton's method to minimize cost function defined in cost module
    Input variable xg is initial guess for location of minimum. When display
    is true, a figure illustrating the convergence of the method should be
    generated

    Output variables: xf -- computed location of minimum, jf -- computed minimum
    Further output can be added to the tuple, output, as needed. It may also
    be left empty.
    """
    #call newton and extract paths
    cost.c_noise = False
    nw = hw2.newton(xg)
    xf = nw[0]
    jf = nw[1]
    jpath_nw = hw2.jpath.copy()
    xpath_nw = hw2.xpath.copy()
    #get the norm distance to [1,1] from each point
    norm_distance = np.linalg.norm((xpath_nw - [1,1]), axis=1)

    #plot
    if display is True:
        fig1 = plt.figure(figsize=(11,7))

        x = list(range(len(jpath_nw)))
        ax = fig1.add_subplot(121)
        ax2 = fig1.add_subplot(122)

        ax.plot(jpath_nw, label='Cost')
        ax.set_xlabel('Iteration number')
        ax.set_ylabel('Value of the cost function')
        ax.set_title('Value of the cost function \n against Newton iterations starting at \n xguess={}.\n Igor Adamski'.format(xg,xg))
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax.grid(linestyle=':', linewidth=0.5)
        ax.set_yscale('log')
        ax.set_xticks(x)

        ax2.plot(norm_distance, label='Distance from [1,1]')
        ax2.set_xlabel('Iteration number')
        ax2.set_ylabel('Norm distance to the minimum at [1,1]')
        ax2.set_title('Distance to the actual minimum at [1,1] \n at each Newton iteration, starting from \n xguess={}.\n Called by newton_test({}, True)  Igor Adamski'.format(xg,xg))
        ax2.set_xticks(x)
        ax2.grid(linestyle=':', linewidth=0.5)

        plt.show()

    return xf,jf


def bracket_descent_test(xg,display=False):
    """ Use bracket-descent to minimize cost function defined in cost module
    Input variable xg is initial guess for location of minimum. When display
    is true, 1-2 figures comparing the B-D and Newton steps should be generated

    Output variables: xf -- computed location of minimum, jf -- computed minimum


    On the first figure we can see the paths that the Newton and bracket descent
    algorithms take to finally converge to the minimum. From that figure we can see
    that Newton tends to do very little steps (at most 6-7) and the steps that it takes
    are big. By big, I mean that Newton takes huge steps in directions that are far from
    the final converged minimum but always manages to return back. On the other hand
    bracket descent does lots of very small steps. The bracket descent steps are small
    but ultimately lead it to the minimum. We can also see that bracket descent doesnt
    do jumps in opposite directions to which its supposed to go in contrast to Newton.

    On the second figure we can see the normed distance between the current point
    and the true minimum at [1,1] plotted against the iterations. We can see that
    Newton at first jumps very far from the actual minimum, only to 2 steps later converge
    to [1,1]. On the other hand, bracket descent tends to stay on the trajectory going down
    except a few times when it goes further from [1,1] for a bit only to return promptly.
    Overall, its visible that newton takes a lot less (a factor of hundreds) iterations to
    converge than bracket descent.

    The behaviour of Newton is not surprising in terms of the number of iterations,
    as it relies on the gradient and hessian to guide it the right way. Gradient and
    hessian give the information about the curvature of the cost surface and this
    makes newton always follow the right track. ON the other hand, bracket descent does not
    possess such high level information about the cost function and it simply iterates
    and finds smaller values of the cost function and proceeds to take little steps adjusting
    the triangles vertices. Bracket descent takes small steps because by construction its only
    allowed to make changes to the vertices of at most 4 heights of the traingle. Newton takes
    big steps because the gradient and hessian can lead him far away only to return to the minimum in the next step.
    """
    #call bracket descent and extract paths
    cost.c_noise = False
    bd = hw2.bracket_descent(xg)
    xf = bd[0]
    jf = bd[1]
    xpath_bd = hw2.xpath.copy()
    #computed normed distance to [1,1]
    norm_distance_bd = np.linalg.norm((xpath_bd - [1,1]), axis=1)

    #purely esthetics of plot
    if np.linalg.norm(xg)<10:
        scale="linear"
    else:
        scale="symlog"

    #plot
    if display is True:
        fig1 = plt.figure()
        fig2 = plt.figure()

        nw = hw2.newton(xg)
        xpath_nw = hw2.xpath
        norm_distance_nw = np.linalg.norm((xpath_nw - [1,1]), axis=1)

        ax = fig1.add_subplot(111)
        ax.plot(xpath_bd[:,0], xpath_bd[:,1], 'r', label='Bracket descent', linewidth = 1)
        ax.scatter(xpath_bd[:,0], xpath_bd[:,1], s=10, c='black')
        ax.scatter(xg[0], xg[1], c='green', label='Initial guess')
        ax.plot(xpath_nw[:,0], xpath_nw[:,1], 'b', label='Newton')
        #ax.plot(xf[0], xf[1],
                #'c*', label='Final convergence of bracket_descent')
        ax.scatter(nw[0][0], nw[0][1],
                c='orange', marker='D', label='Actual minimum')
        ax.set_yscale(scale)
        ax.set_ylabel('y-coordinate')
        ax.set_xlabel('x-coordinate')
        ax.set_title('The paths of bracket descent and Newton \n starting from xguess={}. \n Cmd. bracket_descent_test({},True), Igor Adamski'.format(xg, xg))
        ax.legend()

        ax2 = fig2.add_subplot(111)
        ax2.plot(norm_distance_nw, label='Newton')
        ax2.plot(norm_distance_bd, label='Bracket descent')
        ax2.set_xlabel('Iterations')
        ax2.set_ylabel('Norm distance from [1,1]')
        ax2.set_yscale(scale)
        ax2.set_title('Norm distance from minimum at [1,1] against \n the number of iterations. \n Cmd. bracket_descent_test({},True). Igor Adamski'.format(xg))
        ax2.legend()

        plt.show()

    return xf, jf

#this was once used to get xpath from lbfgs. Can be added to minimize
#function as callback=callback
cb = []
def callback(xk):
    global cb
    cb.append(xk)

def performance(noise=False):
    """ Assess performance of B-D and L-BFGS-B methods. Add input/output as
    needed

    (i) NO-NOISE:

    On the first figure we can see a barplot of number of iterations it takes
    for the l-bfgs-l and bracket descent to converge for different starting points.
    The clear trend here is that the further from [1,1] we are the more iterations
    it takes to converge for both algorithms. We can also see that bracket descent
    takes twice as more iterations when starting at [-1,-3] and the gap grows as we
    decrease the x coordinate of the xguess. In the above questions we have investigated
    more deeply how bracket descent works and its not surprising that it takes a high
    number of iterations. L-bfgs-l is a quasi newton method and approximates the hessian
    so acts kindof like the newton method in that it takes less steps. However it still takes
    more steps than Newton as it only approximates the hessian and does not really have it.

    On the second figure we can see the final values that the cost function takes at the
    points where the two algorithms have converged. We immediately see that the final cost
    values are a lot lower (by a factor of 1e6) for the l-bfgs-l than for bracket descent.
    This is the second performance edge of l-bfgs-l over bracket_descent. It converged to a
    lot lower. However, the final value of the bracket descent is still very small, and close to
    zero. The huge difference in the cost values of l-bfgs-l and BD are due to the fact that near the
    minimum [1,1] the function changes dramatically in the marginal values (going to 0 at [1,1]). Still
    the overall trend we see is that l-bfgs-l converges to a better minimum.

    On the third graph we can see the averge time it takes for the l-bfgs-l and bracket_descent
    to converge to the minimum. The average is taken over 10 consecutive runs. The time is given
    in microseconds (second * 1e-6) so we can see that both algorithms have very fast times. However,
    the bracket descent is significantly faster compared to l-bfgs-l. In fact, bracket descent is
    34 times faster starting from the closest point and the margin increases and for the furthest point
    bracket descent is 40 times faster. This is an interesting performance gain and edge of bracket descent
    over l-bfgs-l. The reason why bracket descent is so much faster is that it does a lot of very small
    and not very heavy on memory computations. The code is written in fortran so the loops are properly optimized
    and the only operations bracket descent is doing is multiplication and division and a bunch of if statements.
    Such computations are a lot quicker than approximating the functions hessian and solving linear equations,
    which is done by l-bfgs-l.

    Overall, we can see that in terms of time bracket descent is a lot faster but needs more iterations to converge
    and converges to weaker minimums. Also, bracket descent does not converge for some values (try [40,240]).

    (ii) NOISE:

    On the first graph we can see the final values of the cost function for the two optimizers
    when the function has had noise=100 introduced. We can see that now the sides flip (in relation to the
    case with no noise), because brackt descent seems to converge to a lower minimum or similar than l-bfgs-l.
    This could be the effect of the fact that l-bfgs-l, because its using a approximation of the hessien,
    has its steps accurately calculated, and is somewhat thrown off when it lands on a spot where noise was added. On the
    other hand bracket descent simply follows the lowering cost and doesnt care and think forward about its walk.

    On the second graph we can see that again the further the starting point from [1,1] the longer
    it takes to converge (logical). We can see again that bracket descent is significantly faster than l-bfgs-l but now
    the margin is smaller. Now bracket descent is 6 times faster for the closest point and 10 times faster for the furthest,
    which comparing to the values without noise on cost is a significant reduction in speed edge. We can see that
    adding cost has slowed bracket descent more than it slowed l-bfgs-l. This might be the effect of ravines and tunnels created
    by the noise, through which the bracket descent had to manouvre through and since it takes more small steps than l-bfgs-l
    it was more affected by that.

    (iii) IMPROVEMENTS:

    One obvious improvement we could add to our code is we could generalize it to work on different cost functions
    rather than just working on the cost provided from the cost.f90 module. This should not be too hard to implement
    we would just have to provide the function as a input. Then this work would have a lot more usability, because someone
    could for example optimize a function which eh wishes to optimize using our module. Then using our python code he could
    compare his results using the newton , bracket descent and l-bfgs-l methods and pick one which suits him best. For that sake
    we could modify our above code so that the user would have the choice of the initial guess, and the graphs and figures would
    produce themselves accordingly. Another improvement we could make is we could add more optimize to which we would compare newton and bracket
    descent. This way a scientifically literate person could find a comrehensive review of the distinct features that
    different optimization gives him and pick the one that works best for him.

    """

    if noise is False:
        cost.c_noise = False
        #initialize matrix of points
        xg_mat = [[-1,-3], [-10,-3], [-50,-3], [-100,-3]]

        #all necessary declarations
        nit_bd = np.zeros([4])
        nit_lbfgs = np.zeros([4])
        idx = 0
        bracket_costs = np.zeros([4])
        lbfgs_costs = np.zeros([4])
        #create figures
        fig1 = plt.figure(figsize=(7,6))
        fig2 = plt.figure(figsize=(7,6))
        fig3 = plt.figure(figsize=(7,6))
        ax = fig1.add_subplot(111)
        ax2 = fig3.add_subplot(111)
        ax3 = fig2.add_subplot(211)
        ax4 = fig2.add_subplot(212)


        #do time average
        lbfgs_av_times = np.zeros([4])
        bracket_av_times = np.zeros([4])
        for i in range(10):
            bracket_times = np.zeros([4])
            lbfgs_times = np.zeros([4])
            k = 0
            for points in xg_mat:
                st_lbfgs = time.time()
                lbfgs = scopt.minimize(cost.costj, points, method='L-BFGS-B',
                                        options={'gtol': 1e-6, 'maxiter': 1000},
                                        )
                lbfgs_times[k] = 1e6 * (time.time() - st_lbfgs)

                st_bracket = time.time()
                bd = hw2.bracket_descent(points)
                bracket_times[k] = 1e6 * (time.time() - st_bracket)
                k += 1
            lbfgs_av_times += lbfgs_times
            bracket_av_times += bracket_times
        lbfgs_av_times = lbfgs_av_times/10
        bracket_av_times = bracket_av_times/10

        #do other graphs
        for points in xg_mat:

            lbfgs = scopt.minimize(cost.costj, points, method='L-BFGS-B',
                                    options={'gtol': 1e-6, 'maxiter': 1000},
                                    )
            bd = hw2.bracket_descent(points)

            #number_of_iterations
            dummy1 = len(hw2.xpath)
            nit_lbfgs[idx] = lbfgs.nit
            nit_bd[idx] = dummy1

            #final costs
            bracket_costs[idx] = bd[1]
            lbfgs_costs[idx] = lbfgs.fun

            idx+=1

        #plot of number of iterations

        x = list(range(1,5))
        ax.bar(np.array(x)-0.1, nit_bd, width=0.2, color='b', align='center', label='bracket_descent')
        ax.bar(np.array(x)+0.1, nit_lbfgs, width=0.2, color='r', align='center', label='l-bfgs-b')
        ax.set_xticks(x)
        ax.set_xticklabels(xg_mat)
        ax.set_xlabel("Initial xguess")
        ax.set_ylabel("Number of algorithm iterations it takes to converge")
        ax.legend()
        ax.set_title("Bar plot showing the number of iterations \n until convergence for bracket descent and l-bfgs-b. \n Cmd. performance(), Igor Adamski")

        #plot of times

        ax2.bar(np.array(x)-0.1, bracket_av_times, width=0.2, color='b', align='center', label='bracket_descent')
        ax2.bar(np.array(x)+0.1, lbfgs_av_times, width=0.2, color='r', align='center', label='l-bfgs-b')
        ax2.set_xticks(x)
        ax2.set_xticklabels(xg_mat)
        ax2.set_xlabel("Initial xguess")
        ax2.set_ylabel("Average time it takes for the algorithm to converge ($\mu$s)")
        ax2.legend()
        ax2.set_title("Bar plot showing the average time it takes \n bracket descent and l-bfgs-b to converge (mean of 10 runs). \n Cmd. performance(), Igor Adamski")
        for i, v in enumerate(bracket_av_times):
            ax2.text(i+0.79, v+0.01, "{0:.1f}".format(v), color='black')
        for i,v in enumerate(lbfgs_av_times):
            ax2.text(i+0.99, v+0.01, "{0:.1f}".format(v), color='black')

        #plot of final costs

        ax3.bar(np.array(x), bracket_costs, width=0.2, color='b', align='center', label='bracket_descent')
        ax3.set_xticks(x)
        ax3.set_xticklabels(xg_mat)
        ax3.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax4.bar(np.array(x), lbfgs_costs, width=0.2, color='r', align='center', label='l-bfgs-b')
        ax3.set_ylabel('Value of cost')
        ax4.set_xlabel('Different starting xguess')
        ax4.set_xticks(x)
        ax4.set_xticklabels(xg_mat)
        ax4.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax4.set_ylabel('Value of cost')
        ax3.legend()
        ax4.legend()
        fig2.suptitle('Final values of the cost function after convergence \n for bracket descent and l-bfgs-l. \n Cmd. performance(), Igor Adamski')

    elif noise is True:
        cost.c_noise = True
        cost.c_noise_amp = 100.
        xg_mat = [[-1,-3], [-10,-3], [-50,-3], [-100,-3]]

        #all the same initializations as above but with noise
        nit_bd = np.zeros([4])
        nit_lbfgs = np.zeros([4])
        xpath_lbfgs = {}
        idx = 0
        bracket_costs = np.zeros([4])
        lbfgs_costs = np.zeros([4])
        fig2 = plt.figure(figsize=(7,6))
        fig3 = plt.figure(figsize=(7,6))
        ax2 = fig3.add_subplot(111)
        ax3 = fig2.add_subplot(211)
        ax4 = fig2.add_subplot(212)


        #do time average
        lbfgs_av_times = np.zeros([4])
        bracket_av_times = np.zeros([4])
        for i in range(10):
            bracket_times = np.zeros([4])
            lbfgs_times = np.zeros([4])
            k = 0
            for points in xg_mat:
                st_lbfgs = time.time()
                lbfgs = scopt.minimize(cost.costj, points, method='L-BFGS-B',
                                        options={'gtol': 1e-6, 'maxiter': 1000},
                                        )
                lbfgs_times[k] = 1e6 * (time.time() - st_lbfgs)

                st_bracket = time.time()
                bd = hw2.bracket_descent(points)
                bracket_times[k] = 1e6 * (time.time() - st_bracket)
                k += 1
            lbfgs_av_times += lbfgs_times
            bracket_av_times += bracket_times
        lbfgs_av_times = lbfgs_av_times/10
        bracket_av_times = bracket_av_times/10

        #do other graphs
        for points in xg_mat:

            lbfgs = scopt.minimize(cost.costj, points, method='L-BFGS-B',
                                    options={'gtol': 1e-6, 'maxiter': 1000},
                                    )
            bd = hw2.bracket_descent(points)

            #final costs
            bracket_costs[idx] = bd[1]
            lbfgs_costs[idx] = lbfgs.fun

            idx+=1

        #plot of times
        x = list(range(1,5))
        ax2.bar(np.array(x)-0.1, bracket_av_times, width=0.2, color='b', align='center', label='bracket_descent')
        ax2.bar(np.array(x)+0.1, lbfgs_av_times, width=0.2, color='r', align='center', label='l-bfgs-b')
        ax2.set_xticks(x)
        ax2.set_xticklabels(xg_mat)
        ax2.set_xlabel("Initial xguess")
        ax2.set_ylabel("Average time it takes for the algorithm to converge ($\mu$s)")
        ax2.legend()
        ax2.set_title("Bar plot showing the average time it takes bracket descent \n and l-bfgs-b to converge (mean of 10 runs) for cost with noise 100. \n Cmd. performance(True), Igor Adamski")
        for i, v in enumerate(bracket_av_times):
            ax2.text(i+0.79, v+0.01, "{0:.1f}".format(v), color='black')
        for i,v in enumerate(lbfgs_av_times):
            ax2.text(i+0.99, v+0.01, "{0:.1f}".format(v), color='black')

        #plot of final costs

        ax3.bar(np.array(x), bracket_costs, width=0.2, color='b', align='center', label='bracket_descent')
        ax3.set_xticks(x)
        ax3.set_xticklabels(xg_mat)
        ax3.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax3.axhline(0, color="black", linewidth=0.9)
        ax3.set_ylabel('Value of cost')
        ax4.bar(np.array(x), lbfgs_costs, width=0.2, color='r', align='center', label='l-bfgs-b')
        ax4.set_xticks(x)
        ax4.set_xticklabels(xg_mat)
        ax4.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax4.axhline(0, color="black",linewidth=0.9)
        ax4.set_ylabel('Value of cost')
        ax4.set_xlabel('Initial xguess')
        ax3.legend()
        ax4.legend()
        fig2.suptitle('Final values of the cost function after convergence \n for bracket descent and l-bfgs-l with added noise of 100. \n Cmd. performance(True), Igor Adamski')




if __name__ == '__main__':
    #Add code here to call newton_test, bracket_descent_test, performance

    #figures for visualize
    visualize(1000,1000)

    #graphs for Newton
    xmat = [[-0.5,0.5], [10,10], [1203918,-1203123]]
    for points in xmat:
        newton_test(points, display=True)

    #figures for bracket descent test
    bracket_descent_test([-7,5],display=True)

    #figures for performance
    performance()
    performance(noise=True)
