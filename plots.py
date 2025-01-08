import matplotlib.pyplot as plt
from utilities import cost_fun, total_distance_matrix
import numpy as np

def plt_cities(cities):
    # plots the cities
    plt.scatter(cities[:,0], cities[:,1], s=200)
    for i, city in enumerate(cities):
        plt.annotate(i,(city[0] + 0.15, city[1] + 0.15), size=16, color='r')
    plt.show()

def plt_solution(cities, solution):
    plt_cities(cities)
    plt.xlim([min(cities[:, 0]) - 1, max(cities[:, 0]) + 1])
    plt.ylim([min(cities[:, 1]) - 1, max(cities[:, 1]) + 1])
    for i in range(len(solution) - 1):
        a = i%len(solution)
        b = (i+1)%len(solution)
        A = solution[a]
        B = solution[b]
        plt.plot([cities[A, 0], cities[B, 0]], [cities[A, 1], cities[B, 1]], c='r')
    cost = cost_fun(total_distance_matrix(cities), solution)
    title_string = "Cost:" + str(cost)
    title_string += "\n" + str(solution)
    plt.title(title_string)


def plt_histogram(states_prob):

    states = np.array(states_prob)[:0]
    probs = np.array(states_prob)[:1].astype(float)
    n = len(states_prob)
    plt.barh(range(n), probs, tick_label=states, align='center')
    plt.show()

