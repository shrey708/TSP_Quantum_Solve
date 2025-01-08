
import numpy as np
import itertools

def assign_cities(N):
    # creates an array for forming cities random on 10x10 grid
    cities = []
    for i in range(N):
        cities.append(np.random.rand(2)*10)
    return np.array(cities)

def distance_points(point1, point2):
    return np.sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2)

def total_distance_matrix(cities):
    no_cities = len(cities)
    matrix = np.zeros((no_cities, no_cities))
    for i in range(no_cities):
        for j in range (i, no_cities):
            matrix[i][j] = distance_points(cities[i], cities[j])
            matrix[j][i] = matrix[i][j]
    return matrix

def cost_fun(cost_matrix, solution):
    cost = 0 
    # -1, as we dont care to return to the initial starting point
    for i in range(len(solution) - 1 ):
        a = i%len(solution)
        b = (i+1)%len(solution)
        cost = cost  + cost_matrix[solution[a]][solution[b]]
    return cost

def solve_trail_error(distance_matrix, initial_city = None, verbose = True):
    
    no_cities = len(distance_matrix)
    initial_order = range(no_cities)

    all_permutations = [list(x) for x in itertools.permutations(initial_order)]
    best_permuation = all_permutations[0]
    best_cost = cost_fun(distance_matrix, best_permuation)*1000

    for permutation in all_permutations:
        if initial_city:
            if permutation[0] != initial_city:
                continue
        current_cost = cost_fun(distance_matrix, permutation)
        if current_cost < best_cost:
            best_permuation = permutation
            best_cost = cost_fun
        if verbose:
            print("Best route to take: ", best_permuation)
            print("Cost: ", best_cost)
        return best_permuation
    
def points_to_binary (points):
        # converts ordered points to binary 
        no_points = len(points)
        binary_conv = np.zeros((len(points))**2)
        for i in range(len(points)):
            p = points[i]
            binary_conv[(no_points) * (i)+ p[1]] = 1
        return binary_conv

def binary_to_points(binary):
    # converts binary to ordered points
    no_points = int(np.sqrt(len(binary)))
    points = []
    for i in range(no_points):
        for j in range(no_points):
            if binary[i*no_points + j] == 1:
                points.append(j)
    return points

def binary_points_fixed_start(binary_point):

    ordered_points = [0]
    no_points = int(np.sqrt(len(binary_point))+1)
    for i in range(no_points-1):
        for j in range(no_points-1):
            if binary_point[(no_points -1)*i+ j] == 1:
                ordered_points.append(j+1)
    return ordered_points



    
