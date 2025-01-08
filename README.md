# Traveling Salesman Problem (TSP) Solver using Quantum Computing

## Problem Description
This implementation solves the Traveling Salesman Problem with the following specifications:

- **Objective**: Minimize total distance traveled
- **Distance Properties**: 
  - Symmetric distances (distance A→B equals B→A)
  - Open-ended tour (final city position doesn't matter)

## Problem Representation

### Distance Matrix
The problem uses a 3x3 distance matrix where each entry (i,j) represents the distance between cities i and j. The matrix is symmetric around the diagonal, with zeros on the diagonal (distance from a city to itself).

### Cost Function
The cost function minimizes the sum of distances between consecutive cities in the path. 

### Solution Encoding
- **Quantum State**: Binary string representing the path through cities
- **Classical Representation**: Arrow notation (e.g., "1->2->0") showing the order of city visits
- **Probability Distribution**: The solver returns multiple possible solutions with their associated probabilities, as shown in the histogram where "1->2->0" has the highest probability (~50%)

The histogram shows that the quantum algorithm finds multiple valid solutions, with varying probabilities reflecting their optimality according to the cost function.
