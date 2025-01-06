#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <queue>
#include <iomanip>
#include <string>
#include <climits>

using namespace std;

// Function to calculate Local Clustering Coefficient (LCC)
double computeLCC(const unordered_map<int, unordered_set<int>>& graph, int node) {
    if (graph.find(node) == graph.end()) return 0.0;

    const auto& neighbors = graph.at(node);
    int degree = neighbors.size();
    if (degree < 2) return 0.0;

    int actualEdges = 0;
    for (auto it1 = neighbors.begin(); it1 != neighbors.end(); ++it1) {
        for (auto it2 = next(it1); it2 != neighbors.end(); ++it2) {
            if (graph.at(*it1).count(*it2)) actualEdges++;
        }
    }

    int possibleEdges = degree * (degree - 1) / 2;
    return static_cast<double>(actualEdges) / possibleEdges;
}

// Function to calculate LCC upper bound
// Equation from 5.3 (ALC)
// double computeLCCUpperBound(int degree, int nv, int k1, int k2) {
//     return (nv + k2 + k1 * degree + k1 * (k1 - 1) / 2) / static_cast<double>((degree + k1) * (degree + k1 - 1) / 2);
// }

double computeLCCUpperBoundOptimized(int degree, int nv, int k) {
    double maxUpperBound = 0.0;

    for (int k1 = 0; k1 <= k; ++k1) {
        int k2 = k - k1;
        double upperBound = (nv + k2 + k1 * degree + k1 * (k1 - 1) / 2) / 
                            static_cast<double>((degree + k1) * (degree + k1 - 1) / 2);
        maxUpperBound = max(maxUpperBound, upperBound);
    }

    return maxUpperBound;
}

// Example implementation for hop distance calculation
int getHopDistance(const unordered_map<int, unordered_set<int>>& graph, int start, int end) {
    unordered_map<int, int> dist;
    queue<int> q;
    q.push(start);
    dist[start] = 0;

    while (!q.empty()) {
        int curr = q.front();
        q.pop();

        for (int neighbor : graph.at(curr)) {
            if (dist.find(neighbor) == dist.end()) {
                dist[neighbor] = dist[curr] + 1;
                if (neighbor == end) {
                    return dist[neighbor];
                }
                q.push(neighbor);
            }
        }
    }
    return INT_MAX; // Return a large value if no path is found
}

// Function to calculate optionality for PONF
int calculateOptionality(const unordered_map<int, unordered_set<int>>& graph, int node, double tau) {
    int optionality = 0;
    if (graph.find(node) == graph.end()) return optionality;

    for (const auto& neighbor : graph.at(node)) {
        // Check if the LCC of the neighbor is within the allowed range
        if (computeLCC(graph, neighbor) <= tau) {
            // Additional hop check (assume a function getHopDistance exists)
            if (getHopDistance(graph, node, neighbor) > 2) {
                optionality++;
            }
        }
    }
    return optionality;
}



// Function to calculate Betweenness Centrality (omega_b)
double computeBetweenness(const unordered_map<int, unordered_set<int>>& graph, int node) {
    cout << "Calculating Betweenness for node " << node << endl;
    unordered_map<int, double> betweenness;
    double betweennessScore = 0.0;

    for (const auto& srcEntry : graph) {
        int src = srcEntry.first;
        if (src == node) continue;

        queue<int> q;
        unordered_map<int, int> dist;
        unordered_map<int, int> paths;
        unordered_map<int, double> delta;

        q.push(src);
        dist[src] = 0;
        paths[src] = 1;

        vector<int> stack;

        while (!q.empty()) {
            int curr = q.front();
            q.pop();
            stack.push_back(curr);

            for (int neighbor : graph.at(curr)) {
                if (dist.find(neighbor) == dist.end()) {
                    dist[neighbor] = dist[curr] + 1;
                    q.push(neighbor);
                }
                if (dist[neighbor] == dist[curr] + 1) {
                    paths[neighbor] += paths[curr];
                }
            }
        }

        while (!stack.empty()) {
            int curr = stack.back();
            stack.pop_back();

            for (int neighbor : graph.at(curr)) {
                if (dist[neighbor] == dist[curr] + 1) {
                    delta[curr] += (static_cast<double>(paths[curr]) / paths[neighbor]) * (1.0 + delta[neighbor]);
                }
            }

            if (curr != src && curr != node) {
                betweennessScore += delta[curr];
            }
        }
    }

    cout << "Node " << node << " Betweenness: " << betweennessScore << endl;
    return betweennessScore;
}

// Function to calculate Closeness Centrality (omega_c)
double computeCloseness(const unordered_map<int, unordered_set<int>>& graph, int node) {
    cout << "Calculating Closeness for node " << node << endl;
    if (graph.find(node) == graph.end()) {
        cerr << "Error: Node " << node << " does not exist in the graph." << endl;
        return 0.0; // Avoid out_of_range error
    }

    queue<int> q;
    unordered_map<int, int> dist;
    int totalDistance = 0;

    q.push(node);
    dist[node] = 0;

    while (!q.empty()) {
        int curr = q.front();
        q.pop();

        for (int neighbor : graph.at(curr)) {
            if (dist.find(neighbor) == dist.end()) {
                dist[neighbor] = dist[curr] + 1;
                totalDistance += dist[neighbor];
                q.push(neighbor);
            }
        }
    }

    if (totalDistance == 0) return 0.0;
    double result = static_cast<double>(dist.size() - 1) / totalDistance;
    cout << "Node " << node << " Closeness: " << result << endl;
    return result;
}

// Function to calculate Degree (omega_d)
int computeDegree(const unordered_map<int, unordered_set<int>>& graph, int node) {
    cout << "Calculating Degree for node " << node << endl;
    if (graph.find(node) == graph.end()) {
        cerr << "Error: Node " << node << " does not exist in the graph." << endl;
        return 0; // Avoid out_of_range error
    }
    int result = graph.at(node).size();
    cout << "Node " << node << " Degree: " << result << endl;
    return result;
}


// Function to implement OISA with EORE, PONF, and ALC
vector<pair<int, int>> oisaAlgorithm(unordered_map<int, unordered_set<int>>& graph, const vector<int>& targets, int k, double tau, double omega_b, double omega_c, double omega_d) {
    vector<pair<int, int>> interventionEdges;

    // Step 1: Precompute LCC upper bounds for all nodes
    unordered_map<int, double> lccUpperBounds;
    for (const auto& entry : graph) {
        int node = entry.first;
        int degree = graph[node].size();
        int nv = computeLCC(graph, node) * (degree * (degree - 1) / 2);
        // lccUpperBounds[node] = computeLCCUpperBound(degree, nv, 1, k - 1);
        lccUpperBounds[node] = computeLCCUpperBoundOptimized(degree, nv, k);

    }

    for (int i = 0; i < k; i++) {
        cout << "Iteration " << i + 1 << ":\n";

        int targetNode = -1;
        double maxLCC = -1.0;

        // Step 2: Find the target node with the highest LCC
        cout << "Finding Target Node...\n";
        for (int node : targets) {
            if (graph.find(node) == graph.end()) continue;
            double lcc = computeLCC(graph, node);
            cout << "Node " << node << " LCC: " << lcc << endl;
            if (lcc > maxLCC) {
                maxLCC = lcc;
                targetNode = node;
            }
        }

        if (targetNode == -1) {
            cout << "No valid target node found. Ending algorithm.\n";
            break;
        }

        cout << "Target Node: " << targetNode << ", Max LCC: " << maxLCC << endl;

        // Step 3: Find the best candidate node to connect to the target node
        int bestNode = -1;
        double bestLCCReduction = 0.0;
        int bestOptionality = INT_MAX;
        double bestMissScore = -1.0;

        cout << "Evaluating Candidate Nodes...\n";
        for (const auto& entry : graph) {
            int candidate = entry.first;
            cout << "candidate: " << candidate << endl;
            if (candidate == targetNode || graph[targetNode].count(candidate)) continue;

            // Skip if LCC upper bound indicates no significant reduction
            if (lccUpperBounds[candidate] > tau) {
                cout << "lccUpperBounds > tau" << endl;   
                continue;
            }
            // Temporarily add the edge
            graph[targetNode].insert(candidate);
            graph[candidate].insert(targetNode);

            double newLCC = computeLCC(graph, targetNode);
            double lccReduction = maxLCC - newLCC;
            int candidateOptionality = calculateOptionality(graph, candidate, tau);

            double candidateBetweenness = computeBetweenness(graph, candidate);
            double candidateCloseness = computeCloseness(graph, candidate);
            int candidateDegree = computeDegree(graph, candidate);

            // Calculate MISS based on the formula
            double missScore = (omega_b * candidateBetweenness) + 
                               (omega_c * candidateCloseness) + 
                               ((1 - omega_b - omega_c) * candidateDegree);

            // Remove the edge
            graph[targetNode].erase(candidate);
            graph[candidate].erase(targetNode);

            // Compare candidates based on LCC reduction, optionality, and MISS
            cout << "Candidate: " << candidate 
                 << ", LCC Reduction: " << lccReduction 
                 << ", Optionality: " << candidateOptionality 
                 << ", MISS Score: " << missScore << endl;

            if (lccReduction > bestLCCReduction || 
                (lccReduction == bestLCCReduction && candidateOptionality < bestOptionality) || 
                (lccReduction == bestLCCReduction && candidateOptionality == bestOptionality && missScore > bestMissScore)) {
                bestLCCReduction = lccReduction;
                bestNode = candidate;
                bestOptionality = candidateOptionality;
                bestMissScore = missScore;
            }
        }

        if (bestNode != -1) {
            cout << "Best Candidate: " << bestNode 
                 << ", Best LCC Reduction: " << bestLCCReduction 
                 << ", Best Optionality: " << bestOptionality 
                 << ", Best MISS Score: " << bestMissScore << endl;

            graph[targetNode].insert(bestNode);
            graph[bestNode].insert(targetNode);
            interventionEdges.emplace_back(targetNode, bestNode);

            cout << "Added Intervention Edge: (" << targetNode << ", " << bestNode << ")\n";
        } else {
            cout << "No valid candidate found for Target Node: " << targetNode << endl;
        }

        cout << "End of Iteration " << i + 1 << "\n\n";
    }

    // Final summary
    cout << "Final Intervention Edges:\n";
    for (const auto& edge : interventionEdges) {
        cout << edge.first << " " << edge.second << endl;
    }

    return interventionEdges;
}


// Function to read input from file
void readInputFile(const string& filename, unordered_map<int, unordered_set<int>>& graph, vector<int>& targets, int& k, double& tau, double& omega_b, double& omega_c, double& omega_d) {
    ifstream file(filename);
    if (!file) {
        cerr << "Error: Unable to open input file." << endl;
        exit(1);
    }

    int n, m;
    file >> n >> m;

    int t;
    file >> t;
    targets.resize(t);
    for (int i = 0; i < t; ++i) {
        file >> targets[i];
    }

    file >> k >> tau >> omega_b >> omega_c >> omega_d;

    for (int i = 0; i < m; ++i) {
        int u, v;
        file >> u >> v;
        graph[u].insert(v);
        graph[v].insert(u);
    }
}

// Function to write output to file
void writeOutputFile(const string& filename, const vector<pair<int, int>>& interventionEdges, double maxLCC) {
    ofstream file(filename);
    if (!file) {
        cerr << "Error: Unable to open output file." << endl;
        exit(1);
    }

    file << interventionEdges.size() << "\n";
    for (const auto& edge : interventionEdges) {
        file << edge.first << " " << edge.second << "\n";
    }
    file << fixed << setprecision(6) << maxLCC << "\n";
}

int main() {
    string inputFile = "../data/example/in.txt";
    string outputFile = "../data/example/output_ver3.txt";

    unordered_map<int, unordered_set<int>> graph;
    vector<int> targets;
    int k;
    double tau, omega_b, omega_c, omega_d;

    // Read input
    readInputFile(inputFile, graph, targets, k, tau, omega_b, omega_c, omega_d);

    // Compute intervention edges
    vector<pair<int, int>> interventionEdges = oisaAlgorithm(graph, targets, k, tau, omega_b, omega_c, omega_d);

    // Compute the final maximum LCC among targets
    double maxLCC = 0.0;
    for (int node : targets) {
        maxLCC = max(maxLCC, computeLCC(graph, node));
    }

    // Write output
    writeOutputFile(outputFile, interventionEdges, maxLCC);

    return 0;
}
