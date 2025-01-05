// without lower bounds of betweenness, closeness, and degree 𝜔𝑏, 𝜔𝑐, 𝜔d

#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <utility>
#include <iomanip>
#include <string>


using namespace std;

// Function to calculate the Local Clustering Coefficient (LCC) for a node
double computeLCC(const unordered_map<int, unordered_set<int>>& graph, int node) {
    const auto& neighbors = graph.at(node);
    int degree = neighbors.size();
    if (degree < 2) return 0.0;

    int actualEdges = 0;
    for (auto it1 = neighbors.begin(); it1 != neighbors.end(); ++it1) {
        for (auto it2 = next(it1); it2 != neighbors.end(); ++it2) {
            if (graph.at(*it1).count(*it2)) {
                actualEdges++;
            }
        }
    }

    int possibleEdges = degree * (degree - 1) / 2;
    return static_cast<double>(actualEdges) / possibleEdges;
}

// Function to read input from file
void readInputFile(const string& filename, unordered_map<int, unordered_set<int>>& graph, vector<int>& targets,
                    int& k, double& tau, double& omega_b, double& omega_c, double& omega_d) {
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

// Function to implement the PONF strategy and select intervention edges
vector<pair<int, int>> ponf(unordered_map<int, unordered_set<int>>& graph, const vector<int>& targets, int k) {
    vector<pair<int, int>> interventionEdges;

    for (int i = 0; i < k; ++i) {
        int targetNode = -1;
        double maxLCC = -1.0;

        // Find the target node with the maximum LCC
        for (int node : targets) {
            double lcc = computeLCC(graph, node);
            if (lcc > maxLCC) {
                maxLCC = lcc;
                targetNode = node;
            }
        }

        int bestNode = -1;
        double bestLCCReduction = 0.0;

        // Find the best candidate node to connect to the target node
        for (const auto& entry : graph) { int candidate = entry.first;
            if (candidate == targetNode || graph[targetNode].count(candidate)) continue;

            // Temporarily add the edge
            graph[targetNode].insert(candidate);
            graph[candidate].insert(targetNode);

            double newLCC = computeLCC(graph, targetNode);
            double lccReduction = maxLCC - newLCC;

            // Remove the edge
            graph[targetNode].erase(candidate);
            graph[candidate].erase(targetNode);

            if (lccReduction > bestLCCReduction) {
                bestLCCReduction = lccReduction;
                bestNode = candidate;
            }
        }

        if (bestNode != -1) {
            graph[targetNode].insert(bestNode);
            graph[bestNode].insert(targetNode);
            interventionEdges.emplace_back(targetNode, bestNode);
        }
    }

    return interventionEdges;
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
    string outputFile = "../data/example/output.txt";

    unordered_map<int, unordered_set<int>> graph;
    vector<int> targets;
    int k;
    double tau, omega_b, omega_c, omega_d;

    // Read input
    readInputFile(inputFile, graph, targets, k, tau, omega_b, omega_c, omega_d);

    // Compute intervention edges
    auto interventionEdges = ponf(graph, targets, k);

    // Compute the final maximum LCC among targets
    double maxLCC = 0.0;
    for (int node : targets) {
        maxLCC = max(maxLCC, computeLCC(graph, node));
    }

    // Write output
    writeOutputFile(outputFile, interventionEdges, maxLCC);

    return 0;
}
