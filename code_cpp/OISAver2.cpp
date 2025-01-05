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
#include <queue>

using namespace std;

// Function to calculate the Local Clustering Coefficient (LCC) for a node
double computeLCC(const unordered_map<int, unordered_set<int>>& graph, int node) {
    cout << "Calculating LCC for node " << node << endl;
    if (graph.find(node) == graph.end()) {
        cerr << "Error: Node " << node << " does not exist in the graph." << endl;
        return 0.0; // Avoid out_of_range error
    }
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
    double result = static_cast<double>(actualEdges) / possibleEdges;
    cout << "Node " << node << " LCC: " << result << endl;
    return result;
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

// Updated PONF function to use omega_b, omega_c, and omega_d
vector<pair<int, int>> ponf(unordered_map<int, unordered_set<int>>& graph, const vector<int>& targets, int k, double omega_b, double omega_c, int omega_d) {
    vector<pair<int, int>> interventionEdges;

    for (int i = 0; i < k; ++i) {
        int targetNode = -1;
        double maxLCC = -1.0;

        // Find the target node with the maximum LCC
        for (int node : targets) {
            if (graph.find(node) == graph.end()) continue; // Avoid invalid nodes
            double lcc = computeLCC(graph, node);
            if (lcc > maxLCC) {
                maxLCC = lcc;
                targetNode = node;
            }
        }

        if (targetNode == -1) {
            cerr << "Error: No valid target node found in iteration " << i + 1 << endl;
            break;
        }

        cout << "Iteration " << i + 1 << ": Target Node = " << targetNode << ", Max LCC = " << maxLCC << endl;

        int bestNode = -1;
        double bestLCCReduction = 0.0;

        // Find the best candidate node to connect to the target node
        for (const auto& entry : graph) {
            int candidate = entry.first;
            if (candidate == targetNode || graph[targetNode].count(candidate)) continue;

            // Temporarily add the edge
            graph[targetNode].insert(candidate);
            graph[candidate].insert(targetNode);

            double newLCC = computeLCC(graph, targetNode);
            double lccReduction = maxLCC - newLCC;

            double candidateBetweenness = computeBetweenness(graph, candidate);
            double candidateCloseness = computeCloseness(graph, candidate);
            int candidateDegree = computeDegree(graph, candidate);

            // Remove the edge
            graph[targetNode].erase(candidate);
            graph[candidate].erase(targetNode);

            cout << "Candidate: " << candidate << ", LCC Reduction: " << lccReduction
                 << ", Betweenness: " << candidateBetweenness
                 << ", Closeness: " << candidateCloseness
                 << ", Degree: " << candidateDegree << endl;

            if (lccReduction > bestLCCReduction &&
                candidateBetweenness >= omega_b &&
                candidateCloseness >= omega_c &&
                candidateDegree >= omega_d) {
                bestLCCReduction = lccReduction;
                bestNode = candidate;
            }
        }

        if (bestNode != -1) {
            graph[targetNode].insert(bestNode);
            graph[bestNode].insert(targetNode);
            interventionEdges.emplace_back(targetNode, bestNode);
            cout << "Selected Edge: (" << targetNode << ", " << bestNode << ")" << endl;
        } else {
            cout << "No suitable candidate found in this iteration." << endl;
        }
    }

    return interventionEdges;
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
        cout << "u: " << u << ", v: " << v << endl;
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
    string inputFile = "../data/test/3980.txt";
    string outputFile = "../data/test/output.txt";

    unordered_map<int, unordered_set<int>> graph;
    vector<int> targets;
    int k;
    double tau, omega_b, omega_c, omega_d;

    // Read input
    readInputFile(inputFile, graph, targets, k, tau, omega_b, omega_c, omega_d);

    // Compute intervention edges
    vector<pair<int, int>> interventionEdges = ponf(graph, targets, k, omega_b, omega_c, omega_d);

    // Compute the final maximum LCC among targets
    double maxLCC = 0.0;
    for (int node : targets) {
        maxLCC = max(maxLCC, computeLCC(graph, node));
    }

    // Write output
    writeOutputFile(outputFile, interventionEdges, maxLCC);

    return 0;
}
