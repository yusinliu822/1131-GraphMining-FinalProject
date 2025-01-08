// claude

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
#include <stdexcept>

using namespace std;

class Graph {
private:
    unordered_map<int, unordered_set<int>> adjacencyList;
    
    void validateNode(int node) const {
        if (adjacencyList.find(node) == adjacencyList.end()) {
            throw runtime_error("Node " + to_string(node) + " does not exist in the graph.");
        }
    }

public:
    void addEdge(int u, int v) {
        adjacencyList[u].insert(v);
        adjacencyList[v].insert(u);
    }

    void removeEdge(int u, int v) {
        validateNode(u);
        validateNode(v);
        adjacencyList[u].erase(v);
        adjacencyList[v].erase(u);
    }

    bool hasEdge(int u, int v) const {
        auto it = adjacencyList.find(u);
        return it != adjacencyList.end() && it->second.count(v);
    }

    const unordered_set<int>& getNeighbors(int node) const {
        auto it = adjacencyList.find(node);
        if (it == adjacencyList.end()) {
            throw runtime_error("Node not found");
        }
        return it->second;
    }

    double computeLCC(int node) const {
        validateNode(node);
        const auto& neighbors = adjacencyList.at(node);
        int degree = neighbors.size();
        if (degree < 2) return 0.0;

        int actualEdges = 0;
        for (auto it1 = neighbors.begin(); it1 != neighbors.end(); ++it1) {
            for (auto it2 = next(it1); it2 != neighbors.end(); ++it2) {
                if (hasEdge(*it1, *it2)) actualEdges++;
            }
        }

        int possibleEdges = degree * (degree - 1) / 2;
        return static_cast<double>(actualEdges) / possibleEdges;
    }

    double computeLCCUpperBound(int node, int k) const {
        validateNode(node);
        int degree = adjacencyList.at(node).size();
        if (degree <= 1) return 0.0;

        double maxUpperBound = 0.0;
        int nv = computeLCC(node) * (degree * (degree - 1) / 2);

        for (int k1 = 0; k1 <= k; ++k1) {
            int k2 = k - k1;
            if (degree + k1 <= 1) continue;  // Avoid division by zero
            
            double upperBound = (nv + k2 + k1 * degree + k1 * (k1 - 1) / 2) / 
                               static_cast<double>((degree + k1) * (degree + k1 - 1) / 2);
            maxUpperBound = max(maxUpperBound, upperBound);
        }

        return maxUpperBound;
    }

    // Floyd-Warshall algorithm for all-pairs shortest paths
    vector<vector<int>> computeAllPairsShortestPaths() const {
        vector<int> nodes;
        for (const auto& pair : adjacencyList) {
            nodes.push_back(pair.first);
        }
        sort(nodes.begin(), nodes.end());
        
        int n = nodes.size();
        vector<vector<int>> dist(n, vector<int>(n, INT_MAX));
        
        // Initialize distances
        for (int i = 0; i < n; i++) {
            dist[i][i] = 0;
            for (int neighbor : adjacencyList.at(nodes[i])) {
                auto it = find(nodes.begin(), nodes.end(), neighbor);
                int j = distance(nodes.begin(), it);
                dist[i][j] = 1;
            }
        }
        
        // Floyd-Warshall algorithm
        for (int k = 0; k < n; k++) {
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    if (dist[i][k] != INT_MAX && dist[k][j] != INT_MAX) {
                        dist[i][j] = min(dist[i][j], dist[i][k] + dist[k][j]);
                    }
                }
            }
        }
        
        return dist;
    }

    double computeBetweenness(int node) const {
        validateNode(node);
        vector<vector<int>> shortestPaths = computeAllPairsShortestPaths();
        double betweennessScore = 0.0;
        
        vector<int> nodes;
        for (const auto& pair : adjacencyList) {
            nodes.push_back(pair.first);
        }
        sort(nodes.begin(), nodes.end());
        
        int nodeIndex = find(nodes.begin(), nodes.end(), node) - nodes.begin();
        int n = nodes.size();
        
        // For each pair of nodes (s,t)
        for (int s = 0; s < n; s++) {
            for (int t = 0; t < n; t++) {
                if (s == t || s == nodeIndex || t == nodeIndex) continue;
                
                // Count number of shortest paths between s and t
                int totalPaths = 0;
                int pathsThroughNode = 0;
                
                // BFS to count paths
                queue<pair<int, vector<int>>> q;
                q.push({s, {s}});
                
                while (!q.empty()) {
                    auto [curr, path] = q.front();
                    q.pop();
                    
                    if (curr == t) {
                        totalPaths++;
                        // Check if current path goes through our target node
                        if (find(path.begin(), path.end(), nodeIndex) != path.end()) {
                            pathsThroughNode++;
                        }
                        continue;
                    }
                    
                    // Add neighbors that are on shortest path
                    for (int i = 0; i < n; i++) {
                        if (shortestPaths[curr][i] == 1 && 
                            shortestPaths[i][t] == shortestPaths[curr][t] - 1) {
                            vector<int> newPath = path;
                            newPath.push_back(i);
                            q.push({i, newPath});
                        }
                    }
                }
                
                if (totalPaths > 0) {
                    betweennessScore += static_cast<double>(pathsThroughNode) / totalPaths;
                }
            }
        }
        
        return betweennessScore;
    }

    double computeCloseness(int node) const {
        validateNode(node);
        vector<vector<int>> shortestPaths = computeAllPairsShortestPaths();
        
        int totalDistance = 0;
        int reachableNodes = 0;
        
        for (const vector<int>& distances : shortestPaths) {
            for (int dist : distances) {
                if (dist != INT_MAX && dist > 0) {
                    totalDistance += dist;
                    reachableNodes++;
                }
            }
        }
        
        return reachableNodes > 0 ? static_cast<double>(reachableNodes) / totalDistance : 0.0;
    }

    int computeDegree(int node) const {
        validateNode(node);
        return adjacencyList.at(node).size();
    }

    const unordered_map<int, unordered_set<int>>& getGraph() const {
        return adjacencyList;
    }
};

class OISA {
private:
    Graph graph;
    vector<int> targets;
    int k;
    double tau, omega_b, omega_c, omega_d;
    
    bool validateParameters() const {
        if (omega_b + omega_c > 1.0) {
            throw invalid_argument("omega_b + omega_c must be <= 1.0");
        }
        if (tau < 0.0 || tau > 1.0) {
            throw invalid_argument("tau must be between 0.0 and 1.0");
        }
        if (k < 0) {
            throw invalid_argument("k must be non-negative");
        }
        return true;
    }

public:
    OISA(const Graph& g, const vector<int>& t, int edges, double threshold, 
         double wb, double wc, double wd) : 
        graph(g), targets(t), k(edges), tau(threshold), 
        omega_b(wb), omega_c(wc), omega_d(wd) {
        validateParameters();
    }

    vector<pair<int, int>> run() {
        vector<pair<int, int>> interventionEdges;
        
        for (int i = 0; i < k; i++) {
            // Find target node with highest LCC
            int targetNode = -1;
            double maxLCC = -1.0;
            
            for (int node : targets) {
                try {
                    double lcc = graph.computeLCC(node);
                    if (lcc > maxLCC) {
                        maxLCC = lcc;
                        targetNode = node;
                    }
                } catch (const runtime_error& e) {
                    cerr << "Warning: " << e.what() << endl;
                    continue;
                }
            }
            
            if (targetNode == -1) break;
            
            // Find best candidate node
            int bestNode = -1;
            double bestScore = -1.0;
            
            for (const auto& entry : graph.getGraph()) {
                int candidate = entry.first;
                if (candidate == targetNode || graph.hasEdge(targetNode, candidate)) {
                    continue;
                }
                
                // Calculate candidate score
                try {
                    double candidateScore = 
                        (omega_b * graph.computeBetweenness(candidate)) +
                        (omega_c * graph.computeCloseness(candidate)) +
                        ((1 - omega_b - omega_c) * graph.computeDegree(candidate));
                        
                    if (candidateScore > bestScore) {
                        bestScore = candidateScore;
                        bestNode = candidate;
                    }
                } catch (const exception& e) {
                    cerr << "Warning: Error calculating score for node " << candidate << endl;
                    continue;
                }
            }
            
            if (bestNode != -1) {
                graph.addEdge(targetNode, bestNode);
                interventionEdges.emplace_back(targetNode, bestNode);
            }
        }
        
        return interventionEdges;
    }
};

void readInputFile(const string& filename, Graph& graph, vector<int>& targets,
                  int& k, double& tau, double& omega_b, double& omega_c, double& omega_d) {
    ifstream file(filename);
    if (!file) {
        throw runtime_error("Unable to open input file: " + filename);
    }

    int n, m, t;
    file >> n >> m;
    if (file.fail()) {
        throw runtime_error("Error reading n and m from input file");
    }

    file >> t;
    targets.resize(t);
    for (int i = 0; i < t; ++i) {
        file >> targets[i];
    }

    file >> k >> tau >> omega_b >> omega_c >> omega_d;
    if (file.fail()) {
        throw runtime_error("Error reading parameters from input file");
    }

    for (int i = 0; i < m; ++i) {
        int u, v;
        file >> u >> v;
        if (file.fail()) {
            throw runtime_error("Error reading edge " + to_string(i));
        }
        graph.addEdge(u, v);
    }
}

void writeOutputFile(const string& filename, const vector<pair<int, int>>& edges, double maxLCC) {
    ofstream file(filename);
    if (!file) {
        throw runtime_error("Unable to open output file: " + filename);
    }

    file << edges.size() << "\n";
    for (const auto& edge : edges) {
        file << edge.first << " " << edge.second << "\n";
    }
    file << fixed << setprecision(6) << maxLCC << "\n";
}

int main() {
    try {
        string inputFile = "../data/example/in.txt";
        string outputFile = "../data/example/output.txt";

        Graph graph;
        vector<int> targets;
        int k;
        double tau, omega_b, omega_c, omega_d;

        readInputFile(inputFile, graph, targets, k, tau, omega_b, omega_c, omega_d);
        
        OISA oisa(graph, targets, k, tau, omega_b, omega_c, omega_d);
        vector<pair<int, int>> interventionEdges = oisa.run();

        double maxLCC = 0.0;
        for (int node : targets) {
            maxLCC = max(maxLCC, graph.computeLCC(node));
        }

        writeOutputFile(outputFile, interventionEdges, maxLCC);
        
        return 0;
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
}