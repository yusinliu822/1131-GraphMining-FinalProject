#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <queue>
#include <iomanip>
#include <climits>

using namespace std;

// Function to calculate Local Clustering Coefficient (LCC)
double computeLCC(const unordered_map<int, unordered_set<int>>& graph, int node) {
    if (graph.find(node) == graph.end() || graph.at(node).size() < 2) {
        cout << "Node " << node << " has less than 2 neighbors. LCC = 0.\n";
        return 0.0;
    } 

    const auto& neighbors = graph.at(node);
    int actualEdges = 0;
    for (auto it1 = neighbors.begin(); it1 != neighbors.end(); ++it1) {
        for (auto it2 = next(it1); it2 != neighbors.end(); ++it2) {
            if (graph.at(*it1).count(*it2)) actualEdges++;
        }
    }

    int possibleEdges = neighbors.size() * (neighbors.size() - 1) / 2;

    cout << "Node " << node << ": Neighbors = " << neighbors.size() 
         << ", ActualEdges = " << actualEdges 
         << ", PossibleEdges = " << possibleEdges << "\n";


    return static_cast<double>(actualEdges) / possibleEdges;
}

// Optimized LCC upper bound computation
double computeLCCUpperBound(int degree, int nv, int k) {
    double maxUpperBound = 0.0;

    for (int k1 = 0; k1 <= k; ++k1) {
        int k2 = k - k1;
        if ((degree + k1) * (degree + k1 - 1) / 2 == 0) continue; // 避免分母為 0

        double upperBound = (nv + k2 + k1 * degree + k1 * (k1 - 1) / 2) /
                            static_cast<double>((degree + k1) * (degree + k1 - 1) / 2);

        maxUpperBound = max(maxUpperBound, upperBound);
    }

    return maxUpperBound;
}

// Compute shortest path (used for hop distance)
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
                if (neighbor == end) return dist[neighbor];
                q.push(neighbor);
            }
        }
    }
    return INT_MAX;
}

// Function to calculate optionality (PONF)
int calculateOptionality(const unordered_map<int, unordered_set<int>>& graph, int node, double tau) {
    int optionality = 0;

    if (graph.find(node) == graph.end()) return optionality;

    for (const auto& neighbor : graph.at(node)) {
        double neighborLCC = computeLCC(graph, neighbor);

        if (neighborLCC <= tau) {
            int hopDistance = getHopDistance(graph, node, neighbor);
            if (hopDistance > 2) optionality++;
        }
    }

    return optionality;
}


// Compute Betweenness Centrality
double computeBetweenness(const unordered_map<int, unordered_set<int>>& graph, int node) {
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

    return betweennessScore;
}

// Compute Closeness Centrality
double computeCloseness(const unordered_map<int, unordered_set<int>>& graph, int node) {
    if (graph.find(node) == graph.end()) return 0.0;

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

    if (dist.size() <= 1) return 0.0;
    return static_cast<double>(dist.size() - 1) / totalDistance;
}

// Compute Degree
int computeDegree(const unordered_map<int, unordered_set<int>>& graph, int node) {
    if (graph.find(node) == graph.end()) return 0;
    return graph.at(node).size();
}

// Main OISA Algorithm
vector<pair<int, int>> oisaAlgorithm(unordered_map<int, unordered_set<int>>& graph, 
                                     const vector<int>& targets, int k, double tau, 
                                     double omega_b, double omega_c, double omega_d) {
    vector<pair<int, int>> interventionEdges;

    unordered_map<int, double> lccUpperBounds;
    for (const auto& entry : graph) {
        int node = entry.first;
        int degree = graph[node].size();
        int nv = computeLCC(graph, node) * (degree * (degree - 1) / 2);
        lccUpperBounds[node] = computeLCCUpperBound(degree, nv, k);
    }

    int j = 1; // 初始化迭代次數
    double lj = 0.0; // 初始化 l_j
    double maxLCC = 0.0;

    // 計算初始最大 LCC
    for (int node : targets) {
        maxLCC = max(maxLCC, computeLCC(graph, node));
    }

    int d2k = 2 * k; // 假設與 k（新增邊數）相關的最大度數
    double comb = d2k * (d2k - 1) / 2.0; // 計算 C(d_{2k}, 2)

    while (lj < maxLCC) { // 當 l_j 小於最大 LCC 時繼續
        lj = static_cast<double>(j) / comb; // 計算 l_j
        cout << "Iteration " << j << ": l_j = " << lj << ", Max LCC = " << maxLCC << "\n";

        int targetNode = -1;
        double currentMaxLCC = -1.0;

        // 找到最大 LCC 的目標節點
        for (int node : targets) {
            if (graph.find(node) == graph.end()) continue;
            double lcc = computeLCC(graph, node);
            cout << "Target Node: " << node << ", LCC = " << lcc << "\n";

            if (lcc > currentMaxLCC) {
                currentMaxLCC = lcc;
                targetNode = node;
            }
        }

        if (targetNode == -1) {
            cout << "No valid target node found.\n";
            break;
        }
        cout << "Selected Target Node: " << targetNode << ", Max LCC = " << currentMaxLCC << "\n";

        int bestNode = -1;
        double bestLCCReduction = 0.0;
        int bestOptionality = INT_MAX;
        double bestMISS = -1.0;

        // 評估候選節點
        for (const auto& entry : graph) {
            int candidate = entry.first;
            if (candidate == targetNode || graph[targetNode].count(candidate)) continue;

            // Skip if LCC upper bound indicates no significant reduction
            if (lccUpperBounds[candidate] > tau) {
                cout << "lccUpperBounds > tau" << endl;   
                continue;
            }

            graph[targetNode].insert(candidate);
            graph[candidate].insert(targetNode);
            
            double newLCC = computeLCC(graph, targetNode);
            double lccReduction = currentMaxLCC - newLCC;
            int optionality = calculateOptionality(graph, candidate, tau);

            double betweenness = computeBetweenness(graph, candidate);
            double closeness = computeCloseness(graph, candidate);
            int degree = computeDegree(graph, candidate);

            double miss = (omega_b * betweenness) + (omega_c * closeness) + 
                          ((1 - omega_b - omega_c) * degree);

            graph[targetNode].erase(candidate);
            graph[candidate].erase(targetNode);

            cout << "Candidate: " << candidate 
                 << ", LCC Reduction: " << lccReduction 
                 << ", Optionality: " << optionality 
                 << ", Betweenness: " << betweenness
                 << ", Closeness: " << closeness
                 << ", Degree: " << degree
                 << ", MISS Score: " << miss << "\n";

            if (lccReduction > bestLCCReduction || 
                (lccReduction == bestLCCReduction && optionality < bestOptionality) || 
                (lccReduction == bestLCCReduction && optionality == bestOptionality && miss > bestMISS)) {
                bestNode = candidate;
                bestLCCReduction = lccReduction;
                bestOptionality = optionality;
                bestMISS = miss;
            }
        }

        if (bestNode != -1) {
            graph[targetNode].insert(bestNode);
            graph[bestNode].insert(targetNode);
            interventionEdges.emplace_back(targetNode, bestNode);
            cout << "Added edge: (" << targetNode << ", " << bestNode << ")\n";
        }

        // 更新最大 LCC
        maxLCC = 0.0;
        for (int node : targets) {
            maxLCC = max(maxLCC, computeLCC(graph, node));
        }

        j++; // 更新迭代次數
    }

    return interventionEdges;
}
// I/O Integration
void readInput(const string& inputFile, unordered_map<int, unordered_set<int>>& graph, 
               vector<int>& targets, int& k, double& tau, double& omega_b, double& omega_c, double& omega_d) {
    ifstream file(inputFile);
    if (!file) {
        cerr << "Error: Unable to open input file.\n";
        exit(1);
    }

    int n, m;
    file >> n >> m;

    int t;
    file >> t;
    targets.resize(t);
    for (int i = 0; i < t; ++i) file >> targets[i];

    file >> k >> tau >> omega_b >> omega_c >> omega_d;

    for (int i = 0; i < m; ++i) {
        int u, v;
        file >> u >> v;
        graph[u].insert(v);
        graph[v].insert(u);
    }
}

void writeOutput(const string& outputFile, const vector<pair<int, int>>& interventionEdges, double maxLCC) {
    ofstream outFile(outputFile);
    if (!outFile) {
        cerr << "Error: Unable to open output file.\n";
        exit(1);
    }

    outFile << interventionEdges.size() << "\n";
    for (const auto& edge : interventionEdges) {
        outFile << edge.first << " " << edge.second << "\n";
    }
    outFile << fixed << setprecision(6) << maxLCC << "\n";
}

int main() {
    string inputFile = "../data/example/in.txt";
    string outputFile = "../data/example/output_ver5.txt";

    unordered_map<int, unordered_set<int>> graph;
    vector<int> targets;
    int k;
    double tau, omega_b, omega_c, omega_d;

    readInput(inputFile, graph, targets, k, tau, omega_b, omega_c, omega_d);

    auto interventionEdges = oisaAlgorithm(graph, targets, k, tau, omega_b, omega_c, omega_d);

    double maxLCC = 0.0;
    for (int node : targets) {
        // maxLCC = max(maxLCC, computeLCC(graph, node));
        double lcc = computeLCC(graph, node);
        maxLCC = max(maxLCC, lcc);
        cout << "Final Target Node: " << node << ", LCC = " << lcc << "\n";
    }

    writeOutput(outputFile, interventionEdges, maxLCC);

    return 0;
}
