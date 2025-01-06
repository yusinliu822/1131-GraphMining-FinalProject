#include <bits/stdc++.h>
#define pii pair<int, int>
#define ADJ_LIST unordered_map<int, unordered_set<int>>
using namespace std;

ADJ_LIST calculateAdjList(vector<pair<int, int>>& edges, int N) {
    ADJ_LIST adjList(N + 1);
    for (auto [u, v] : edges) {
        adjList[u].insert(v);
        adjList[v].insert(u);
    }
    return adjList;
}
void updateAdjList(ADJ_LIST& adjList, int u, int v) {
    adjList[u].insert(v);
    adjList[v].insert(u);
}

void restoreAdjList(ADJ_LIST& adjList, int u, int v) {
    adjList[u].erase(v);
    adjList[v].erase(u);
}

int countNeighborEdges(int v, ADJ_LIST& adjList) {
    int count = 0;
    auto neighbors = adjList[v];
    for (auto u : neighbors) {
        for (auto w : neighbors) {
            if (u != w && adjList[u].find(w) != adjList[u].end()) { // u and w are neighbors
                count++;
            }
        }
    }
    count /= 2; // each edge is counted twice
    return count;
}

unordered_map<int, int> calculateNeighborEdgeCounts(ADJ_LIST& adjList, int N) {
    unordered_map<int, int> neighborEdgeCounts;
    for (int v = 1; v <= N; v++) {
        neighborEdgeCounts[v] = countNeighborEdges(v, adjList);
    }
    return neighborEdgeCounts;
}

double calculateLCCForNode(int v, ADJ_LIST& adjList) {
    int degree = adjList[v].size();
    if (degree < 2) return 0.0; // LCC is not defined for nodes with degree less than 2
    int neighborEdgeCount = countNeighborEdges(v, adjList);
    return (double)neighborEdgeCount / (degree * (degree - 1) / 2);
}

vector<double> calculateLCCForAllNodes(ADJ_LIST& adjList, int N) {
    vector<double> lccList(N + 1);
    for (int v = 1; v <= N; v++) {
        lccList[v] = calculateLCCForNode(v, adjList);
    }
    return lccList;
}

void updateLCCForNode(int v, vector<double>& lccList, ADJ_LIST& adjList, unordered_map<int, int>& neighborEdgeCounts) {
    int degree = adjList[v].size();
    if (degree < 2) {
        lccList[v] = 0.0;
        return;
    }
    int neighborEdgeCount = neighborEdgeCounts[v];
    lccList[v] = (double)neighborEdgeCount / (degree * (degree - 1) / 2);
}

void updateLCC(vector<double>& lccList, ADJ_LIST& adjList, unordered_map<int, int>& neighborEdgeCounts, int u, int v) {
    lccList[u] = calculateLCCForNode(u, adjList);
    lccList[v] = calculateLCCForNode(v, adjList);

    // if both u and v are neighbors of w, then update LCC(w)
    for (auto w : adjList[u]) {
        if (w != v && adjList[v].find(w) != adjList[v].end()) {
            neighborEdgeCounts[w]++;
            updateLCCForNode(w, lccList, adjList, neighborEdgeCounts);
        }
    }
}

void restoreLCC(vector<double>& lccList, ADJ_LIST& adjList, unordered_map<int, int>& neighborEdgeCounts, int u, int v) {
    lccList[u] = calculateLCCForNode(u, adjList);
    lccList[v] = calculateLCCForNode(v, adjList);

    // if both u and v are neighbors of w, then update LCC(w)
    for (auto w : adjList[u]) {
        if (w != v && adjList[v].find(w) != adjList[v].end()) {
            neighborEdgeCounts[w]--;
            updateLCCForNode(w, lccList, adjList, neighborEdgeCounts);
        }
    }
}

double getMaxLCC(vector<double>& lccList, vector<int>& targetNodes) {
    double maxLCC = 0.0;
    for (auto targetNode : targetNodes) {
        maxLCC = max(maxLCC, lccList[targetNode]);
    }
    return maxLCC;
}


// input: an vector of n elements, k elements to choose
// output: all possible combinations of k elements from n elements
// use template to support any type of elements
template <typename T>
vector<vector<T>> combinations(vector<T> elements, int k) {
    vector<vector<T>> result;
    vector<bool> mask(elements.size());
    fill(mask.begin(), mask.begin() + k, true);

    do {
        vector<T> combination;
        for (int i = 0; i < elements.size(); i++) {
            if (mask[i]) {
                combination.push_back(elements[i]);
            }
        }
        result.push_back(combination);
    } while (prev_permutation(mask.begin(), mask.end()));

    return result;
}

vector<pii> getCandidateEdges(vector<pii>& edges, int N) {
    vector<pii> candidateEdges;
    for (int u = 1; u <= N; u++) {
        for (int v = u + 1; v <= N; v++) {
            if (find(edges.begin(), edges.end(), make_pair(u, v)) == edges.end()) {
                candidateEdges.push_back({ u, v });
            }
        }
    }
    return candidateEdges;
}

pair<vector<pii>, double> enumeration(vector<pii>& edges, ADJ_LIST& adjList, unordered_map<int, int>& neighborEdgeCounts, vector<double>& lccList, vector<int>& targetNodes, int N, int K) {
    vector<pii> candidateEdges = getCandidateEdges(edges, N);

    // choose K edges from candidateEdges
    vector<vector<pii>> edgeCombinations = combinations(candidateEdges, K);
    double maxLCC = getMaxLCC(lccList, targetNodes);
    vector<pii> interventionEdges;
    for (auto edgeCombination : edgeCombinations) {

        for (auto edge : edgeCombination) {
            updateAdjList(adjList, edge.first, edge.second);
            updateLCC(lccList, adjList, neighborEdgeCounts, edge.first, edge.second);
        }
        double localMaxLCC = getMaxLCC(lccList, targetNodes);
        if (localMaxLCC < maxLCC) {
            maxLCC = localMaxLCC;
            interventionEdges = edgeCombination;
        }

        for (auto edge : edgeCombination) {
            restoreAdjList(adjList, edge.first, edge.second);
            restoreLCC(lccList, adjList, neighborEdgeCounts, edge.first, edge.second);
        }
    }
    return { interventionEdges, maxLCC };
}

void readInput(string fileName, int& N, int& M, int& T, vector<int>& targetNodes, int& K, double& TAU, double& OMEGA_B, double& OMEGA_C, double& OMEGA_D, vector<pair<int, int>>& edges) {
    ifstream file(fileName);

    if (!file.is_open()) {
        cout << "File not found!" << endl;
        return;
    }

    file >> N >> M >> T;
    for (int i = 0; i < T; i++) {
        int targetNode;
        file >> targetNode;
        targetNodes.push_back(targetNode);
    }
    file >> K;
    file >> TAU >> OMEGA_B >> OMEGA_C >> OMEGA_D;
    for (int i = 0; i < M; i++) {
        int u, v;
        file >> u >> v;
        edges.push_back({ u, v });
    }
}

void writeLCC(vector<double>& lcc) {
    string lccFileName = "../data/test/lcc.txt";
    ofstream lccFile(lccFileName);
    for (int i = 1; i <= lcc.size(); i++) {
        lccFile << lcc[i] << endl;
    }
    lccFile.close();
    printf("LCC values are written to %s\n", lccFileName.c_str());
}

void writeOutput(string fileName, vector<pii> interventionEdges, double maxLCC) {
    ofstream file(fileName);
    file << interventionEdges.size() << endl;
    for (auto edge : interventionEdges) {
        file << edge.first << " " << edge.second << endl;
    }
    file << maxLCC << endl;
    file.close();
    printf("Output is written to %s\n", fileName.c_str());
}

int main() {
    // string inputFileName = "../data/example/in.txt";
    // string outputFileName = "../data/example/out.txt";
    string inputFileName = "../data/test/3980_small.txt";
    string outputFileName = "../data/test/3980_small_out.txt";


    int N, M, T, K;
    double TAU, OMEGA_B, OMEGA_C, OMEGA_D;
    vector<int> targetNodes;
    vector<pair<int, int>> edges;

    readInput(inputFileName, N, M, T, targetNodes, K, TAU, OMEGA_B, OMEGA_C, OMEGA_D, edges);
    printf("N = %d, M = %d, T = %d, K = %d, TAU = %f, OMEGA_B = %f, OMEGA_C = %f, OMEGA_D = %f\n", N, M, T, K, TAU, OMEGA_B, OMEGA_C, OMEGA_D);


    ADJ_LIST adjList = calculateAdjList(edges, N);
    unordered_map<int, int> neighborEdgeCounts = calculateNeighborEdgeCounts(adjList, N);
    vector<double> lccList = calculateLCCForAllNodes(adjList, N);

    writeLCC(lccList);

    // time start
    auto start = chrono::high_resolution_clock::now();
    auto [interventionEdges, maxLCC] = enumeration(edges, adjList, neighborEdgeCounts, lccList, targetNodes, N, K);
    // time end
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;
    printf("Elapsed time: %f seconds\n", elapsed.count());


    writeOutput(outputFileName, interventionEdges, maxLCC);

    return 0;
}