#include <bits/stdc++.h>
#define pii pair<int, int>
using namespace std;

vector<vector<int>> calculateAdjList(vector<pair<int, int>>& edges, int N) {
    vector<vector<int>> adjList(N + 1);
    for (auto [u, v] : edges) {
        adjList[u].push_back(v);
        adjList[v].push_back(u);
    }
    return adjList;
}

double calculateLCCForNode(int v, vector<vector<int>>& adjList) {
    int count = 0;  // count edges between nodes which are neighbors of v
    int degree = adjList[v].size();
    if (degree < 2) return 0.0; // LCC is not defined for nodes with degree less than 2
    for (auto u : adjList[v]) {
        for (auto w : adjList[u]) {
            if (find(adjList[v].begin(), adjList[v].end(), w) != adjList[v].end()) { // w is neighbor of v
                count++;
            }
        }
    }
    // printf("count of edges between neighbors of %d is %d\n", v, count);
    double lcc = (double)count / (degree * (degree - 1));
    return lcc;
}

vector<double> calculateLCCForAllNodes(vector<vector<int>>& adjList, int N) {
    vector<double> lccList(N + 1);
    for (int v = 1; v <= N; v++) {
        lccList[v] = calculateLCCForNode(v, adjList);
    }
    return lccList;
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
    vector<int> nodes;
    for (int i = 1; i <= N; i++) {
        nodes.push_back(i);
    }
    vector<vector<int>> allEdgesVec = combinations(nodes, 2);
    vector<pii> allEdges;
    for (auto edge : allEdgesVec) {
        allEdges.push_back({ edge[0], edge[1] });
    }
    for (auto edge : allEdges) {
        if (find(edges.begin(), edges.end(), edge) == edges.end()) {
            candidateEdges.push_back(edge);
        }
    }
    return candidateEdges;
}

pair<vector<pii>, double> enumeration(vector<pii>& edges, vector<double>& lccList, vector<int>& targetNodes, int N, int K) {
    vector<pii> candidateEdges = getCandidateEdges(edges, N);

    // choose K edges from candidateEdges
    vector<vector<pii>> edgeCombinations = combinations(candidateEdges, K);
    double maxLCC = getMaxLCC(lccList, targetNodes);
    vector<pii> interventionEdges;
    for (auto edgeCombination : edgeCombinations) {
        puts("New combination");
        vector<pair<int, int>> newEdges = edges;
        for (auto edge : edgeCombination) {
            newEdges.push_back(edge);
            // printf("edge (%d, %d) is added\n", edge.first, edge.second);
        }
        vector<vector<int>> newAdjList = calculateAdjList(newEdges, N);
        vector<double> newLCCList = calculateLCCForAllNodes(newAdjList, newAdjList.size() - 1);
        double localMaxLCC = getMaxLCC(newLCCList, targetNodes);
        if (localMaxLCC < maxLCC) {
            maxLCC = localMaxLCC;
            interventionEdges = edgeCombination;
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
    string lccFileName = "../data/example/lcc.txt";
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
}

int main() {
    string inputFileName = "../data/example/in.txt";
    string outputFileName = "../data/example/out.txt";


    int N, M, T, K;
    double TAU, OMEGA_B, OMEGA_C, OMEGA_D;
    vector<int> targetNodes;
    vector<pair<int, int>> edges;

    readInput(inputFileName, N, M, T, targetNodes, K, TAU, OMEGA_B, OMEGA_C, OMEGA_D, edges);
    printf("N = %d, M = %d, T = %d, K = %d, TAU = %f, OMEGA_B = %f, OMEGA_C = %f, OMEGA_D = %f\n", N, M, T, K, TAU, OMEGA_B, OMEGA_C, OMEGA_D);


    vector<vector<int>> adjList = calculateAdjList(edges, N);
    vector<double> lccList = calculateLCCForAllNodes(adjList, N);

    writeLCC(lccList);
    auto [interventionEdges, maxLCC] = enumeration(edges, lccList, targetNodes, N, K);
    writeOutput(outputFileName, interventionEdges, maxLCC);

    return 0;
}