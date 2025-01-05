#include <bits/stdc++.h>
using namespace std;

int N, M;
int T, K;
double TAU, OMEGA_B, OMEGA_C, OMEGA_D;
vector<int> targetNodes;
vector<pair<int, int>> edges;
vector<vector<int>> adjList;
vector<double> lcc;

void calculateAdjList() {
    adjList.resize(N + 1);
    for (int i = 0; i < M; i++) {
        adjList[edges[i].first].push_back(edges[i].second);
        adjList[edges[i].second].push_back(edges[i].first);
    }
}

double calculateLCCForNode(int v) {
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
    double lcc = count; // an edge is counted twice
    lcc /= (degree * (degree - 1));
    return lcc;
}

void calculateLCCForAllNodes() {
    lcc.resize(N + 1);
    for (int v = 1; v <= N; v++) {
        lcc[v] = calculateLCCForNode(v);
    }
    return;
}

void readInput(string fileName) {
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

void writeLCC() {
    string lccFileName = "../data/example/lcc.txt";
    ofstream lccFile(lccFileName);
    for (int i = 1; i <= N; i++) {
        lccFile << lcc[i] << endl;
    }
    lccFile.close();
    printf("LCC values are written to %s\n", lccFileName.c_str());
}

int main() {
    string inputFileName = "../data/example/in.txt";
    string outputFileName = "../data/example/out.txt";

    readInput(inputFileName);
    printf("N = %d, M = %d, T = %d, K = %d, TAU = %f, OMEGA_B = %f, OMEGA_C = %f, OMEGA_D = %f\n", N, M, T, K, TAU, OMEGA_B, OMEGA_C, OMEGA_D);

    calculateAdjList();
    calculateLCCForAllNodes();

    writeLCC();

    return 0;
}