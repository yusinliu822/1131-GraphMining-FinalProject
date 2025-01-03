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

void calculateLCC() {
    lcc.resize(N + 1);
    for (int v = 1; v <= N; v++) {
        int count = 0;  // count edges between nodes which are neighbors of v
        for (auto u : adjList[v]) {
            for (auto w : adjList[u]) {
                if (find(adjList[v].begin(), adjList[v].end(), w) != adjList[v].end()) { // w is neighbor of v
                    count++;
                }
            }
        }
        // printf("count of edges between neighbors of %d is %d\n", v, count);
        lcc[v] = count; // an edge is counted twice
        lcc[v] /= adjList[v].size() * (adjList[v].size() - 1);
    }
    return;
}

int main() {
    cin >> N >> M;
    cin >> T;
    for (int i = 0; i < T; i++) {
        int targetNode;
        cin >> targetNode;
        targetNodes.push_back(targetNode);
    }
    cin >> K;
    cin >> TAU >> OMEGA_B >> OMEGA_C >> OMEGA_D;
    for (int i = 0; i < M; i++) {
        int u, v;
        cin >> u >> v;
        edges.push_back({ u, v });
    }
    printf("N = %d, M = %d, T = %d, K = %d, TAU = %f, OMEGA_B = %f, OMEGA_C = %f, OMEGA_D = %f\n", N, M, T, K, TAU, OMEGA_B, OMEGA_C, OMEGA_D);
    calculateAdjList();
    calculateLCC();

    // write lcc to file ../data/lcc.txt
    ofstream lccFile;
    lccFile.open("../data/lcc.txt");
    for (int i = 1; i <= N; i++) {
        lccFile << lcc[i] << endl;
    }
    lccFile.close();
    puts("lcc is written to file ../data/lcc.txt");


    return 0;
}