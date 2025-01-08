#include <bits/stdc++.h>
#define pii pair<int, int>
#define vi vector<int>
#define unmapii unordered_map<int, int>
#define unmapid unordered_map<int, double>
#define ADJ_LIST unordered_map<int, unordered_set<int>>
using namespace std;

string USING_ALGORITHM = "ENUM";
class InterventionTarget {
public:
    vi targetNodes;
    int K; // number of edges to intervene
    double TAU, OMEGA_B, OMEGA_C, OMEGA_D;

    InterventionTarget() {}

    InterventionTarget(vi targetNodes, int K, double TAU, double OMEGA_B, double OMEGA_C, double OMEGA_D) {
        this->targetNodes = targetNodes;
        this->K = K;
        this->TAU = TAU;
        this->OMEGA_B = OMEGA_B;
        this->OMEGA_C = OMEGA_C;
        this->OMEGA_D = OMEGA_D;
    }

    ~InterventionTarget() {
        targetNodes.clear();
    }
};
class Graph {
public:
    ADJ_LIST adj_list;
    unmapii neighbor_edge_counts; // number of edges between neighbors for each node
    unmapid lcc_list; // local clustering coefficient for each node
    vector<pii> candidate_edges;

    Graph(vector<pii> edges, int num_nodes) {
        for (auto& [u, v] : edges) {
            adj_list[u].insert(v);
            adj_list[v].insert(u);
        }
        if (adj_list.size() != num_nodes) {
            puts("Warning: some nodes are not connected to any other nodes");
            for (int i = 1; i <= num_nodes; i++) {
                if (adj_list.find(i) == adj_list.end()) {
                    adj_list[i] = {};
                }
            }
        }

        for (auto& [u, neighbors] : adj_list) {
            neighbor_edge_counts[u] = count_neighbor_edges(u);
            lcc_list[u] = cal_lcc_of_node(u);
        }

        if (USING_ALGORITHM == "ENUM") {
            cal_candidate_edges();
        }
    }

    ~Graph() {
        adj_list.clear();
        neighbor_edge_counts.clear();
        lcc_list.clear();
    }

    void add_edge(int u, int v) {
        adj_list[u].insert(v);
        adj_list[v].insert(u);

        update_lcc_and_neighbor_edge_counts(u);
        update_lcc_and_neighbor_edge_counts(v);

        // if both u and v are neighbors of w, then update LCC(w)
        for (auto& w : adj_list[u]) {
            if (w != v && adj_list[v].find(w) != adj_list[v].end()) {
                update_lcc_and_neighbor_edge_counts(w, 1);
            }
        }
    }

    void remove_edge(int u, int v) {
        adj_list[u].erase(v);
        adj_list[v].erase(u);

        update_lcc_and_neighbor_edge_counts(u);
        update_lcc_and_neighbor_edge_counts(v);

        // if both u and v are neighbors of w, then update LCC(w)
        for (auto& w : adj_list[u]) {
            if (w != v && adj_list[v].find(w) != adj_list[v].end()) {
                update_lcc_and_neighbor_edge_counts(w, -1);
            }
        }
    }


private:
    int count_neighbor_edges(int v) {
        int count = 0;
        for (auto u : adj_list[v]) {
            for (auto w : adj_list[v]) {
                if (u != w && adj_list[u].find(w) != adj_list[u].end()) { // u and w are neighbors
                    count++;
                }
            }
        }
        count /= 2; // each edge is counted twice
        return count;
    }

    double cal_lcc_of_node(int v) {
        int degree = adj_list[v].size();
        if (degree < 2) return 0.0; // LCC is not defined for nodes with degree less than 2
        return (double)neighbor_edge_counts[v] / (degree * (degree - 1) / 2);
    }

    void update_lcc_and_neighbor_edge_counts(int v) {
        neighbor_edge_counts[v] = count_neighbor_edges(v);
        lcc_list[v] = cal_lcc_of_node(v);
    }

    void update_lcc_and_neighbor_edge_counts(int v, int delta) {
        neighbor_edge_counts[v] += delta;
        lcc_list[v] = cal_lcc_of_node(v);
    }

    // get candidate edges that can be added to the graph
    void cal_candidate_edges() {
        for (auto& [u, neighbors] : adj_list) {
            for (auto& [v, _] : adj_list) {
                if (u < v && neighbors.find(v) == neighbors.end()) {
                    candidate_edges.push_back({ u, v });
                }
            }
        }
        return;
    }
};

double get_max_lcc_of_targets(const Graph& graph, const InterventionTarget& target) {
    double max_lcc = 0.0; // Initialize maximum LCC value to 0
    for (const auto& node : target.targetNodes) {
        max_lcc = max(max_lcc, graph.lcc_list.at(node));
    }
    return max_lcc;
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


pair<vector<pii>, double> enumeration(const Graph& graph, const InterventionTarget& target) {

    // choose K edges from candidateEdges
    vector<vector<pii>> edgeCombinations = combinations(graph.candidate_edges, target.K);
    double maxLCC = get_max_lcc_of_targets(graph, target);
    vector<pii> interventionEdges;
    Graph newGraph = graph;
    for (auto& edgeCombination : edgeCombinations) {

        for (auto& edge : edgeCombination) {
            newGraph.add_edge(edge.first, edge.second);
        }

        double localMaxLCC = get_max_lcc_of_targets(newGraph, target);
        if (localMaxLCC < maxLCC) {
            maxLCC = localMaxLCC;
            interventionEdges = edgeCombination;
        }

        for (auto edge : edgeCombination) {
            newGraph.remove_edge(edge.first, edge.second);
        }
    }
    return { interventionEdges, maxLCC };
}

void readInput(string fileName, int& N, InterventionTarget& target, vector<pii>& edges) {
    ifstream file(fileName);

    if (!file.is_open()) {
        cout << "File not found!" << endl;
        return;
    }

    int M, T;
    file >> N >> M >> T;

    for (int i = 0; i < T; i++) {
        int targetNode;
        file >> targetNode;
        target.targetNodes.push_back(targetNode);
    }

    file >> target.K >> target.TAU >> target.OMEGA_B >> target.OMEGA_C >> target.OMEGA_D;

    for (int i = 0; i < M; i++) {
        int u, v;
        file >> u >> v;
        if (u > v) swap(u, v);
        edges.push_back({ u, v });
    }

    file.close();
    return;
}

void writeOutput(string fileName, vector<pii> interventionEdges, double maxLCC) {
    ofstream file(fileName);
    file << interventionEdges.size() << endl;
    for (auto edge : interventionEdges) {
        file << edge.first << " " << edge.second << endl;
    }
    file << fixed << setprecision(2) << maxLCC << endl;
    file.close();
    printf("Output is written to %s\n", fileName.c_str());
}

int main() {

    string inputFileName = "../data/example/in2.txt";
    string outputFileName = "../data/example/out2.txt";

    int num_nodes = 0;
    InterventionTarget target;
    vector<pii> edges;

    readInput(inputFileName, num_nodes, target, edges);
    printf("N = %d, M = %d, T = %d, K = %d\n", num_nodes, edges.size(), target.targetNodes.size(), target.K);
    printf("TAU = %f, OMEGA_B = %f, OMEGA_C = %f, OMEGA_D = %f\n", target.TAU, target.OMEGA_B, target.OMEGA_C, target.OMEGA_D);

    Graph graph(edges, num_nodes);

    // time start
    auto start = chrono::high_resolution_clock::now();
    auto [interventionEdges, maxLCC] = enumeration(graph, target);
    // time end
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;
    printf("Elapsed time: %f seconds\n", elapsed.count());


    writeOutput(outputFileName, interventionEdges, maxLCC);

    return 0;
}