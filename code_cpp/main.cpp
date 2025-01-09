#include <bits/stdc++.h>
#define pii pair<int, int>
#define vi vector<int>
#define unmapii unordered_map<int, int>
#define unmapid unordered_map<int, double>
#define ADJ_LIST unordered_map<int, unordered_set<int>>
#define ADJ_MATRIX unordered_map<int, unordered_map<int, int>>
using namespace std;

string USING_ALGORITHM = "OISA"; // ENUM or OISA
class InterventionTarget {
public:
    vi targetNodes;
    int K; // number of edges to intervene
    double TAU, OMEGA_B, OMEGA_C;
    int OMEGA_D;

    InterventionTarget() {}

    InterventionTarget(vi targetNodes, int K, double TAU, double OMEGA_B, double OMEGA_C, int OMEGA_D) {
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
    ADJ_MATRIX shortest_path_matrix;    // length of shortest path between nodes

    unmapid betweenness_score;
    unmapid closeness_score;

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

    int get_hop_distance(int u, int v) {
        if (shortest_path_matrix[u].find(v) == shortest_path_matrix[u].end()) {
            return -1;
        }
        return shortest_path_matrix[u][v];
    }

    vi get_candidate_nodes(int v) {
        vi candidate_nodes;
        for (auto& [u, _] : adj_list) {
            if (u != v && adj_list[v].find(u) == adj_list[v].end()) {
                candidate_nodes.emplace_back(u);
            }
        }
        return candidate_nodes;
    }

    double get_betweenness_score(int v) {
        cal_betweenness_score();
        return betweenness_score[v];
    }

    double get_closeness_score(int v) {
        cal_closeness_score(v);
        return closeness_score[v];
    }

    int get_degree(int v) {
        return adj_list[v].size();
    }

    double get_lcc(int v) {
        return lcc_list[v];
    }

    double get_lcc_reduction(int target_node, int candidate_node) {
        double lcc_before = get_lcc(target_node);
        add_edge(target_node, candidate_node);
        double lcc_after = get_lcc(target_node);
        remove_edge(target_node, candidate_node);
        return lcc_before - lcc_after;
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

    void cal_betweenness_score() {
        // use Brandes' algorithm to calculate betweenness score
        // complexity: O(V * (V + E))

        for (auto& [u, _] : adj_list) {
            betweenness_score[u] = 0.0;
        }

        for (auto& [s, _] : adj_list) {
            stack<int> S;
            unordered_map<int, vector<int>> prev;   // predecessors
            unordered_map<int, int> sigma;          // number of shortest paths 
            unordered_map<int, int>& dist = shortest_path_matrix[s];    // distance from s
            unordered_map<int, double> delta;       // dependency of s on v

            for (auto& [v, _] : adj_list) {
                prev[v] = {};
                sigma[v] = 0;
                dist[v] = -1;
                delta[v] = 0.0;
            }

            sigma[s] = 1;
            dist[s] = 0;

            queue<int> Q;
            Q.push(s);

            // single source shortest path
            while (!Q.empty()) {
                int u = Q.front();
                Q.pop();
                S.push(u);

                for (auto& v : adj_list[u]) {
                    if (dist[v] < 0) {
                        Q.push(v);
                        dist[v] = dist[u] + 1;
                    }

                    if (dist[v] == dist[u] + 1) {
                        sigma[v] += sigma[u];
                        prev[v].push_back(u);
                    }
                }
            }

            // back propagation of dependencies
            while (!S.empty()) {
                int v = S.top();
                S.pop();

                for (auto& u : prev[v]) {
                    delta[u] += (sigma[u] / sigma[v]) * (1 + delta[v]);

                    if (u != s) {
                        betweenness_score[u] += delta[u];
                    }
                }
            }
        }

        for (auto& [u, _] : adj_list) {
            betweenness_score[u] /= 2;  // divide by 2 because each shortest path is counted twice
        }

        return;
    }

    void cal_closeness_score(int s) {
        int total_distance = 0;
        for (auto& [_, dist] : shortest_path_matrix)  total_distance += dist[s];
        if (total_distance == 0) {
            closeness_score[s] = 0.0;
            return;
        }
        closeness_score[s] = (double)(adj_list.size() - 1) / total_distance;
        return;
    }

};

class Score {
public:
    double betweenness;
    double closeness;
    double degree;

    Score(Graph graph, InterventionTarget target, int node) {
        betweenness = graph.get_betweenness_score(node);
        closeness = graph.get_closeness_score(node);
        degree = graph.get_degree(node);
        set_weights(target.OMEGA_B, target.OMEGA_C);
    }

    double get_score() {
        return w_b * betweenness + w_c * closeness + w_d * degree;
    }

private:
    double w_b, w_c, w_d;

    void set_weights(double omega_b, double omega_c) {
        w_b = (omega_b - betweenness) / (omega_b + omega_c);
        w_c = (omega_c - closeness) / (omega_b + omega_c);
        w_d = 1 - w_b - w_c;
    }
};

struct Node {
    int node = -1;
    double lcc_reduction = 0.0;
    int optionality = INT_MAX;
    double miss_score = -1.0;

    Node() {}

    Node(int node, double lcc_reduction, int optionality, double miss_score) {
        this->node = node;
        this->lcc_reduction = lcc_reduction;
        this->optionality = optionality;
        this->miss_score = miss_score;
    }

    // the smaller the better
    bool operator<(const Node& other) const {
        if (lcc_reduction != other.lcc_reduction) {
            return lcc_reduction > other.lcc_reduction;
        }
        if (optionality != other.optionality) {
            return optionality < other.optionality;
        }
        return miss_score > other.miss_score;
    }
};

bool check_lcc_less_than_tau(const Graph& graph, const InterventionTarget& target) {
    for (const auto& [node, lcc] : graph.lcc_list) {
        if (lcc > target.TAU) {
            return false;
        }
    }
    return true;
}

double get_max_lcc_of_targets(const Graph& graph, const InterventionTarget& target) {
    double max_lcc = 0.0; // Initialize maximum LCC value to 0
    for (const auto& node : target.targetNodes) {
        max_lcc = max(max_lcc, graph.lcc_list.at(node));
    }
    return max_lcc;
}

pair<int, double> get_target_with_max_lcc(const Graph& graph, const InterventionTarget& target) {
    int target_node = -1;
    double max_lcc = 0.0;
    for (const auto& node : target.targetNodes) {
        if (graph.lcc_list.at(node) > max_lcc) {
            max_lcc = graph.lcc_list.at(node);
            target_node = node;
        }
    }
    return { target_node, max_lcc };
}

int get_max_degree_of_lcc_top_2k_nodes(Graph& graph, const InterventionTarget& target) {
    /*
        d_2k is the maximum degree among all top-2k nodes with the largest LCCs in T.
    */

    vector<pair<double, int>> lcc_nodes;
    for (const auto& [node, lcc] : graph.lcc_list) {
        lcc_nodes.emplace_back(lcc, node);
    }
    sort(lcc_nodes.rbegin(), lcc_nodes.rend());

    int d_2k = 0;
    for (int i = 0; i < 2 * target.K; i++) {
        d_2k = max(d_2k, graph.get_degree(lcc_nodes[i].second));
    }

    return d_2k;
}

int cal_k_t(Graph& graph, const InterventionTarget& target, int node, double l) {
    /*
        Lemma 1. Given a nodet of degree dG(t), to reduce LCC_G(t) to any targeted LCC l,
        the minimum number of intervention edges k_t is the smallest number satisfying:
        LCC_G(t) * dG(t) * (dG(t) - 1) <= l * (dG(t) + k_t) * (dG(t) + k_t - 1)
    */
    // printf("cal_k_t: node = %d\n", node);
    int degree = graph.get_degree(node);
    double lcc = get_max_lcc_of_targets(graph, target);
    int k_t = 1;
    double left = lcc * degree * (degree - 1);
    // printf(" l = %f, lcc = %f, degree = %d\n", l, lcc, degree);
    double right = l * (degree + k_t) * (degree + k_t - 1);
    assert(right > 0);

    while (left > l * (degree + k_t) * (degree + k_t - 1)) {
        // printf("k_t = %d\n", k_t);
        // printf("left = %f, right = %f\n", left, l * (degree + k_t) * (degree + k_t - 1));
        k_t++;
    }

    return k_t;
}

int cal_optionality(Graph& graph, const InterventionTarget& target, int node) {
    /*
        Definition 4. Optionality.
        For a target t, the optionality of t denotes the number of nodes in the option set U_t ⊆ T such that for every u_t ∈ U_t, either
        1) the hop number from ut to t is no smaller than 3, or
        2) ut is two-hop away from t and adding an edge (t,u_t) does not increase the LCC of any common neighbor by more than τ.
    */
    int optionality = 0;
    for (const auto& [v, _] : graph.adj_list) {
        int hop_distance = graph.get_hop_distance(node, v);
        if (hop_distance > 2) {
            optionality++;
        }
        else if (hop_distance == 2) {
            graph.add_edge(node, v);
            optionality += check_lcc_less_than_tau(graph, target) ? 1 : 0;
            graph.remove_edge(node, v);
        }
    }
    return optionality;
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
    printf("Using ENUM algorithm\n");
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

pair<vector<pii>, double> OISA(Graph& graph, const InterventionTarget& target) {
    printf("Using OISA algorithm\n");
    vector<pii> globalInterventionEdges;
    double maxLCC = get_max_lcc_of_targets(graph, target);

    int j = 1;
    int d_2k = get_max_degree_of_lcc_top_2k_nodes(graph, target);
    double l_j = 2.0 * j / (d_2k * (d_2k - 1));
    for (; l_j < maxLCC; j++, d_2k = get_max_degree_of_lcc_top_2k_nodes(graph, target), l_j = 2.0 * j / (d_2k * (d_2k - 1))) {
        printf("start iteration %d\n", j);
        // double k_G = 0; // sum of max{k_t, omega_d - d_G(t)} for all t in T
        // for (const auto& node : target.targetNodes) {
        //     double k_t = cal_k_t(graph, target, node, l_j);
        //     k_G += 0.5 * max(k_t, (double)target.OMEGA_D - graph.get_degree(node));
        //     printf("k_t = %d, omega_d = %d, d_G(t) = %d\n", k_t, target.OMEGA_D, graph.get_degree(node));
        // }
        // printf("k_G = %d, target.K = %d\n", k_G, target.K);
        // if (k_G < target.K) {
        vector<pii> interventionEdges;
        for (int k = 0; k < target.K; k++) {
            // Step 2: Find the target node with the highest LCC
            auto [targetNode, maxLCC] = get_target_with_max_lcc(graph, target);
            if (targetNode == -1) break;

            // Step 3: Find the best candidate node to connect to the target node
            Node bestNode;
            Node candidateNode;
            vi candidate_nodes = graph.get_candidate_nodes(targetNode);
            for (auto& candidate : candidate_nodes) {
                int optionality = cal_optionality(graph, target, candidate);
                double miss_score = Score(graph, target, candidate).get_score();
                double lcc_reduction = graph.get_lcc_reduction(targetNode, candidate);

                candidateNode = Node(candidate, lcc_reduction, optionality, miss_score);

                if (candidateNode < bestNode) {
                    bestNode = candidateNode;
                }
            }

            if (bestNode.node != -1) {
                graph.add_edge(targetNode, bestNode.node);
                interventionEdges.emplace_back(targetNode, bestNode.node);
            }
        }

        double localMaxLCC = get_max_lcc_of_targets(graph, target);
        printf("localMaxLCC = %f, maxLCC = %f\n", localMaxLCC, maxLCC);
        if (localMaxLCC < maxLCC) {
            maxLCC = localMaxLCC;
            globalInterventionEdges = interventionEdges;
        }
        // }
    }

    maxLCC = get_max_lcc_of_targets(graph, target);

    return { globalInterventionEdges, maxLCC };
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

    string inputFileName = "../data/example/in3.txt";
    string outputFileName = "../data/example/out3.txt";
    int num_nodes = 0;
    InterventionTarget target;
    vector<pii> edges;

    readInput(inputFileName, num_nodes, target, edges);
    printf("N = %d, M = %d, T = %d, K = %d\n", num_nodes, edges.size(), target.targetNodes.size(), target.K);
    printf("TAU = %f, OMEGA_B = %f, OMEGA_C = %f, OMEGA_D = %d\n", target.TAU, target.OMEGA_B, target.OMEGA_C, target.OMEGA_D);

    Graph graph(edges, num_nodes);

    // time start
    auto start = chrono::high_resolution_clock::now();
    pair<vector<pii>, double> result;
    if (USING_ALGORITHM == "ENUM") {
        result = enumeration(graph, target);
    }
    else {
        result = OISA(graph, target);
    }
    // time end
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;
    printf("Elapsed time: %f seconds\n", elapsed.count());


    writeOutput(outputFileName, result.first, result.second);

    return 0;
}