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
    ADJ_MATRIX shortest_path_matrix;    // length of shortest path between nodes

    mutable unmapid betweenness_score;
    mutable unmapid closeness_score;

    Graph(vector<pii> edges, int num_nodes) {
        for (auto& [u, v] : edges) {
            adj_list.at(u).insert(v);
            adj_list.at(v).insert(u);
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
            neighbor_edge_counts.at(u) = count_neighbor_edges(u);
            lcc_list.at(u) = cal_lcc_of_node(u);
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
        adj_list.at(u).insert(v);
        adj_list.at(v).insert(u);

        update_lcc_and_neighbor_edge_counts(u);
        update_lcc_and_neighbor_edge_counts(v);

        // if both u and v are neighbors of w, then update LCC(w)
        for (auto& w : adj_list.at(u)) {
            if (w != v && adj_list.at(v).find(w) != adj_list.at(v).end()) {
                update_lcc_and_neighbor_edge_counts(w, 1);
            }
        }
    }

    void remove_edge(int u, int v) {
        adj_list.at(u).erase(v);
        adj_list.at(v).erase(u);

        update_lcc_and_neighbor_edge_counts(u);
        update_lcc_and_neighbor_edge_counts(v);

        // if both u and v are neighbors of w, then update LCC(w)
        for (auto& w : adj_list.at(u)) {
            if (w != v && adj_list.at(v).find(w) != adj_list.at(v).end()) {
                update_lcc_and_neighbor_edge_counts(w, -1);
            }
        }
    }

    int get_hop_distance(int u, int v) const {
        if (shortest_path_matrix.at(u).find(v) == shortest_path_matrix.at(u).end()) {
            return -1;
        }
        return shortest_path_matrix.at(u).at(v);
    }

    vi get_candidate_nodes(int v) const {
        vi candidate_nodes;
        for (auto& [u, _] : adj_list) {
            if (u != v && adj_list.at(v).find(u) == adj_list.at(v).end()) {
                candidate_nodes.emplace_back(u);
            }
        }
        return candidate_nodes;
    }

    double get_betweenness_score(int v) const {
        cal_betweenness_score();
        return betweenness_score.at(v);
    }

    double get_closeness_score(int v) const {
        cal_closeness_score(v);
        return closeness_score.at(v);
    }

    int get_degree(int v) const {
        return adj_list.at(v).size();
    }

    double get_lcc(int v) const {
        return lcc_list.at(v);
    }

    double get_lcc_reduction(int target_node, int candidate_node) const {
        // double lcc_before = get_lcc(target_node);
        // add_edge(target_node, candidate_node);
        // double lcc_after = get_lcc(target_node);
        // remove_edge(target_node, candidate_node);
        // return lcc_before - lcc_after;

        // 創建當前圖的副本
        Graph graph_copy = *this;
        // 在副本上模擬新增邊
        graph_copy.add_edge(target_node, candidate_node);
        // 計算新增邊後的 LCC
        double lcc_after = graph_copy.get_lcc(target_node);
        // 返回 LCC 減少的量
        return get_lcc(target_node) - lcc_after;
    }

private:
    int count_neighbor_edges(int v) {
        int count = 0;
        for (auto u : adj_list.at(v)) {
            for (auto w : adj_list.at(v)) {
                if (u != w && adj_list.at(u).find(w) != adj_list.at(u).end()) { // u and w are neighbors
                    count++;
                }
            }
        }
        count /= 2; // each edge is counted twice
        return count;
    }

    double cal_lcc_of_node(int v) {
        int degree = adj_list.at(v).size();
        if (degree < 2) return 0.0; // LCC is not defined for nodes with degree less than 2
        return (double)neighbor_edge_counts.at(v) / (degree * (degree - 1) / 2);
    }

    void update_lcc_and_neighbor_edge_counts(int v) {
        neighbor_edge_counts.at(v) = count_neighbor_edges(v);
        lcc_list.at(v) = cal_lcc_of_node(v);
    }

    void update_lcc_and_neighbor_edge_counts(int v, int delta) {
        neighbor_edge_counts.at(v) += delta;
        lcc_list.at(v) = cal_lcc_of_node(v);
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

    // void cal_betweenness_score() const {
    //     // use Brandes' algorithm to calculate betweenness score
    //     // complexity: O(V * (V + E))

    //     for (auto& [u, _] : adj_list) {
    //         betweenness_score.at(u) = 0.0;
    //     }

    //     for (auto& [s, _] : adj_list) {
    //         stack<int> S;
    //         unordered_map<int, vector<int>> prev;   // predecessors
    //         unordered_map<int, int> sigma;          // number of shortest paths 
    //         // unordered_map<int, int>& dist = shortest_path_matrix.at(s);    // distance from s
    //         const std::unordered_map<int, int>& dist = shortest_path_matrix.at(s); // const 引用
    //         unordered_map<int, double> delta;       // dependency of s on v

    //         for (auto& [v, _] : adj_list) {
    //             prev.at(v) = {};
    //             sigma.at(v) = 0;
    //             dist.at(v) = -1;
    //             delta.at(v) = 0.0;
    //         }

    //         sigma[s] = 1;
    //         dist[s] = 0;

    //         queue<int> Q;
    //         Q.push(s);

    //         // single source shortest path
    //         while (!Q.empty()) {
    //             int u = Q.front();
    //             Q.pop();
    //             S.push(u);

    //             for (auto& v : adj_list.at(u)) {
    //                 if (dist.at(v) < 0) {
    //                     Q.push(v);
    //                     dist.at(v) = dist.at(u) + 1;
    //                 }

    //                 if (dist.at(v) == dist.at(u) + 1) {
    //                     sigma.at(v) += sigma.at(u);
    //                     prev.at(v).push_back(u);
    //                 }
    //             }
    //         }

    //         // back propagation of dependencies
    //         while (!S.empty()) {
    //             int v = S.top();
    //             S.pop();

    //             for (auto& u : prev.at(v)) {
    //                 delta.at(u) += (sigma.at(u) / sigma.at(v)) * (1 + delta.at(v));

    //                 if (u != s) {
    //                     betweenness_score.at(u) += delta.at(u);
    //                 }
    //             }
    //         }
    //     }

    //     for (auto& [u, _] : adj_list) {
    //         betweenness_score.at(u) /= 2;  // divide by 2 because each shortest path is counted twice
    //     }

    //     return;
    // }

    void cal_betweenness_score() const {
        // 初始化所有節點的 betweenness score 為 0
        for (auto& [u, _] : adj_list) {
            betweenness_score[u] = 0.0;
        }

        // 對每個節點作為起點計算 betweenness score
        for (auto& [s, _] : adj_list) {
            // 創建本地副本來存儲距離、前驅節點、以及其他信息
            std::unordered_map<int, std::vector<int>> predecessors; // 前驅節點
            std::unordered_map<int, int> sigma; // 到達每個節點的最短路徑數
            std::unordered_map<int, int> dist; // 每個節點的距離
            std::unordered_map<int, double> delta; // 依賴度

            // 初始化
            for (auto& [v, _] : adj_list) {
                predecessors[v] = {};
                sigma[v] = 0;
                dist[v] = -1; // 設置初始距離為 -1（表示未訪問）
                delta[v] = 0.0;
            }

            sigma[s] = 1; // 起點到自己的路徑數為 1
            dist[s] = 0; // 起點到自己的距離為 0

            // BFS 計算最短路徑數和距離
            std::queue<int> Q;
            Q.push(s);
            std::stack<int> S;

            while (!Q.empty()) {
                int u = Q.front();
                Q.pop();
                S.push(u); // 將訪問過的節點壓入棧

                for (auto& v : adj_list.at(u)) {
                    // 發現未訪問的節點
                    if (dist[v] == -1) {
                        Q.push(v);
                        dist[v] = dist[u] + 1; // 更新距離
                    }

                    // 如果找到更短的路徑
                    if (dist[v] == dist[u] + 1) {
                        sigma[v] += sigma[u]; // 更新到達 v 的路徑數
                        predecessors[v].push_back(u); // 添加前驅節點
                    }
                }
            }

            // 反向計算依賴度
            while (!S.empty()) {
                int v = S.top();
                S.pop();

                for (auto& u : predecessors[v]) {
                    delta[u] += (static_cast<double>(sigma[u]) / sigma[v]) * (1 + delta[v]);
                }

                // 除去起點 s，本節點的貢獻加入到 betweenness score
                if (v != s) {
                    betweenness_score[v] += delta[v];
                }
            }
        }

        // 將每個節點的 betweenness score 除以 2，因為每條路徑被計算了兩次
        for (auto& [u, score] : betweenness_score) {
            betweenness_score[u] /= 2.0;
        }
    }

    void cal_closeness_score(int s) const {
        int total_distance = 0;
        for (auto& [_, dist] : shortest_path_matrix)  total_distance += dist.at(s);
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

int cal_optionality(Graph& graph, const InterventionTarget& target, int node) {
    /*
        Definition 4. Optionality.
        For a target t, the optionality of t denotes the number of nodes in the option set U_t ⊆ T such that for every u_t ∈ U_t, either
        1) the hop number from ut to t is no smaller than 3, or
        2) ut is two-hop away from t and adding an edge (t,u_t) does not increase the LCC of any common neighbor by more than τ.
    */
    int optionality = 0;
    for (const auto& neighbor : graph.adj_list.at(node)) {
        int hop_distance = graph.get_hop_distance(node, neighbor);
        if (hop_distance > 2) {
            optionality++;
        }
        else if (hop_distance == 2) {
            graph.add_edge(node, neighbor);
            optionality += check_lcc_less_than_tau(graph, target) ? 1 : 0;
            graph.remove_edge(node, neighbor);
        }
    }
    return optionality;
}

int calculate_minimum_edges_to_reduce_lcc(const Graph& graph, const InterventionTarget& target, int node) {
    // 當前節點的 LCC 和度數
    double current_lcc = graph.get_lcc(node);
    int current_degree = graph.get_degree(node);

    // 如果當前度數小於 2，無法降低 LCC，返回 0
    if (current_degree < 2) {
        return 0;
    }

    // 設定目標 LCC
    double target_lcc = target.TAU;

    // 初始化所需的最小新增邊數
    int kt = 0;

    // 使用迴圈計算所需的最小新增邊數
    while (true) {
        // 計算分母與分子根據 Lemma 1 的公式
        double numerator = current_lcc * current_degree * (current_degree - 1);
        double denominator = (current_degree + kt) * (current_degree + kt - 1);

        // 如果滿足條件，返回當前的 kt
        if (numerator / denominator <= target_lcc) {
            return kt;
        }

        // 否則，增加 kt 並繼續檢查
        kt++;
    }
}

bool check_remaining_edges(const Graph& graph, const InterventionTarget& target, int current_k, int used_edges) {
    // 剩餘的新增邊數
    int remaining_edges = target.K - used_edges;

    // 遍歷每個目標節點，檢查剩餘邊數是否足夠
    for (const int& node : target.targetNodes) {
        // 計算節點的當前度數
        int current_degree = graph.get_degree(node);

        // 計算將該節點的 LCC 降低到目標 LCC 所需的最小邊數
        int kt = calculate_minimum_edges_to_reduce_lcc(graph, target, node);

        // 確保剩餘邊數至少能滿足：
        // - kt（降低 LCC 所需邊數）
        // - ωd - current_degree（滿足度數限制所需邊數）
        // int required_edges = std::max(kt, target.OMEGA_D - current_degree);
        int required_edges = std::max(kt, static_cast<int>(target.OMEGA_D - current_degree));


        // 如果剩餘邊數不足以滿足需求，則返回 false
        if (remaining_edges < required_edges) {
            return false;
        }
    }

    // 如果所有目標節點的需求都能滿足，則返回 true
    return true;
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

pair<vector<pii>, double> OISA(Graph& graph, const InterventionTarget& target) {
    vector<pii> interventionEdges;
    double maxLCC = get_max_lcc_of_targets(graph, target);

    for (int k = 0; k < target.K; k++) {
        // Step 2: Find the target node with the highest LCC
        auto [targetNode, maxLCC] = get_target_with_max_lcc(graph, target);
        if (targetNode == -1) break;

        // // Step 3.1: Update the remaining edge count and check feasibility
        // if (!check_remaining_edges(graph, target, k, interventionEdges.size())) {
        //     cout << "!check_remaining_edges" << endl;
        //     break; 
        // }
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

    maxLCC = get_max_lcc_of_targets(graph, target);

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
    string outputFileName = "../data/example/out_ver6.txt";

    int num_nodes = 0;
    InterventionTarget target;
    vector<pii> edges;

    readInput(inputFileName, num_nodes, target, edges);
    printf("N = %d, M = %d, T = %d, K = %d\n", num_nodes, edges.size(), target.targetNodes.size(), target.K);
    printf("TAU = %f, OMEGA_B = %f, OMEGA_C = %f, OMEGA_D = %f\n", target.TAU, target.OMEGA_B, target.OMEGA_C, target.OMEGA_D);

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

    printf("Intervention edges: %lu\n", result.first.size());
    for (auto edge : result.first) {
        printf("Edge: %d - %d\n", edge.first, edge.second);
    }
    printf("Max LCC after intervention: %.2f\n", result.second);

    writeOutput(outputFileName, result.first, result.second);

    return 0;
}