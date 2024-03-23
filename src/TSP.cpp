#include "TSP.hpp"

#include <algorithm>
#include <stack>
#include <optional>

std::vector<size_t> find_path(std::vector<vertex_t> edges){
    std::vector<size_t> path {edges[0].row, edges[0].col};
    edges.erase(edges.begin());

    while(!edges.empty()){
        bool found = false;
        for(size_t i=0; i<edges.size(); i++){
            if(path.back() == edges[i].row){
                path.push_back(edges[i].col);
                auto it = edges.begin();
                std::advance(it, i);
                edges.erase(it);
                found = true;
            }
            else if(path.front() == edges[i].col){
                path.insert(path.begin(),edges[i].row);
                auto it = edges.begin();
                std::advance(it, i);
                edges.erase(it);
                found = true;
            }
        }
        if (!found){
            return path;
        }
    }
    return path;
}



std::ostream& operator<<(std::ostream& os, const CostMatrix& cm) {
    for (std::size_t r = 0; r < cm.size(); ++r) {
        for (std::size_t c = 0; c < cm.size(); ++c) {
            const auto& elem = cm[r][c];
            os << (is_inf(elem) ? "INF" : std::to_string(elem)) << " ";
        }
        os << "\n";
    }
    os << std::endl;

    return os;
}


path_t StageState::get_path() {
    reduce_cost_matrix();
    NewVertex ver=choose_new_vertex();
    matrix_[ver.coordinates.row][ver.coordinates.col]=INF;
    append_to_path(ver.coordinates);
    reduce_cost_matrix();
    ver=choose_new_vertex();
    update_cost_matrix(ver.coordinates);
    append_to_path(ver.coordinates);

    std::vector<std::size_t> path = {unsorted_path_[0].row, unsorted_path_[0].col};
    std::size_t s=unsorted_path_[0].col, i=1;
    while (path.size()<unsorted_path_.size()){
        if (s==unsorted_path_[i].row){
            s=unsorted_path_[i].col;
            path.push_back(s);
        }
        if (i==unsorted_path_.size()-1)
            i=1;
        else
            i+=1;
    }
    return path;
}

std::vector<cost_t> CostMatrix::get_min_values_in_rows() const {
    std::vector<cost_t> min_values;
    for (std::size_t i=0; i<matrix_.size(); i++){
        cost_t value = *std::min_element(matrix_[i].begin(), matrix_[i].end());
        if (value==INF) value=0;
        min_values.push_back(value);
    }
    return min_values;
}

cost_t CostMatrix::reduce_rows() {
    std::vector<cost_t> to_reduce = get_min_values_in_rows();
    for (std::size_t i=0; i<matrix_.size(); i++)
        for(std::size_t j=0; j<matrix_.size(); j++){
            if (!is_inf(matrix_[i][j]))
                matrix_[i][j]-=to_reduce[i];
        }
    return std::accumulate(to_reduce.begin(), to_reduce.end(), 0);
}

std::vector<cost_t> CostMatrix::get_min_values_in_cols() const {
    std::vector<cost_t> min_values;
    for (std::size_t i=0; i<matrix_.size(); i++){
        cost_t value = INF;
        for (std::size_t j=0; j<matrix_.size(); j++){
            if (value>matrix_[j][i]){value = matrix_[j][i];}
        }
        if (value==INF) value=0;
        min_values.push_back(value);
    }
    return min_values;
}

cost_t CostMatrix::reduce_cols() {
    std::vector<cost_t> to_reduce = get_min_values_in_cols();
    for (std::size_t i=0; i<matrix_.size(); i++)
        for(std::size_t j=0; j<matrix_.size(); j++){
            if (!is_inf(matrix_[i][j]))
                matrix_[i][j]-=to_reduce[j];
        }
    return std::accumulate(to_reduce.begin(), to_reduce.end(), 0);
}

cost_t CostMatrix::get_vertex_cost(std::size_t row, std::size_t col) const {
    cost_t min_r=INF;
    cost_t min_c=INF;
    for (std::size_t i=0;i<matrix_.size();i++) {
        if (i != col && min_r>matrix_[row][i])
            min_r=matrix_[row][i];
        if(i!=row && min_c>matrix_[i][col])
            min_c=matrix_[i][col];
    }
    if (min_r==INF) min_r=0;
    if (min_c==INF) min_c=0;
    return min_r + min_c;
}

NewVertex StageState::choose_new_vertex() {
    std::size_t r,c;
    vertex_t v(0,0);
    NewVertex ver(v,-1);
    for (std::size_t i=0; i<matrix_.size(); i++)
        for(std::size_t j=0; j<matrix_.size(); j++)
            if (matrix_[i][j]==0) {
                cost_t cost = matrix_.get_vertex_cost(i, j);
                if (cost > ver.cost) {
                    ver.coordinates.row = i;
                    ver.coordinates.col = j;
                    ver.cost=cost;
                }
            }
    return ver;
}

void StageState::update_cost_matrix(vertex_t new_vertex) {
    for (std::size_t i=0; i<matrix_.size(); i++){
        matrix_[new_vertex.row][i]=INF;
        matrix_[i][new_vertex.col]=INF;
    }
    matrix_[new_vertex.col][new_vertex.row]=INF;
    // zabronienie cyklu
    std::vector<vertex_t> edges(unsorted_path_);
    while (!edges.empty()){
        std::vector<std::size_t> path = find_path(edges);
        matrix_[path.back()][path[0]] = INF;
        edges.erase(edges.begin());
    }
}

cost_t StageState::reduce_cost_matrix() {
    cost_t reduced = matrix_.reduce_rows() + matrix_.reduce_cols();
    return reduced;
}

cost_t get_optimal_cost(const path_t& optimal_path, const cost_matrix_t& m) {
    cost_t cost = 0;

    for (std::size_t idx = 1; idx < optimal_path.size(); ++idx) {
        cost += m[optimal_path[idx - 1]][optimal_path[idx]];
    }

    cost += m[optimal_path[optimal_path.size() - 1]][optimal_path[0]];

    return cost;
}

StageState create_right_branch_matrix(cost_matrix_t m, vertex_t v, cost_t lb) {
    CostMatrix cm(m);
    cm[v.row][v.col] = INF;
    return StageState(cm, {}, lb);
}

tsp_solutions_t filter_solutions(tsp_solutions_t solutions) {
    cost_t optimal_cost = INF;
    for (const auto& s : solutions) {
        optimal_cost = (s.lower_bound < optimal_cost) ? s.lower_bound : optimal_cost;
    }

    tsp_solutions_t optimal_solutions;
    std::copy_if(solutions.begin(), solutions.end(),
                 std::back_inserter(optimal_solutions),
                 [&optimal_cost](const tsp_solution_t& s) { return s.lower_bound == optimal_cost; }
    );

    return optimal_solutions;
}

tsp_solutions_t solve_tsp(const cost_matrix_t& cm) {

    StageState left_branch(cm);

    std::stack<StageState> tree_lifo;

    std::size_t n_levels = cm.size() - 2;

    tree_lifo.push(left_branch);

    cost_t best_lb = INF;
    tsp_solutions_t solutions;

    while (!tree_lifo.empty()) {

        left_branch = tree_lifo.top();
        tree_lifo.pop();

        while (left_branch.get_level() != n_levels && left_branch.get_lower_bound() <= best_lb) {

            if (left_branch.get_level() == 0) {
                left_branch.reset_lower_bound();
            }

            cost_t new_cost = left_branch.reduce_cost_matrix();

            left_branch.update_lower_bound(new_cost);
            if (left_branch.get_lower_bound() > best_lb) {
                break;
            }

            NewVertex new_vertex = left_branch.choose_new_vertex();

            left_branch.append_to_path(new_vertex.coordinates);

            left_branch.update_cost_matrix(new_vertex.coordinates);

            cost_t new_lower_bound = left_branch.get_lower_bound() + new_vertex.cost;
            tree_lifo.push(create_right_branch_matrix(cm, new_vertex.coordinates,
                                                      new_lower_bound));
        }

        if (left_branch.get_lower_bound() <= best_lb) {
            best_lb = left_branch.get_lower_bound();
            path_t new_path = left_branch.get_path();
            solutions.push_back({get_optimal_cost(new_path, cm), new_path});
        }
    }

    return filter_solutions(solutions);
}
