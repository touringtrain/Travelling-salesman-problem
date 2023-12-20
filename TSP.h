#include <iostream>
#include <string>
#include <queue>
#include <algorithm>
#include <getopt.h>
#include <cmath>
#include <iomanip>
#include <cstdint>

#ifndef TSP_TSP_H
#define TSP_TSP_H

using namespace std;

void get_options(int argc, char** argv, string& mode_string) {
    int option_index = 0;
    int option;

    opterr = false;

    struct option longOpts[] = {{ "mode", required_argument, nullptr, 'm' },
                                { "help", no_argument, nullptr, 'h' },
                                { nullptr, 0, nullptr, '\0' }};

    while ((option = getopt_long(argc, argv, "m:h", longOpts, &option_index)) != -1) {
        switch (option) {
            case 'm':
                mode_string = optarg;
                break;
            case 'h':
                std::cout << "Zookeeper.\n";
                exit(0);
            default:
                std::cout << "Error: Invalid command line option." << endl;
                exit(1);
        }
    }
}

class Opttsp{
public:
    Opttsp(vector<int>& p, vector<long long>& x, vector<long long>& y, double m){
        path = p;
        vertex_x_vector = x;
        vertex_y_vector = y;
        min_cost = m;
        current_cost = 0;
        optimal_path = path;
        point_num = (int)x.size();
    }

    [[nodiscard]] double get_min() const{
        return min_cost;
    }

    vector<int> get_route(){
        return optimal_path;
    }

    void genPerms(size_t permLength) {

        if (permLength == path.size()) {
            double temp_distance = sqrt((vertex_x_vector[path[permLength - 1]] - vertex_x_vector[0]) *
                                        (vertex_x_vector[path[permLength - 1]] - vertex_x_vector[0]) +
                                        (vertex_y_vector[path[permLength - 1]] - vertex_y_vector[0]) *
                                        (vertex_y_vector[path[permLength - 1]] - vertex_y_vector[0]));
            current_cost += temp_distance;
            if (current_cost < min_cost){
                optimal_path = path;
                min_cost = current_cost;
            }
            current_cost -= temp_distance;
            return;
        }

        if (!promising(permLength)) {
            return;
        }

        for (size_t i = permLength; i < path.size(); ++i) {
            swap(path[permLength], path[i]);
            double temp_distance = sqrt((vertex_x_vector[path[permLength - 1]] - vertex_x_vector[path[permLength]]) *
                                        (vertex_x_vector[path[permLength - 1]] - vertex_x_vector[path[permLength]]) +
                                        (vertex_y_vector[path[permLength - 1]] - vertex_y_vector[path[permLength]]) *
                                        (vertex_y_vector[path[permLength - 1]] - vertex_y_vector[path[permLength]]));
            current_cost += temp_distance;
            genPerms(permLength + 1);
            current_cost -= temp_distance;
            swap(path[permLength], path[i]);
        }
    }

    bool promising(size_t permLength){
        if (path.size() - permLength <= 3) return true;
        int64_t min_distance1 = numeric_limits<int64_t>::max(), min_distance2 = numeric_limits<int64_t>::max();
        for (size_t i = permLength; i < point_num; i++){
            int64_t distance1 = (vertex_x_vector[path[i]] - vertex_x_vector[0]) * (vertex_x_vector[path[i]] - vertex_x_vector[0]) +
                                (vertex_y_vector[path[i]] - vertex_y_vector[0]) * (vertex_y_vector[path[i]] - vertex_y_vector[0]);
            int64_t distance2 = (vertex_x_vector[path[i]] - vertex_x_vector[path[permLength - 1]]) *
                                (vertex_x_vector[path[i]] - vertex_x_vector[path[permLength - 1]]) +
                                (vertex_y_vector[path[i]] - vertex_y_vector[path[permLength - 1]]) *
                                (vertex_y_vector[path[i]] - vertex_y_vector[path[permLength - 1]]);
            if (distance1 < min_distance1) min_distance1 = distance1;
            if (distance2 < min_distance2) min_distance2 = distance2;
        }

        double estimate = current_cost + sqrt(min_distance1) + sqrt(min_distance2) + find_MST(permLength);
        if (estimate > min_cost) return false;
        return true;
    }

    double find_MST(size_t permLength){
        vector<pair<int, int64_t>> join_sequence(point_num - permLength, {0, numeric_limits<int64_t>::max()});
        for (size_t i = 0; i < point_num - permLength; i++){
            join_sequence[i].first = path[i + permLength];
        }
        vector<int64_t> total_distance;
        join_sequence[0].second = 0;
        for (size_t i = 1; i < point_num - permLength; i++){
            join_sequence[i].second = (vertex_x_vector[join_sequence[i].first] - vertex_x_vector[join_sequence[0].first]) *
                                      (vertex_x_vector[join_sequence[i].first] - vertex_x_vector[join_sequence[0].first]) +
                                      (vertex_y_vector[join_sequence[i].first] - vertex_y_vector[join_sequence[0].first]) *
                                      (vertex_y_vector[join_sequence[i].first] - vertex_y_vector[join_sequence[0].first]);

        }
        for (size_t i = 1; i < point_num - permLength; i++){
            int64_t min_distance2 = numeric_limits<int64_t>::max();
            size_t temp_index = 0;
            for (size_t j = i; j < point_num - permLength; j++){
                if (join_sequence[j].second < min_distance2){
                    min_distance2 = join_sequence[j].second;
                    temp_index = j;
                }
            }
            swap(join_sequence[temp_index], join_sequence[i]);
            total_distance.push_back(min_distance2);
            for (size_t k = i + 1; k < point_num - permLength; k++){
                int64_t new_distance2 = (vertex_x_vector[join_sequence[k].first] - vertex_x_vector[join_sequence[i].first]) *
                                        (vertex_x_vector[join_sequence[k].first] - vertex_x_vector[join_sequence[i].first]) +
                                        (vertex_y_vector[join_sequence[k].first] - vertex_y_vector[join_sequence[i].first]) *
                                        (vertex_y_vector[join_sequence[k].first] - vertex_y_vector[join_sequence[i].first]);
                if (new_distance2 < join_sequence[k].second){
                    join_sequence[k].second = new_distance2;
                }

            }
        }
        double total_distance_root = 0;
        for (auto i : total_distance){
            total_distance_root += sqrt(i);
        }
        return total_distance_root;
    }

private:
    double min_cost;
    double current_cost;
    size_t point_num;
    vector<int> path;
    vector<int> optimal_path;
    vector<long long> vertex_x_vector;
    vector<long long> vertex_y_vector;
};

double euclideanDistance(long long x1, long long y1, long long x2, long long y2) {
    double xDiff = (double)(x1) - (double)(x2);
    double yDiff = (double)(y1) - (double)(y2);
    double sumOfSquares = xDiff * xDiff + yDiff * yDiff;
    return sqrt(sumOfSquares);
}

double calculateTourDistance(const vector<long long>& vertex_x_vector, const vector<long long>& vertex_y_vector,
                             const vector<int>& tour) {
    double totalDistance = 0.0;
    for (size_t i = 0; i < tour.size() - 1; ++i) {
        totalDistance += sqrt((vertex_x_vector[tour[i]] - vertex_x_vector[tour[i + 1]]) *
                              (vertex_x_vector[tour[i]] - vertex_x_vector[tour[i + 1]]) +
                              (vertex_y_vector[tour[i]] - vertex_y_vector[tour[i + 1]]) *
                              (vertex_y_vector[tour[i]] - vertex_y_vector[tour[i + 1]]));
    }
    totalDistance += sqrt((vertex_x_vector[tour[0]] - vertex_x_vector[tour.back()]) *
                          (vertex_x_vector[tour[0]] - vertex_x_vector[tour.back()]) +
                          (vertex_y_vector[tour[0]] - vertex_y_vector[tour.back()]) *
                          (vertex_y_vector[tour[0]] - vertex_y_vector[tour.back()]));
    return totalDistance;
}


#endif //TSP_TSP_H
