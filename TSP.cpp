#include <iostream>
#include <string>
#include <queue>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include "TSP.h"

using namespace std;

int main(int argc, char** argv){
    std::ios::sync_with_stdio(false);
    cout << std::setprecision(2);
    cout << std::fixed;
    string mode_string;
    get_options(argc, argv, mode_string);
    int point_num;
    cin >> point_num;
    vector<long long> vertex_x_vector(point_num);
    vector<long long> vertex_y_vector(point_num);
    bool edge = false;
    bool all_safe = true;
    bool all_wild = true;
    for (int i = 0; i < point_num; i++){
        int x, y;
        cin >> x >> y;
        if ((x <= 0 && y == 0) || (x == 0 && y <= 0)){
            edge = true;
        }
        if (x > 0 || y > 0){
            all_wild = false;
        }
        if (x < 0 && y < 0){
            all_safe = false;
        }
        vertex_x_vector[i] = x;
        vertex_y_vector[i] = y;
    }
    if (mode_string == "MST" && !all_safe && !all_wild && !edge) {
        cout << "Cannot construct MST\n";
        return 0;
    }
    if (mode_string == "FASTTSP"){
        vector<int> route(2, 0);
        vector<bool> visited(point_num, false);
        vector<double> total_increase(point_num - 1);
        visited[0] = true;
        for (int i = 1; i < point_num; i++){
            double min_increase2 = numeric_limits<double>::max();
            int index = -1;
            for (int j = 0; j < i; j++){
                double increase2 = euclideanDistance(vertex_x_vector[i], vertex_y_vector[i],
                                                     vertex_x_vector[route[j]], vertex_y_vector[route[j]]) +
                                   euclideanDistance(vertex_x_vector[i],
                                                     vertex_y_vector[i], vertex_x_vector[route[j + 1]], vertex_y_vector[route[j + 1]]) -
                                   euclideanDistance(vertex_x_vector[route[j + 1]], vertex_y_vector[route[j + 1]],
                                                     vertex_x_vector[route[j]], vertex_y_vector[route[j]]);
                if (increase2 < min_increase2){
                    index = j;
                    min_increase2 = increase2;
                }
            }
            route.insert(route.begin() + index + 1, i);
        }
        double curr_distance = calculateTourDistance(vertex_x_vector, vertex_y_vector, route);
        route.pop_back();
        cout << curr_distance << "\n";
        for (int i : route){
            cout << i << " ";
        }
        cout << "\n";
    }
    else if(mode_string == "OPTTSP"){
        vector<int> route(point_num);
        vector<int64_t> total_distance;
        for (int i = 0; i < point_num; i++){
            route[i] = i;
        }
        for (int i = 1; i < point_num; i++){
            int temp_index = -1;
            int64_t min_distance2 = numeric_limits<int64_t>::max();
            for (int j = i; j < point_num; j++){
                int64_t distance2 = (vertex_x_vector[route[j]] - vertex_x_vector[route[i - 1]]) *
                                    (vertex_x_vector[route[j]] - vertex_x_vector[route[i - 1]]) +
                                    (vertex_y_vector[route[j]] - vertex_y_vector[route[i - 1]]) *
                                    (vertex_y_vector[route[j]] - vertex_y_vector[route[i - 1]]);
                if (distance2 < min_distance2){
                    temp_index = j;
                    min_distance2 = distance2;
                }
            }
            total_distance.push_back(min_distance2);
            swap(route[temp_index], route[i]);
        }
        double curr_distance = 0;
        for (auto i : total_distance){
            curr_distance += sqrt(i);
        }
        curr_distance += sqrt((vertex_x_vector[route[0]] - vertex_x_vector[route[point_num - 1]]) *
                              (vertex_x_vector[route[0]] - vertex_x_vector[route[point_num - 1]]) +
                              (vertex_y_vector[route[0]] - vertex_y_vector[route[point_num - 1]]) *
                              (vertex_y_vector[route[0]] - vertex_y_vector[route[point_num - 1]]));
        Opttsp partc = Opttsp(route, vertex_x_vector, vertex_y_vector, curr_distance);
        partc.genPerms(1);
        cout << partc.get_min() << "\n";
        route = partc.get_route();
        for (int i : route){
            cout << i << " ";
        }
        cout << "\n";
    }
    else{
        cout << "Error: Invalid mode\n";
        exit(1);
    }

}

