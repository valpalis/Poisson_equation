#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <omp.h>
#include "parall_sor.cpp"

double ex_fun(double x, double y){
    return -2 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y);
    }      

double u_real(double x, double y){
    return sin(M_PI * x) * sin(M_PI * y);
    }

int main() {
    double eps = 0.000001, w = 1.25;

    std::vector<int> N = {50, 100, 150};
	std::vector<int> num_threads = {1, 2, 3, 4};
    std::vector<double> h(N.size(), 0.0);
    for (int i = 0; i< h.size(); ++i){
        h[i] = 1.0/N[i];
    }
    int K = 100000;


    printf("%-4s | %-4s | %-12s | %-10s | %-12s\n",
	       "N", "Потоки", "Время (мс)", "Итерации", "Погрешность");
    printf("=======================================================================\n");

    for (int i = 0; i < N.size(); ++i){
        std::vector<std::vector<double>> f(N[i]+1, std::vector<double>(N[i]+1, 0.0));
        std::vector<int> k = {};

        for (int q = 0; q < f.size(); ++q){
            for (int t = 0; t < f[0].size(); ++t){
                f[q][t] = ex_fun(q * h[i], t * h[i]);
                }
            } 

        for (int j = 0; j < num_threads.size(); ++j){
            double time = 0.0;
            std::vector<std::vector<double>> u(N[i]+1, std::vector<double>(N[i]+1, 0.0));

            std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

            PARALL_SOR(u, f, w, h[i], eps, K, num_threads[j], k);

            std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();

            auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
            time = elapsed_ms.count();
            
            double sum = 0.0;
            #pragma omp parallel num_threads(4)
            #pragma omp for reduction (+:sum)
            for (int m = 0; m < u.size(); ++m){
                for (int n = 0; n < u[0].size(); ++n) {
                    sum += std::pow(u[m][n] - u_real(m * h[i], n * h[i]), 2);                    
                }
            }

            printf("%-4d | %-6d | %-10.2f | %-8d | %-9.3f\n",
			       N[i],
			       num_threads[j],
			       time,
			       k.back(),
                    std::sqrt(sum));
                }
        printf("-----------------------------------------------------------------------\n");
    }    
    return 0;
}

