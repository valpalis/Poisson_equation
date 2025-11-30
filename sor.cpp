#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>



void SOR(std::vector<std::vector<double>>& u, 
        const std::vector<std::vector<double>>& f,
        double w, double h, double eps, double K) {
            int m = u.size();
            int n = u[0].size();

            for (int k = 1; k <= K; ++k) {
                double sum = 0.0;

                for (int i = 1; i < m - 1; ++i) {
                    for (int j = 1; j < n - 1; ++j) {
                        double u_k = u[i][j];
                        double u_k1 = 0.25 * 
                        (u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] - h*h*f[i][j]);                    
                        u[i][j] = (1-w) * u_k + w * u_k1;
                        sum += std::pow(u[i][j] - u_k, 2);
                    }
                }

                if (sum < eps * eps){
                    std::cout << "iterations: " << k << "\n" << "delta: " << std::sqrt(sum) << "\n";
                    break;
                }
            }
        }


 double ex_fun(double x, double y){
    return -2 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y);
    }      

double u_real(double x, double y){
    return sin(M_PI * x) * sin(M_PI * y);
    }

int main() {
    int N = 100, K = 100000;
    double eps = 0.000001, h = 1.0/N, w = 1.25;
    std::vector<std::vector<double>> f(N+1, std::vector<double>(N+1, 0.0));
    std::vector<std::vector<double>> u(N+1, std::vector<double>(N+1, 0.0));

    for (int i = 0; i < f.size(); ++i){
        for (int j = 0; j < f[0].size(); ++j){
            f[i][j] = ex_fun(i * h, j * h);
        }
    }   
    auto start_time = std::chrono::steady_clock::now();

    SOR(u, f, w, h, eps, K);

    auto end_time = std::chrono::steady_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    std::cout << "elapsed time: " << elapsed_ms.count() << " ms\n";

    double sum = 0.0;

    for (int i = 0; i < u.size(); ++i){
        for (int j = 0; j < u[0].size(); ++j) {
            sum += std::pow(u[i][j] - u_real(i * h, j * h), 2);
        }
    }
    std::cout << "L_2 norm: " << std::sqrt(sum) << "\n";

    return 0;
}


