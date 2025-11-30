#include <iostream>
#include <cmath>
#include <vector>
#include <omp.h>
#include "parall_sor.h"

void PARALL_SOR(std::vector<std::vector<double>>& u, 
        const std::vector<std::vector<double>>& f,
        double w, double h, double eps, double K, int num_threads,
        std::vector<int>& v) {
            int m = u.size();
            int n = u[0].size();

            for (int k = 1; k <= K; ++k) {
                double sum = 0.0;
                
                #pragma omp parallel num_threads(num_threads)
                #pragma omp for reduction (+:sum)
                for (int i = 1; i < m - 1; ++i) {
                    for (int j = 1; j < n - 1; ++j) {
                        if (std::abs(i-j) % 2 == 0){
                            double u_k = u[i][j];
                            double u_k1 = 0.25 * 
                            (u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] - h*h*f[i][j]);                    
                            u[i][j] = (1-w) * u_k + w * u_k1;
                            sum += std::pow(u[i][j] - u_k, 2);
                        }
                    }
                
                }

                #pragma omp parallel num_threads(num_threads)
                #pragma omp for reduction (+:sum)
                for (int i = 1; i < m - 1; ++i) {
                    for (int j = 1; j < n - 1; ++j) {
                        if (std::abs(i-j) % 2 == 1){
                            double u_k = u[i][j];
                            double u_k1 = 0.25 * 
                            (u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] - h*h*f[i][j]);                    
                            u[i][j] = (1-w) * u_k + w * u_k1;
                            sum += std::pow(u[i][j] - u_k, 2);
                        }
                    }
                }

                if (sum < eps * eps){
                    v.push_back(k);
                    break;
                }
            }
        }
