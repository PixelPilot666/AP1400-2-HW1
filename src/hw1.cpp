#include "hw1.h"

#include <iostream>
#include <vector>
#include <random>
#include <iomanip>
#include <stdexcept>

using Matrix = std::vector<std::vector<double>>;

void getRowsCols(const Matrix& matrix, size_t& rows, size_t& cols){
    rows = matrix.size();
    if (rows == 0){
        throw std::logic_error("The matrix is null.");
    }
    cols = matrix[0].size();
}


Matrix algebra::zeros(size_t n, size_t m){
    std::vector<std::vector<double>> matrix (n, std::vector<double>(m, 0.));
    return matrix;
}

Matrix algebra::ones(size_t n, size_t m){
    std::vector<std::vector<double>> matrix (n, std::vector<double>(m, 1.));
    return matrix;
}

Matrix algebra::random(size_t n, size_t m, double min, double max){
    if (max < min){
        throw std::logic_error("[algebra::random] The maximum must be bigger than minimum.");
    }

    std::vector<std::vector<double>> matrix (n, std::vector<double>(m));

    static std::random_device rd;
    static std::mt19937 generator(rd());
    std::uniform_real_distribution<double> distribution(min, max);
    
    for (size_t row = 0; row < n; ++row){
        for (size_t col = 0; col < m; ++col){
            matrix[row][col] = distribution(generator);
        }
    }
    return matrix;
}

void algebra::show(const Matrix& matrix){
    size_t rows = 0, cols = 0;
    getRowsCols(matrix, rows, cols);
    for (size_t row = 0; row < rows; ++row){
        for (size_t col = 0; col < cols; ++col){
            std::cout << std::fixed << std::setprecision(3) << matrix[row][col];
        }
        std::cout << std::endl;
    }
    // 重置输出格式到默认状态
    std::cout.unsetf(std::ios_base::floatfield); // 取消设置浮点格式
}

Matrix algebra::multiply(const Matrix& matrix, double c){
    size_t rows = 0, cols = 0;
    getRowsCols(matrix, rows, cols);

    Matrix new_matrix = zeros(rows, cols);
    
    for (size_t row = 0; row < rows; ++row){
        for (size_t col = 0; col < cols; ++col){
            new_matrix[row][col] = (matrix[row][col] * c);
        }
    }

    return new_matrix;
}

Matrix algebra::multiply(const Matrix& matrix1, const Matrix& matrix2){
    size_t rows_1 = 0, cols_1 = 0;
    size_t rows_2 = 0, cols_2 = 0;
    getRowsCols(matrix1, rows_1, cols_1);
    getRowsCols(matrix2, rows_2, cols_2);

    if(cols_1 != rows_2){
        throw std::logic_error("[algebra::multiply] matrix size is not match.");
    }

    Matrix new_matrix = zeros(rows_1, cols_2);

    for (size_t row_1 = 0; row_1 < rows_1; ++row_1) {
        for (size_t col_2 = 0; col_2 < cols_2; ++col_2) {  // 遍历新矩阵的每一列
            for (size_t k = 0; k < cols_1; ++k) {  // 中间变量 k 用于遍历共同维度
                new_matrix[row_1][col_2] += matrix1[row_1][k] * matrix2[k][col_2];
            }
        }
    }
    return new_matrix;
}

