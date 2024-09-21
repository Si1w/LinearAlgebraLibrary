#include "hw1.h"

static inline void check_size(size_t n, size_t m) {
    if (n == 0 || m == 0) {
        throw std::logic_error("Matrix size cannot be zero");
    }
}

static inline void check_empty(const Matrix& matrix) {
    for (const auto& row : matrix) {
        if (row.empty()) {
            throw std::logic_error("Matrix cannot be empty");
        }
    }
}

static inline void check_range(double min, double max) {
    if (min >= max) {
        throw std::logic_error("Invalid range");
    }
}

static inline void check_square(const Matrix& matrix) {
    if(matrix.size() != matrix[0].size()) {
        throw std::logic_error("Matrix is not square");
    }
}

namespace algebra{

    Matrix zeros(size_t n, size_t m) {
        check_size(n, m);
        return Matrix(n, std::vector<double>(m, 0));
    }

    Matrix ones(size_t n, size_t m) {
        check_size(n, m);
        return Matrix(n, std::vector<double>(m, 1));
    }

    Matrix random(size_t n, size_t m, double min, double max) {
        Matrix mat = zeros(n, m);
        check_size(n, m);
        check_range(min, max);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> num(min, max);
        for (auto& row : mat) {
            for (auto& elem : row) {
                elem = num(gen);
            }
        }
        return mat;
    }

    void show(const Matrix& matrix) {
        check_empty(matrix);
        std::cout << std::showpoint << std::setprecision(3);
        for (auto& v : matrix) {
            for (auto& elem : v) {
                std::cout << elem << "\t";
            }
            std::cout << std::endl;
        }
        std::cout << std::noshowpoint << std::endl;
    }

    Matrix multiply(const Matrix& matrix, double c) {
        check_empty(matrix);
        Matrix mat = matrix;
        for (auto& row : mat) {
            for (auto& elem : row) {
                elem *= c;
            }
        }
        return mat;
    }

    Matrix multiply(const Matrix& matrix1, const Matrix& matrix2) {
        if (matrix1.empty()) return matrix2;
        if (matrix2.empty()) return matrix1;

        size_t n1, m1, n2, m2;
        n1 = matrix1.size();
        n2 = matrix2.size();
        m1 = matrix1[0].size();
        m2 = matrix2[0].size();

        if (m1 != n2) {
            throw std::logic_error("Matrices cannot be multiplied");
        }
        Matrix mat = zeros(n1, m2);
        for (size_t i = 0; i < n1; ++i) {
            for (size_t j = 0; j < m2; ++j) {
                for (size_t k = 0; k < m1; ++k) {
                    mat[i][j] += matrix1[i][k] * matrix2[k][j];
                }
            }
        }
        return mat;
    }

    Matrix sum(const Matrix& matrix, double c) {
        if (matrix.empty()) return Matrix{};
        Matrix mat = matrix;
        for (auto& row : mat) {
            for (auto& elem : row) {
                elem += c;
            }
        }
        return mat;
    }

    Matrix sum(const Matrix& matrix1, const Matrix& matrix2) {
        if (matrix1.empty() && matrix2.empty()) return Matrix{};
        if (matrix1.empty() || matrix2.empty()) throw std::logic_error("One of matrix is empty");
        size_t n1, m1, n2, m2;
        n1 = matrix1.size();
        n2 = matrix2.size();
        m1 = matrix1[0].size();
        m2 = matrix2[0].size();
        if (n1 != n2 || m1 != m2) {
            throw std::logic_error("Matrices cannot be summed");
        }

        Matrix mat = matrix1;
        for (size_t i = 0; i < n1; ++i) {
            for (size_t j = 0; j < m1; ++j) {
                mat[i][j] += matrix2[i][j];
            }
        }
        return mat;
    }

    Matrix transpose(const Matrix& matrix) {
        if (matrix.empty()) return Matrix{};
        size_t n, m;
        n = matrix.size();
        m = matrix[0].size();
        Matrix mat = zeros(m, n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                mat[j][i] = matrix[i][j];
            }
        }
        return mat;
    }

    Matrix minor(const Matrix& matrix, size_t n, size_t m) {
        if (matrix.empty()) return Matrix{};
        check_square(matrix);
        size_t n1 = matrix.size();
        if (n >= n1 || m >= n1) {
            throw std::logic_error("Invalid index");
        }
        Matrix mat;
        for (size_t i = 0; i < n1; ++i) {
            if (i == n) continue;
            std::vector<double> row;
            for (size_t j = 0; j < n1; j++) {
                if (j == m) continue;
                row.push_back(matrix[i][j]);
            }
            mat.push_back(row);
        }
        return mat;
    }

    double determinant(const Matrix& matrix) {
        if (matrix.empty()) return 1;
        check_square(matrix);
        size_t n = matrix.size();
        if (n == 1) return matrix[0][0];
        if (n == 2) return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        double det = 0;
        for (size_t i = 0; i < n; ++i) {
            Matrix minor_mat = minor(matrix, 0, i);
            int sign = (i & 1) ? -1 : 1;
            det += matrix[0][i] * sign * determinant(minor_mat);
        }
        return det;
    }

    Matrix inverse(const Matrix& matrix) {
        if (matrix.empty()) return Matrix{};
        check_square(matrix);
        double det = determinant(matrix);
        if (det == 0) throw std::logic_error("Matrix is singular");
        size_t n = matrix.size();
        Matrix adjugate = zeros(n, n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                int sign = ((i + j) % 2 == 0) ? 1 : -1;
                Matrix minor_mat = minor(matrix, i, j);
                adjugate[j][i] = sign * determinant(minor_mat);
            }
        }
        return multiply(adjugate, 1 / det);
    }

    Matrix concatenate(const Matrix& matrix1, const Matrix& matrix2, int axis=0) {
        size_t n1, m1, n2, m2;
        n1 = matrix1.size();
        n2 = matrix2.size();
        if (n1 == 0) return matrix2;
        if (n2 == 0) return matrix1;
        if (matrix1.empty() && matrix2.empty()) return Matrix{};
        m1 = matrix1[0].size();
        m2 = matrix2[0].size();
        if (axis == 0) {
            if (m1 != m2) throw std::logic_error("Column number doesn't match");
            Matrix mat = matrix1;
            mat.insert(mat.end(), matrix2.begin(), matrix2.end());
            return mat;
        }else if (axis == 1) {
            if (n1 != n2) throw std::logic_error("Row number doesn't match");
            Matrix mat = matrix1;
            for (size_t i = 0; i < n1; ++i) {
                mat[i].insert(mat[i].end(), matrix2[i].begin(), matrix2[i].end());
            }
            return mat;
        }
        else {
            throw std::logic_error("axis must be 0 or 1");
        }
    }

    Matrix ero_swap(const Matrix& matrix, size_t r1, size_t r2) {
        size_t n = matrix.size();
        if (n == 0) return matrix;
        if (r1 >= n || r2 >= n) {
            throw std::logic_error("r1 or r2 is out of range");
        }
        Matrix mat = matrix;
        std::swap(mat[r1], mat[r2]);
        return mat;
    }

    Matrix ero_multiply(const Matrix& matrix, size_t r, double c) {
        size_t n = matrix.size();
        if (n == 0) return matrix;
        if (r >= n) {
            throw std::logic_error("r is out of range");
        }
        Matrix mat = matrix;
        auto& row = mat[r];
        for (auto& elem : row) {
            elem *= c;
        }
        return mat;
    }

    Matrix ero_sum(const Matrix& matrix, size_t r1, double c, size_t r2) {
        size_t n = matrix.size();
        if (n == 0) return matrix;
        if (r1 >= n || r2 >= n) {
            throw std::logic_error("r1 or r2 is out of range");
        }
        Matrix mat = matrix;
        auto& row1 = mat[r1];
        auto& row2 = mat[r2];
        for (size_t i = 0; i < row1.size(); ++i) {
            row2[i] += row1[i] * c;
        }
        return mat;
    }
    
    Matrix upper_triangular(const Matrix& matrix) {
        size_t n = matrix.size();
        if (n == 0) return matrix;
        check_square(matrix);
        if (n == 1) return matrix;
        Matrix mat = matrix;
        for (size_t i = 0; i < n; ++i) {
            if (mat[i][i] == 0) {
                bool swapped = false;
                for (size_t k = i + 1; k < n; ++k) {
                    if (mat[k][i] != 0) {
                        mat = ero_swap(mat, i, k);
                        swapped = true;
                        break;
                    }
                }
                if (!swapped) {
                    continue;
                }
            }
            double diag = mat[i][i];
            for (size_t j = i + 1; j < n; ++j) {
                double elem = mat[j][i];
                mat = ero_sum(mat, i, -elem / diag, j);
            }
        }
        return mat;
    }

}
