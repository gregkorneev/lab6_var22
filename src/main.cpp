#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <string>

using std::cout;
using std::endl;
using std::vector;

const double EPS = 1e-9;

// Печать системы в виде "x + 3y - 2z = 6"
static void printSystem(const vector<vector<double>>& A, const vector<double>& b) {
    const int n = (int)A.size();
    auto varName = [&](int j)->std::string{
        static const char* xyz[] = {"x","y","z"};
        return (n==3 && j<3) ? xyz[j] : "x"+std::to_string(j+1);
    };

    cout << "Система уравнений:\n";
    cout << std::fixed << std::setprecision(0);
    for (int i = 0; i < n; ++i) {
        bool first = true;
        for (int j = 0; j < n; ++j) {
            double c = A[i][j];
            if (std::fabs(c) < EPS) continue;

            if (!first) cout << (c >= 0 ? " + " : " - ");
            else if (c < 0) cout << "-";

            double ac = std::fabs(c);
            if (std::fabs(ac - 1.0) >= EPS) cout << ac;
            cout << varName(j);
            first = false;
        }
        cout << " = " << b[i] << "\n";
    }
    cout << endl;
}

// Жордан–Гаусс с частичным выбором главного элемента
vector<double> gaussJordan(vector<vector<double>> A, vector<double> b) {
    const int n = (int)A.size();
    vector<vector<double>> M(n, vector<double>(n+1));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) M[i][j] = A[i][j];
        M[i][n] = b[i];
    }
    for (int col = 0, row = 0; col < n && row < n; ++col, ++row) {
        int pivot = row;
        for (int i = row; i < n; ++i)
            if (std::fabs(M[i][col]) > std::fabs(M[pivot][col])) pivot = i;
        if (std::fabs(M[pivot][col]) < EPS)
            throw std::runtime_error("Нулевой столбец — система вырождена/многосвязна.");
        if (pivot != row) std::swap(M[pivot], M[row]);

        double lead = M[row][col];
        for (int j = col; j <= n; ++j) M[row][j] /= lead;

        for (int i = 0; i < n; ++i) if (i != row) {
            double f = M[i][col];
            if (std::fabs(f) < EPS) continue;
            for (int j = col; j <= n; ++j) M[i][j] -= f * M[row][j];
        }
    }
    vector<double> x(n);
    for (int i = 0; i < n; ++i) x[i] = M[i][n];
    return x;
}

int main() {
    vector<vector<double>> A = {
        {1,  3, -2},
        {2, -1,  1},
        {4,  2, -3}
    };
    vector<double> b = {6, 3, 11};

    printSystem(A, b);

    try {
        vector<double> x = gaussJordan(A, b);

        cout << std::fixed << std::setprecision(6);
        cout << "Решение:\n";
        const char* xyz[] = {"x","y","z"};
        for (int i = 0; i < (int)x.size(); ++i) {
            std::string name = (x.size()==3 && i<3) ? xyz[i] : ("x"+std::to_string(i+1));
            cout << name << " = " << x[i] << "\n";
        }
        // Ожидаемо: x≈2.133333 (32/15), y≈1.333333 (4/3), z≈0.0666667 (1/15)
    } catch (const std::exception& e) {
        std::cerr << "Ошибка: " << e.what() << endl;
        return 1;
    }
    return 0;
}
