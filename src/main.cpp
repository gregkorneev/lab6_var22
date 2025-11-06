// Лабораторная №6. Метод Жордана–Гаусса (полного исключения)
// -----------------------------------------------------------
// Вариант 22:
//   x + 3y - 2z = 6
//   2x - y + z  = 3
//   4x + 2y - 3z = 11
//
// Печатаем:
// 1) систему уравнений;
// 2) единичную матрицу E и столбец правых частей (решений) рядом: | E | x |;
// 3) значения x, y, z.
//
// По просьбе: НЕ печатаем полную расширенную матрицу по шагам.

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <string>

using std::cout;
using std::endl;
using std::vector;

const double EPS = 1e-9; // порог «практического нуля» для double

// ---------- Печать системы вида "x + 3y - 2z = 6" ----------
static void printSystem(const vector<vector<double>>& A, const vector<double>& b) {
    const int n = (int)A.size();
    auto varName = [&](int j) -> std::string {
        static const char* xyz[] = {"x","y","z"};
        return (n == 3 && j < 3) ? xyz[j] : "x" + std::to_string(j + 1);
    };

    cout << "Система уравнений (вариант 22):\n";
    cout << std::fixed << std::setprecision(0);

    for (int i = 0; i < n; ++i) {
        bool first = true;
        for (int j = 0; j < n; ++j) {
            double c = A[i][j];
            if (std::fabs(c) < EPS) continue;

            if (!first)  cout << (c >= 0 ? " + " : " - ");
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

// ---------- Печать блока "| E | x |" (левая часть — единичная матрица, справа — решения) ----------
static void printIdentityWithRight(const vector<vector<double>>& M) {
    // Ожидаем, что на вход подана приведённая расширенная матрица [E | x]
    // из функции Жордана–Гаусса. Печатаем только финальный вид: | E | x |
    const int n = (int)M.size();
    if (n == 0) return;

    cout << "Единичная матрица с правыми частями (финальный вид):\n";
    cout << std::fixed << std::setprecision(0);

    for (int i = 0; i < n; ++i) {
        cout << "| ";
        // Левая часть — E (округляем к 0/1 для красоты)
        for (int j = 0; j < n; ++j) {
            double v = M[i][j];
            if (std::fabs(v) < 5e-10) v = 0.0;
            else if (std::fabs(v - 1.0) < 5e-10) v = 1.0;
            cout << std::setw(2) << v << " ";
        }
        cout << "| ";

        // Правая часть — решения (вектор x), печатаем в том же «столбцовом» формате
        double rhs = M[i][n];
        cout << std::setprecision(6) << std::setw(10) << rhs << " ";
        cout << "|\n";
        cout << std::setprecision(0); // вернуть формат для следующих строк E
    }
    cout << endl;
}

// ---------- Метод Жордана–Гаусса с частичным выбором главного элемента ----------
vector<double> gaussJordan(const vector<vector<double>>& A,
                           const vector<double>& b,
                           vector<vector<double>>& Mout) // вернём финальную [E | x] для печати
{
    const int n = (int)A.size();
    vector<vector<double>> M(n, vector<double>(n + 1, 0.0)); // расширенная [A|b]

    // Собираем [A|b]
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) M[i][j] = A[i][j];
        M[i][n] = b[i];
    }

    // Полное исключение: нормировка ведущей строки + обнуление и сверху, и снизу
    for (int col = 0, row = 0; col < n && row < n; ++col, ++row) {
        // Выбор pivot по максимуму модуля в столбце col (для устойчивости)
        int pivot = row;
        for (int i = row; i < n; ++i)
            if (std::fabs(M[i][col]) > std::fabs(M[pivot][col])) pivot = i;

        if (std::fabs(M[pivot][col]) < EPS)
            throw std::runtime_error("Столбец почти нулевой: система вырождена/многосвязна.");

        if (pivot != row) std::swap(M[pivot], M[row]);

        // Нормировка ведущей строки → ведущий элемент становится 1
        double lead = M[row][col];
        for (int j = col; j <= n; ++j) M[row][j] /= lead;

        // Обнуляем во всём столбце (и над, и под) — ключ метода Жордана–Гаусса
        for (int i = 0; i < n; ++i) if (i != row) {
            double f = M[i][col];
            if (std::fabs(f) < EPS) continue;
            for (int j = col; j <= n; ++j) M[i][j] -= f * M[row][j];
        }
    }

    // Финальный вид: [E | x]
    Mout = M;

    // Решение — последний столбец
    vector<double> x(n);
    for (int i = 0; i < n; ++i) x[i] = M[i][n];
    return x;
}

int main() {
    // ===== Вариант 22 =====
    // x + 3y - 2z = 6
    // 2x - y + z  = 3
    // 4x + 2y - 3z = 11
    vector<vector<double>> A = {
        {1,  3, -2},
        {2, -1,  1},
        {4,  2, -3}
    };
    vector<double> b = {6, 3, 11};

    // 1) Печать исходной системы (для отчёта)
    printSystem(A, b);

    try {
        // 2) Решаем; получаем финальный вид [E | x] для печати
        vector<vector<double>> Mred;
        vector<double> x = gaussJordan(A, b, Mred);

        // 3) Печатаем | E | x |
        printIdentityWithRight(Mred);

        // 4) Выводим численные значения решений
        cout << std::fixed << std::setprecision(6);
        cout << "Решение:\n";
        const char* xyz[] = {"x","y","z"};
        for (int i = 0; i < (int)x.size(); ++i) {
            std::string name = (x.size() == 3 && i < 3) ? xyz[i] : ("x" + std::to_string(i+1));
            cout << name << " = " << x[i] << "\n";
        }
        // Ожидаемо: x ≈ 2.133333 (32/15), y ≈ 1.333333 (4/3), z ≈ 0.0666667 (1/15)

    } catch (const std::exception& e) {
        std::cerr << "Ошибка: " << e.what() << endl;
        return 1;
    }

    return 0;
}