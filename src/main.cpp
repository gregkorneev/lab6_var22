#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

using std::cout;
using std::endl;
using std::vector;

const double EPS = 1e-9; // числовая "погрешность нуля"

// Решает A * x = b методом Жордана–Гаусса.
// На вход подаём квадратную матрицу A и вектор b одинаковой размерности.
// Возвращает вектор решений x. Бросает исключение, если система вырождена.
vector<double> gaussJordan(vector<vector<double>> A, vector<double> b) {
    const int n = (int)A.size();

    // ----- 1) Собираем расширенную матрицу [A | b] -----
    vector<vector<double>> M(n, vector<double>(n + 1));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) M[i][j] = A[i][j];
        M[i][n] = b[i];
    }

    // ----- 2) По каждому столбцу делаем "ведущий элемент" = 1
    // и обнуляем остальные элементы столбца (сверху и снизу) -----
    for (int col = 0, row = 0; col < n && row < n; ++col, ++row) {

        // 2.1) Частичный выбор главного элемента — ищем строку с максимальным |элементом| в текущем столбце
        int pivot = row;
        for (int i = row; i < n; ++i)
            if (std::fabs(M[i][col]) > std::fabs(M[pivot][col]))
                pivot = i;

        // Если в столбце все элементы ~0, система вырождена или имеет бесконечно много решений
        if (std::fabs(M[pivot][col]) < EPS)
            throw std::runtime_error("Нулевой столбец — система вырождена или имеет бесконечное число решений.");

        // 2.2) Меняем строки местами, чтобы pivot оказался на текущей "рабочей" строке
        if (pivot != row) std::swap(M[pivot], M[row]);

        // 2.3) Нормируем строку так, чтобы ведущий элемент стал 1
        double lead = M[row][col];
        for (int j = col; j <= n; ++j) M[row][j] /= lead;

        // 2.4) Обнуляем все остальные элементы в этом столбце (и выше, и ниже)
        for (int i = 0; i < n; ++i) if (i != row) {
            double factor = M[i][col];
            if (std::fabs(factor) < EPS) continue; // уже почти ноль — пропускаем
            for (int j = col; j <= n; ++j) M[i][j] -= factor * M[row][j];
        }
    }

    // ----- 3) Читаем решения из последнего столбца -----
    vector<double> x(n);
    for (int i = 0; i < n; ++i) x[i] = M[i][n];
    return x;
}

int main() {
    vector<vector<double>> A = {
        {1, 3,  -2},
        {2,  -1, 1},
        {4,  2,  -3}
    };
    vector<double> b = {6, 3, 11};

    try {
        vector<double> x = gaussJordan(A, b);

        // Красивый вывод с фиксированной точностью
        cout << std::fixed << std::setprecision(6);
        for (int i = 0; i < (int)x.size(); ++i)
            cout << "x" << (i + 1) << " = " << x[i] << endl;

        // Ожидаем: x1 = 2.333333, x2 = 1.000000, x3 = 1.333333
        // что соответствует 7/3, 1, 4/3.

    } catch (const std::exception& e) {
        std::cerr << "Ошибка: " << e.what() << endl;
        return 1;
    }

    return 0;
}
