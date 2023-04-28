#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

void arrOut(double** arr, int n);
double** inversion(double** a, int num);
void reverse_iterations(double** arr, int n, double eps);
double Determinant(double **mas, int m);
void arrIn(double** arr, int n)
{
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) cin >> arr[i][j];
}

int main() {
    int n;
    cin >> n;
    auto** arr = new double* [n];
    for (int i = 0; i < n; i++) arr[i] = new double[n];
    double eps = 0.00000001;
    arrIn(arr, n);

    reverse_iterations(arr, n, eps);

    arrOut(arr, n);
    return 0;
}

void arrOut(double** arr, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (abs(arr[i][j]) < 1e-10) cout << "0.0 ";
            else cout << setprecision(8) << arr[i][j] << " ";
        }
        cout << "\n";
    }
    cout << "\n";
}

// СОЗДАЁМ МАТРИЦУ ОБРАТНУЮ МАТРИЦЕ ARR МЕТОДОМ ГАУССА, ПРИСОЕДИНЕНИЕМ ЕДИНИЧНОЙ МАТРИЦЫ
void swap(double* a, double* b); //declare subroutine
double** inversion(double** a, int num)
{
    int* is, * js, i, j, k;
    int n = num;
    double temp, fmax;
    auto tp = new double* [num];
    for (i = 0; i < num; i++) tp[i] = new double[num];
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            tp[i][j] = a[i][j];
        }
    }
    is = new int[n];
    js = new int[n];
    for (k = 0; k < n; k++)
    {
        fmax = 0.0;
        for (i = k; i < n; i++) {
            for (j = k; j < n; j++)
            {
                temp = abs(tp[i][j]);//Find the maximum
                if (temp > fmax)
                {
                    fmax = temp;
                    is[k] = i; js[k] = j;
                }
            }
        }

        if ((fmax + 1.0) == 1.0)
        {
            delete[] is;
            delete[] js;
            exit(1);
        }
        if ((i = is[k]) != k)
            for (j = 0; j < n; j++)
                swap(&tp[k][j], &tp[i][j]);//Exchange pointer
        if ((j = js[k]) != k)
            for (i = 0; i < n; i++)
                swap(&tp[i][k], &tp[i][j]); //Exchange pointer
        tp[k][k] = 1.0 / tp[k][k];
        for (j = 0; j < n; j++)
            if (j != k)
                tp[k][j] *= tp[k][k];
        for (i = 0; i < n; i++)
            if (i != k)
                for (j = 0; j < n; j++)
                    if (j != k)
                        tp[i][j] = tp[i][j] - tp[i][k] * tp[k][j];
        for (i = 0; i < n; i++)
            if (i != k)
                tp[i][k] *= -tp[k][k];
    }
    for (k = n - 1; k >= 0; k--)
    {
        if ((j = js[k]) != k)
            for (i = 0; i < n; i++)
                swap(&tp[j][i], &tp[k][i]);
        if ((i = is[k]) != k)
            for (j = 0; j < n; j++)
                swap(&tp[j][i], &tp[j][k]);
    }
    delete[] is;
    delete[] js;
    return tp;
}
void swap(double* a, double* b)
{
    double c;
    c = *a;
    *a = *b;
    *b = c;
}

double vector_scalar(double const * vector1, double const * vector2, int n)
{
    double vector_scalar = 0;
    for(int i = 0; i < n; i++)
    {
        vector_scalar += vector1[i] * vector2[i];
    }
    return vector_scalar;
}

// бесконечная норма или макс. норма
double max_vector_norma(double* x, int n)
{
    double max = 0;
    for(int i = 0; i < n; i++)
    {
        if(abs(x[i]) > max) max = abs(x[i]);
    }
    return max;
}

double vector_norma_2(double* x, int n)
{
    double norma = 0;
    for (int i = 0; i < n; i++) norma += pow(x[i], 2);
    norma = sqrt(norma);
    return norma;
}

double vector_norma_1(double* x, int n)
{
    double norma = 0;
    for (int i = 0; i < n; i++) norma += abs(x[i]);
    return norma;
}

void reverse_iterations(double** arr, int n, double eps)
{
    // создаём копию матрицы arr - a, в неё записываем результат (a - lambda*E) и нахожим матрицу ей обратную
    auto** a = new double* [n];
    for (int i = 0; i < n; i++)
    {
        a[i] = new double[n];
        for (int j = 0; j < n; j++) a[i][j] = arr[i][j];
    }
    double approx_lambda;
    cin >> approx_lambda; // приблизительное собственное значение
    for (int i = 0; i < n; i++) a[i][i] -= approx_lambda; // a = (a - lambda*E) - матрица сдвига
    double det = Determinant(a, n);
    if(det == 0.0) { cout << "det = 0" << endl; exit(0); }
    a = inversion(a, n); // находим обратную для матрицы сдвига - (a - lambda*E)^(-1)

    int iterations = 0;
    auto* x_n = new double[n]();
    x_n[0] = 1; // пусть начальный вектор будет единичным
    auto* y_n = new double[n];

    double real_lambda; // искомое собственное значение
    double sigma; // чиселка, которую и нужно найти для с.з.
    double temp_eps; // значение остонова
    double norma;
    double temp_sigma = 0;

    // main iteration process
    do
    {
        for(int i = 0; i < n; i++) y_n[i] = 0;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                y_n[i] += a[i][j] * x_n[j]; // y_n = (a - lambda*E)^(-1) * x_n - основа метода
            }
        }
        sigma = vector_scalar(y_n, x_n, n) / vector_scalar(x_n, x_n, n);

        norma = max_vector_norma(y_n, n);
        for(int j = 0; j < n; j++)
        {
            x_n[j] = (y_n[j] / norma);
        }

        temp_eps = abs(temp_sigma - sigma);
        temp_sigma = sigma;
        iterations++;
    } while (temp_eps > eps);

    real_lambda = 1/sigma + approx_lambda;
    cout << "Iterations: " << iterations << endl;
    cout << "Eigenvalue = " << setprecision(8) << real_lambda << "\n";
}

void GetMatr(double **mas, double **p, int i, int j, int m)
{
    int ki, kj, di, dj;
    di = 0;
    for (ki = 0; ki<m - 1; ki++)
    { // проверка индекса строки
        if (ki == i) di = 1;
        dj = 0;
        for (kj = 0; kj<m - 1; kj++)
        { // проверка индекса столбца
            if (kj == j) dj = 1;
            p[ki][kj] = mas[ki + di][kj + dj];
        }
    }
}

double Determinant(double **mas, int m) {
    int i, k, n;
    double d;
    double **p;
    p = new double*[m];
    for (i = 0; i<m; i++)
        p[i] = new double[m];
    d = 0;
    k = 1; //(-1) в степени i
    n = m - 1;
    if (m<1) cout << "Определитель вычислить невозможно!";
    if (m == 1) {
        d = mas[0][0];
        return d;
    }
    if (m == 2) {
        d = mas[0][0] * mas[1][1] - (mas[1][0] * mas[0][1]);
        return d;
    }
    if (m>2) {
        for (i = 0; i<m; i++) {
            GetMatr(mas, p, i, 0, m);
            d = d + k * mas[i][0] * Determinant(p, n);
            k = -k;
        }
    }
    return(d);
}