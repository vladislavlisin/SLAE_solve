//Метод Ричардсона с чебышевскими параметрами для трехчленной формулы. Работает для всех матриц
// Нужно знать минимальное и максимальное собственное значение. их ввод после ввода размерности
// самой матрицы и вектора столбца свободных коэф.

#include <iomanip>
#include <iostream>

using namespace std;

void multiply(int N, double **M, const double *V, double *R)
// умножение матрицы на вектор-столбец
//N- размерность, M- матрица, V- вектор, R- результат
{
    for(int i =0 ; i < N; i++)
    {
        R[i] = 0;
        for(int j = 0; j < N; j++)
            R[i] += M[i][j] * V[j];
    }
}

double norm(int n, const double* a)
{
    double norma = 0;
    for (int i = 0; i < n; i++)
        norma = norma + abs(a[i]);
    return norma;
}

int Richarson_method(int N, double **A, const double *b, double *x, double eps, double min, double max)
// N - размерность матрицы A, A - основа, b - вектор-столбец свободных коэф-ов, x - вектор решений, eps - точность
// min, max - минимальное и максимальное собственное значение
{
    int count = 1;//  количество итераций
    auto* x_k = new double [N](); // x^k = 0
    auto* x_k_minus_1 = new double [N](); // x^(k-1) = 0
    auto *tempX = new double [N](); // x^k - x^(k-1)
    auto* delta = new double [N](); // delta = Ax^k - b
    double temp_eps;
    double omega_k_plus_1; // w^k+1
    double omega_k; // w^k

    // x^1 = 0
    // x^2 = - (2/(max + min))*delta[i];
    // первое приближение мы выбираем сами т.е x^(k-1) = 0 (нулевой вектор)
    // alpha = min, betta = max

    double opt = 2/(max + min); // r
    double omega_1 = -(max - min)/(max + min); // значение первого w_1
    omega_k = omega_1;

    multiply(N, A, x_k_minus_1, delta);
    for(int i = 0; i < N; i++) // второе приближение
    {
        x_k[i] = x_k_minus_1[i] - opt * (delta[i] - b[i]);
    }

    // основной циклы
    do
    {
        omega_k_plus_1 = 1/(2/omega_1 - omega_k); // w^(k+1) = 1/(2*t_1 - w^k)
        multiply(N, A, x_k, delta); // delta = A * x^k
        double temp = omega_k*omega_k_plus_1;
        for(int i = 0; i < N; i++)
        {
            // основная итерационная трёхчленная формула
            // x^(k+1) = x^k + w_k * w_(k+1) * (x^k - x^(k-1)) - 2/(max + min)*(1 + w_k * w_(k+1)) * (A*x^k - b)
            x[i] = x_k[i] + temp*(x_k[i] - x_k_minus_1[i]) - opt*(1 + temp)*(delta[i] - b[i]);
        }

        // условие завершение итерационного процесса
        for(int i = 0; i < N; i++)
        {
            tempX[i] = x[i]-x_k[i];
        }
        temp_eps = norm(N, tempX);
        count++;

        // переобозначение переменных после вычисления приближения
        omega_k = omega_k_plus_1; // w_k = w_(K+1)

        for(int i = 0; i < N; i++)
        {
            x_k_minus_1[i] = x_k[i]; // x^(k-1) = x^k
            x_k[i] = x[i]; // x^k = x^(k+1)
        }
    }
    while (temp_eps>=eps);

    delete[] tempX;
    delete[] delta;
    delete[] x_k;
    delete[] x_k_minus_1;
    return count;
}

void print_matrix(double** A, double* b, int n)
{
    for (int q = 0; q < n; q++) // вывод матрицы A+b
    {
        for (int c = 0; c < n; c++)
        {
            cout << fixed << setprecision(8) << A[q][c] << " ";
        }
        cout << " " << fixed << setprecision(8) << b[q] << endl;
    }
}

void init_matrix(double** matrix, int N)
{
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            cin >> matrix[i][j];
        }
    }
}

void init_vector(double* b, int N)
{
    for(int i = 0; i < N; i++)
    {
        cin >> b[i];
    }
}

int main()
{
    double min, max; // минимальное, максимальное с.з
    double** A;
    double* b;
    double* x;

    int N;
    int iteration_count;
    double eps = 0.0000000001;

    cin >> N;

    x = new double[N]();
    A = new double*[N];
    for(int i = 0; i < N; i++)
    {
        A[i] = new double[N];
    }
    b = new double[N]; // СТОЛБИК СВОБОДНЫХ

    init_matrix(A, N);
    init_vector(b, N);
    print_matrix(A, b, N);

    cin >> min >> max;

    iteration_count = Richarson_method(N, A, b, x, eps, min, max);

    cout << "The final answer is:" << endl;
    for(int i = 0; i < N; i++ )
    {
        cout << "x[" << i+1 << "] = " << x[i] << endl;
    }

    cout << endl << "Iterations' number is " << iteration_count;
    return 0;
}