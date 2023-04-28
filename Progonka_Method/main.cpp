// МЕТОД ПРОГОНКИ

#include <iostream>
#include <iomanip>

using namespace std;

void init_matrix(double** array, int N);
void init_vector(double* vector, int N);
void solveMatrix(int n, double** a, const double *b, double *x);
void print_matrix(double ** A, double* b, int n);

int main() {
    int N; // DIMENSION
    cin >> N;
    double** matrix;
    double *x, *b; // нижняя, верхняя, главная диагонали матрицы

    b = new double[N]; // СТОЛБИК СВОБОДНЫХ
    x = new double[N]; // ГДЕ БУДУТ ЛЕЖАТЬ РЕШЕНИЯ
    matrix = new double*[N];
    for(int i = 0; i < N; i++)
    {
        matrix[i] = new double[N];
    }

    init_matrix(matrix, N);
    init_vector(b, N);
    print_matrix(matrix, b, N);

    solveMatrix(N, matrix, b, x);

    for(int i = 0; i < N; i++ )
    {
        cout << "x[" << i+1 << "] = " << x[i] << endl;
    }

    return 0;
}

void isMatrixCorrect(double** a, int n)
{
    for(int i = 1; i < n-1; i++)
    {
        if(abs(a[i][i]) < (abs(a[i][i-1]) + abs(a[i][i+1])) )
        {
            cout << "Condition is not met1";
            exit(1);
        }
    }

    if((abs(a[0][0]) < abs(a[0][1])) || (abs(a[n-1][n-1]) < abs(a[n-1][n-2])))
   {
        cout << "Condition is not met2";
        exit(1);
    }

    for(int i = 0; i < n; i++)
    {
        if(a[i][i] == 0)
        {
            cout << "Zero elements on main diag";
            exit(1);
        }
    }
}

void solveMatrix(int n, double** a, const double *b, double *x)
{
    isMatrixCorrect(a, n);
    auto* v = new double[n];
    auto* u = new double[n];

    v[0] = a[0][1]/(-a[0][0]); // прогон коэф
    u[0] = (-b[0]) / (-a[0][0]); // прогон коэф

    for (int i = 1; i < n; i++)
    {
        v[i] = a[i][i+1]/(-a[i][i] - a[i][i-1]*v[i-1]);
        u[i] = (a[i][i-1]*u[i - 1] - b[i]) / (-a[i][i] - a[i][i-1]*v[i-1]);
    }

    v[n-1] = 0;
    u[n-1] = (a[n-1][n-2]*u[n-2] - b[n-1]) / (-a[n-1][n-1] - a[n-1][n-2]*v[n-2]);

    x[n-1] = u[n-1];
    for (int i = n - 1; i >= 0; i--)
    {
        x[i-1] = u[i-1] + v[i-1]*x[i];
    }

    delete[] v;
    delete[] u;
}

void init_matrix(double** array, int N)
{
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            cin >> array[i][j];
        }
    }
}

void init_vector(double* vector, int N)
{
    for(int i = 0; i < N; i++)
    {
        cin >> vector[i];
    }
}

void print_matrix(double ** A, double* b, int n)
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

