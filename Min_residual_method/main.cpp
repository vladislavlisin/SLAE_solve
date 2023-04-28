// Итерационный метод минимальных невязок для решения слау
// для симметричных положительно определенных матриц

#include <iostream>
#include <cstdlib>
#include <iomanip>

using namespace std;

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

void transpos(double** A, int N)
{
    double** matrix;
    matrix = new double*[N];
    for(int i = 0; i < N; i++)
    {
        matrix[i] = new double[N];
    }

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            matrix[i][j] = A[i][j];
        }
    }

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            A[i][j] = matrix[j][i];
        }
    }

    for(int i = 0; i < N; i++)
    {
        delete[] matrix[i];
    }
    delete[] matrix;
}

// Получение матрицы без i-й строки и j-го столбца
void GetMatr(double **mas, double **p, int i, int j, int m) {
    int ki, kj, di, dj;
    di = 0;
    for (ki = 0; ki<m - 1; ki++) { // проверка индекса строки
        if (ki == i) di = 1;
        dj = 0;
        for (kj = 0; kj<m - 1; kj++) { // проверка индекса столбца
            if (kj == j) dj = 1;
            p[ki][kj] = mas[ki + di][kj + dj];
        }
    }
}

// Рекурсивное вычисление определителя
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
    if (m < 1) cout << "Det cant be calculated!";
    if (m == 1) {
        d = mas[0][0];
        return (d);
    }
    if (m == 2) {
        d = mas[0][0] * mas[1][1] - (mas[1][0] * mas[0][1]);
        return (d);
    }
    if (m>2) {
        for (i = 0; i<m; i++) {
            GetMatr(mas, p, i, 0, m);
            d = d + k * mas[i][0] * Determinant(p, n);
            k = -k;
        }
    }
    return (d);
}

int is_matrix_positive(double** A, double** A_t, int N)
{

    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            A_t[i][j] = A[i][j] + A_t[i][j];
        }
    }

    auto* b = new double[N]();
    print_matrix(A_t, b, N);

    double det;
    for(int i = 1; i <= N; i++)
    {
        det = Determinant(A_t, i);
        if(det <= 0)
        {
            cout << "Corner Minors less than zero" << endl;
            exit(1);
        }
    }
    return 0;
}

void MatrVekt(int N, double **M, const double *V, double *R)
//N- размерность, M- матрица, V- вектор, R- результат
{
    for(int i=0; i<N; i++)
    {
        R[i]=0;
        for(int j=0; j<N; j++)
            R[i]+= M[i][j]*V[j];
    }
}

int MinNev(int N, double **A, const double *b, double *X, double eps)
//N- размерность, A- матрица, F- вектор свободных членов, X-вектор результат eps- точность
{
    int count=0;//  количество итераций
    auto *R = new double [N];
    auto *Delta = new double [N]; // delta = Ax^k - b
    auto *TempX = new double[N](); // x^k
    double maxi, Tau, TempTau;

    do
    {
        MatrVekt(N, A, TempX, R); //произведение матрицы коэффициентов на вектор приближения: R = A*x^k
        for(int i=0; i<N; i++)
        {
            Delta[i] = R[i] - b[i]; //Вектор невязок: delta = Ax^k - b
        }
        MatrVekt(N, A, Delta, R); // R = A*Delta

        Tau=0.0; // наш коэф
        TempTau=0.0;
        for(int i=0; i<N; i++) // скалярное произведение векторов
        {
            Tau+=R[i]*Delta[i]; // (Ar^k, r^k)
            TempTau+=R[i]*R[i]; // (Ar^k, Ar^k)
        }
        Tau=Tau/TempTau; // (Ar^k, r^k) / (Ar^k, Ar^k) = t

        for(int i=0; i<N; i++) // x^k+1 = x^k - t*delta
            X[i]=TempX[i]-Tau*Delta[i];

        maxi = abs(X[0] - TempX[0]); // max = |x^k+1[0] - x^k[0]|
        for(int i=0; i<N; i++)
        {
            if(abs(X[i]-TempX[i])>maxi) //вычисление текущего приближения
                maxi=abs(X[i]-TempX[i]);
            TempX[i]=X[i];
        }
        count++;
    }
    while (maxi>=eps);

    delete[] R;
    delete[] Delta;
    delete[] TempX;
    return count;
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
    double** A, **A_t;
    double* b;
    double* x;
    int N;
    int it;
    double eps = 0.0000000001;
    cin >> N;
    x = new double[N]();
    A = new double*[N];
    A_t = new double*[N];
    for(int i = 0; i < N; i++)
    {
        A_t[i] = new double[N];
        A[i] = new double[N];
    }
    b = new double[N]; // СТОЛБИК СВОБОДНЫХ
    init_matrix(A,N);
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            A_t[i][j] = A[i][j];
        }
    }
    transpos(A_t, N);
    init_vector(b,N);
    is_matrix_positive(A, A_t, N);
    print_matrix(A,b,N);
    it = MinNev(N, A, b, x, eps);

    cout << "The final answer is:" << endl;
    for(int i = 0; i < N; i++ )
    {
        cout << "x[" << i+1 << "] = " << x[i] << endl;
    }

    cout << endl << it;
    return 0;
}