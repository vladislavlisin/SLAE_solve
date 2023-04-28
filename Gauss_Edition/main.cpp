#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

double Determinant(double **mas, int m);
void make_triangle(double** a, double* b, int n);
void find_x(double **a, const double *b, double *x, int n);
void print_matrix(double ** A, double* b, int n);
double find_det_diagM(double** A, int n);
void find_inverse_m(double** matrix, int order);

int main()
{
    int n;
    cout << "Enter matrix's dimension" << endl;
    cin >> n;
    double **A, **A1;
    double *b;
    double *x;
    double det, det1;
    x = new double[n];
    b = new double[n];
    A = new double *[n]; // матрица n*n
    A1 = new double *[2*n]; // матрица 2n*2n

    for (int i = 0; i < n; i++) {
        A[i] = new double[n];
    }
    for (int i = 0; i < 2*n; i++) {
        A1[i] = new double[2*n];
    }

    for (int q = 0; q < 2*n; q++) // заполняем A1[] нулями
    {
        for (int c = 0; c < 2*n; c++) {
            A1[q][c] = 0;
        }
    }

    for (int i = 0; i < n; i++) // инициализируем x[]
    {
        x[i] = 0;
    }
    for(int i = 0; i < n; i++) // вводим A[]
    {
        for (int j = 0; j < n; j++)
        {
            cin >> A[i][j];
            A1[i][j] = A[i][j];
        }
    }
    det1 = Determinant(A, n);
    for(int i = 0; i < n; i++) // вводим b[]
    {
        cin >> b[i];
    }

    make_triangle(A, b, n); // приводим к верхнетреуг виду матрицу A
    det = find_det_diagM(A, n);
    print_matrix(A, b, n); // печатать расширенную матрицу A
    find_x(A,b,x,n); // НАХОДИМ РЕШЕНИЕ СЛАУ
    cout << "Roots: " << endl;
    for(int q = 0; q < n; q++)
    {
        cout << "x" << q+1 << " = " << fixed << setprecision(8) << x[q] << endl;
    }
    find_inverse_m(A1, n); // НАХОДИМ ОБРАТ
    det1 = Determinant(A1, n);
    cout << "Det = " << fixed << setprecision(8) << det1 << endl;
    return 0;
}

void make_triangle(double** a, double* b, int n) // схема единственного деления, прямой ход
{
    double v;
    for(int k = 0,i,j,im; k < n - 1; k++)
    {
        im = k;
        for(i = k + 1; i < n; i++)
        {
            if(fabs(a[im][k]) < fabs(a[i][k]))
            {
                im = i;
            }
        }
        if(im != k)
        {
            for(j = 0; j < n; j++)
            {
                v = a[im][j];
                a[im][j] = a[k][j];
                a[k][j] = v;
            }
            v = b[im];
            b[im] = b[k];
            b[k]  = v;
        }
        for(i = k + 1; i < n; i++)
        {
            v = a[i][k]/a[k][k];
            a[i][k] = 0;
            b[i] = b[i] - v*b[k];
            if(v != 0)
                for(j = k + 1; j < n; j++)
                {
                    a[i][j] = a[i][j] - v*a[k][j];
                }
        }
    }
}

void find_x(double **a, const double *b, double *x, int n)
{
    double s;
    x[n - 1] = b[n - 1]/a[n - 1][n - 1];
    for(int i = n - 2, j; 0 <= i; i--)
    {
        s = 0;
        for(j = i + 1; j < n; j++)
        {
            s = s+a[i][j]*x[j];
        }
        x[i] = (b[i] - s)/a[i][i];
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

// Функция для печати обратной матрицы

void PrintMatrix(double** ar, int n, int m)
{
    for (int q = 0; q < n; q++) // вывод матрицы A
    {
        for (int c = 0; c < m; c++) {
            cout << fixed << setprecision(8) << ar[q][c] << " ";
        }
        cout << endl;
    }
}

// Функция для выполнения обратной операции над матрицей.
void find_inverse_m(double **matrix, int order)
{

    printf("[A]\n");
    PrintMatrix(matrix, order, order);

    // Создаем матрицу [A|E]
    for (int i = 0; i < order; i++)
    {
        for (int j = 0; j < 2 * order; j++)
        {
            if (j == (i + order)) {
                matrix[i][j] = 1;
            }
        }
    }

    double temp;
    // Преобразования
    for (int i = order - 1; i > 0; i--)
    {
        if (matrix[i - 1][0] < matrix[i][0])
            for(int j = 0; j < 2 * order; j ++)
            {
                temp = matrix[i][j];
                matrix[i][j] = matrix[i - 1][j];
                matrix[i - 1][j] = temp;
            }
        if (matrix[i - 1][0] < matrix[i][0])
        {
            double *temp1 = matrix[i];
            *matrix[i] = *matrix[i - 1];
            *matrix[i - 1] = *temp1;
        }
    }

    // Печать матрицы после преобразований.
    printf("\n[A|E]\n");
    PrintMatrix(matrix, order, 2*order);

    // приведение к диагональному виду слева
    for (int i = 0; i < order; i++)
    {
        for (int j = 0; j < order; j++)
        {
            if (j != i)
            {
                temp = matrix[j][i] / matrix[i][i];
                for (int k = 0; k < 2 * order; k++)
                {
                    matrix[j][k] -= matrix[i][k] * temp;
                }
            }
        }
    }

    // Делим элемент строки по диагональному элементу
    for (int i = 0; i < order; i++)
    {
        temp = matrix[i][i];
        for (int j = 0; j < 2 * order; j++)
        {
            matrix[i][j] = matrix[i][j] / temp;
        }
    }

    // выводим обратную матрицу.
    printf("\n[A^-1]\n");
    PrintMatrix(matrix, order, 2 * order);
}

double find_det_diagM(double** A, int n)
{
    double det = 1;
    for(int i = 0; i < n; i++)
    {
        det = det*A[i][i];
        if(abs(A[i][i]) <= 0.000000001)
        {
            cout << "No solution!!!" << endl;
            exit(1);
        }
    }
    return det;
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


void solve_3diag_Matrix (int n, double *a, double *c, double *b, double *f, double *x)
{
    double m;
    for (int i = 1; i < n; i++)
    {
        m = a[i]/c[i-1];
        c[i] = c[i] - m*b[i-1];
        f[i] = f[i] - m*f[i-1];
    }

    x[n-1] = f[n-1]/c[n-1];

    for (int i = n - 2; i >= 0; i--)
    {
        x[i]=(f[i]-b[i]*x[i+1])/c[i];
    }
}