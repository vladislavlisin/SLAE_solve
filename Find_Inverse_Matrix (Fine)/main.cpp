#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;

void arrIn(double** arr, int n)
{
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) cin >> arr[i][j];
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
double** inv(double** a, int num);

int main()
{
    int n;
    double eps = 0.00000001;
    cin >> n;
    auto arr = new double* [n];
    auto inv_arr = new double* [n];
    for (int i = 0; i < n; i++)
    {
        arr[i] = new double[n];
        inv_arr[i] = new double[n];
    }

    arrIn(arr, n);
    inv_arr = inv(arr, n);
    arrOut(inv_arr, n);


    return 0;
}

void swap(double* a, double* b); //declare subroutine
double** inv(double** a, int num)
{
    int* is, * js, i, j, k;
    int n = num;
    double temp, fmax;
    double** tp = new double* [num];
    for (int i = 0; i < num; i++) tp[i] = new double[num];
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
            return NULL;
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
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
        }
    }
    return tp;
}
void swap(double* a, double* b)
{
    double c;
    c = *a;
    *a = *b;
    *b = c;
}