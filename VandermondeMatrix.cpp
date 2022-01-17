#include <iostream>
#include <algorithm>
#include <iterator>
#include <iomanip>
#include <vector>
#include <cmath>
#include <math.h>
using namespace std;

double three_sig(double N) // returns the given number rounded to 3 significant digits
{
    double a, b, i;

    // counting the number of digits to the left of decimal point
    a = N;
    for (i = 0; a >= 1; ++i)
    {
        a /= 10;
    }

    // modifies value for rounding
    a = N * pow(10, 3 - i);
    b = a + 0.5;

    // checking whether to round up
    if (b == ceil(a))
    {
        b--;
    }

    return floor(b) / pow(10, 3 - i);
}

void print_matrix(vector<vector<double>> M) // displays formatted matrix
{
    for (int i = 0; i < M.size(); i++)
    {
        for (int j = 0; j < M[i].size(); j++)
        {
            cout << setw(7) << three_sig(M[i][j]) << " ";
        }
        cout << "\n";
    }
}

vector<vector<double>> transpose(vector<vector<double>> M) //transposes matrix
{
    vector<vector<double>> trans(M[0].size(), vector<double>(M.size()));

    for (int i = 0; i < M[0].size(); i++)
    {
        for (int j = 0; j < M.size(); j++)
        {
            trans[i][j] = M[j][i];
        }
    }

    return trans;
}

double dot_product(vector<double> v1, vector<double> v2) // returns dot product of 2 vectors
{
    if (v1.size() == v2.size()) // checks condition to avoid error
    {
        double sum = 0;
        for (int i = 0; i < v1.size(); i++)
        {
            sum += v1[i] * v2[i];
        }
        return sum;
    }
    return -1;
}

vector<vector<double>> product(vector<vector<double>> M1, vector<vector<double>> M2) //calculates product of 2 matrices
{
    vector<vector<double>> P(M1.size(), vector<double>(M2[0].size())); //product matrix
    vector<vector<double>> M2T = transpose(M2);                        // using transpose to multiply vectors within second matrix

    for (int i = 0; i < M1.size(); i++)
    {
        for (int j = 0; j < M2[0].size(); j++)
        {
            P[i][j] = dot_product(M1[i], M2T[j]);
        }
    }

    return P;
}

double determinant(vector<vector<double>> M) //calculates determinant of matrix
{
    int size = M.size();

    if (size == 1) // base cases
    {
        return M[0][0];
    }
    else if (size == 2)
    {
        return M[0][0] * M[1][1] - M[0][1] * M[1][0];
    }
    else if (size == 3)
    {
        return M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
    }

    double det = 0;

    //creates minor matrix
    for (int i = 0; i < size; i++)
    {
        vector<vector<double>> temp(size - 1, vector<double>(size - 1));
        for (int m = 1; m < size; m++)
        {
            int z = 0;
            for (int n = 0; n < size; n++)
            {
                if (n != i) //only adds elements not in the same column
                {
                    temp[m - 1][z] = M[m][n];
                    z++;
                }
            }
        }

        //recursive call adds determinant of minor matrix to determinant
        det += pow(-1, i) * M[0][i] * determinant(temp);
    }

    return det;
}

vector<vector<double>> cofactor(vector<vector<double>> M) //calculates cofactor matrix
{
    vector<vector<double>> C(M.size(), vector<double>(M.size()));
    vector<vector<double>> temp(M.size() - 1, vector<double>(M.size() - 1));

    // creates minor matrix for each elements and fills in cofactor matrix with determinant of each
    for (int i = 0; i < M.size(); i++)
    {
        for (int j = 0; j < M[0].size(); j++)
        {
            int p = 0;
            for (int x = 0; x < M.size(); x++)
            {
                if (x != i)
                {
                    int q = 0;
                    for (int y = 0; y < M.size(); y++)
                    {
                        if (y != j)
                        {
                            temp[p][q] = M[x][y];
                            q++;
                        }
                    }
                    p++;
                }
            }
            C[i][j] = pow(-1, i + j) * determinant(temp);
        }
    }
    return C;
}

vector<vector<double>> inverse(vector<vector<double>> M) //calculates inverse of matrix
{
    int det = determinant(M);
    vector<vector<double>> inv = transpose(cofactor(M)); //adjoint matrix

    for (int i = 0; i < M.size(); i++)
    {
        for (int j = 0; j < M.size(); j++)
        {
            inv[i][j] /= det;
        }
    }

    return inv;
}

int main()
{
    // data vectors
    vector<double> x;
    vector<double> y;
    int size;
    int order;
    double input;

    // user input for data
    cout << "Enter number of data points: ";
    cin >> size;

    cout << "Enter all x values:\n";
    for (int i = 0; i < size; i++)
    {
        cin >> input;
        x.push_back(input);
    }

    cout << "Enter all y values:\n";
    for (int i = 0; i < size; i++)
    {
        cin >> input;
        y.push_back(input);
    }
    
    cout << "Enter the order of polynomial fit:\n";
    cin >> order;

    vector<vector<double>> V(size, vector<double>(order + 1)); // Vandermonde Matrix

    // filling in Vandermonde Matrix with cubic values of data
    for (int i = 0; i < order + 1; i++)
    {
        for (int j = 0; j < size; j++)
        {
            V[j][i] = pow(x[j], order - i);
        }
    }

    // Transpose of Vandermonde Matrix
    vector<vector<double>> T = transpose(V);

    // Product of Transpose and Vandermonde Matrices
    vector<vector<double>> TV = product(T, V); //Product of Transpose and Vandermonde Matrices

    //Inverse of Product
    vector<vector<double>> I = inverse(TV);

    //converting y to a column vector
    vector<vector<double>> columny(size, vector<double>(1));
    for (int i = 0; i < size; i++)
    {
        columny[i][0] = y[i];
    }

    //Coefficient Vector
    vector<vector<double>> a = product(product(I, T), columny);

    // displaying data and results
    cout << "Data:\nx: ";
    copy(x.begin(), x.end(), ostream_iterator<double>(cout, " "));
    cout << "\ny: ";
    copy(y.begin(), y.end(), ostream_iterator<double>(cout, " "));
    cout << "\nVandermonde Matrix:\n";
    print_matrix(V);
    cout << "Transpose Vandermonde Matrix:\n";
    print_matrix(T);
    cout << "Product of Transpose and Vandermonde Matrix:\n";
    print_matrix(TV);
    cout << "Inverse of Product:\n";
    print_matrix(I);
    cout << "Polynomial Coefficients:\n";
    print_matrix(a);
}