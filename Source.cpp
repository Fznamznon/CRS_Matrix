#include <omp.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <time.h>

using namespace std;
const double eps = 0.000001;
struct crsmtx 
{
	int n, m;
	int notzeros;
	double* values;
	int* cols;
	int* pointer;
};

void initmtx(int n, int m, int nz, crsmtx &mtx) 
{
	mtx.n = n;
	mtx.m = m;
	mtx.notzeros = nz;
	mtx.values = new double[nz];
	mtx.cols = new int[nz];
	mtx.pointer = new int[n + 1];
	for (int i = 0; i < nz; ++i) 
	{
		mtx.cols[i] = 0;
		mtx.values[i] = 0.0;
 	}
	for (int i = 0; i < n + 1; ++i)
	{
		mtx.pointer[i] = 0;
	}
}

void deletemtx(crsmtx &mtx) 
{
	delete[] mtx.values;
	delete[] mtx.cols;
	delete[] mtx.pointer;
}

void setmtxfromfile(string file, crsmtx &mtx) 
{
	double a;
	int n, m, nz;
	ifstream fin;
	fin.open(file);
	fin >> n >> m >> nz;
	initmtx(n, m, nz, mtx);
	int cnt = 0;
	mtx.pointer[0] = 0;
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			fin >> a;
			if (abs(a) > eps)
			{
				mtx.values[cnt] = a;
				mtx.cols[cnt] = j;
				cnt++;
			}
		}
		mtx.pointer[i + 1] = cnt;
	}
	fin.close();


}

int comparemtx(crsmtx a, crsmtx b, double& diff)
{
	if (a.n != b.n || a.m != b.m)
		return 1;
	int n = a.n;
	int m = a.m;
	
	vector<vector<double>> p(n, vector<double>(m));


	for (int i = 0; i < n; ++i) 
	{
		for (int j = a.pointer[i]; j < a.pointer[i + 1]; ++j) 
		{
			p[i][a.cols[j]] = a.values[j];
		}
	}

	
	double max = -1;
	double tmp, tmp2;
	for (int i = 0; i < n; ++i)
	{
		for (int j = b.pointer[i]; j < b.pointer[i + 1]; ++j)
		{
			tmp = abs(p[i][b.cols[j]] - b.values[j]);
			if (tmp > max) max = tmp;
		}
	}

	diff = max;
	return 0;
}

void transposemtx(crsmtx &a, crsmtx &at) 
{
	int n = a.n;
	int m = a.m;
	int nz = a.notzeros;
	initmtx(m, n, nz, at);

	for (int i = 0; i < nz; ++i) 
	{
		at.pointer[a.cols[i] + 1]++;
	}

	int ind = 0;
	int tmp;
	for (int i = 1; i < m + 1; ++i) 
	{
		tmp = at.pointer[i];
		at.pointer[i] = ind;
		ind += tmp;
	}
	int col, row;
	double val;
	for (int i = 0; i < n; ++i) 
	{
		col = i;
		for (int j = a.pointer[i]; j < a.pointer[i + 1]; ++j) 
		{
			val = a.values[j];
			row = a.cols[j];
			ind = at.pointer[row + 1];
			at.values[ind] = val;
			at.cols[ind] = col;
			at.pointer[row + 1]++;
		}
	}

}

int mtxmultiply(crsmtx &a, crsmtx &b, crsmtx &c, double &time) 
{
	if (a.n != b.n)
		return 0;

	vector<double> values;
	vector<int> cols;
	vector<int> pointer;
	double s = 0;
	int ind;
	int nz = 0;

	clock_t start = clock();


	for (int i = 0; i < a.n; ++i) 
	{
		vector<int> x(a.m, -1);

		for (int j = a.pointer[i]; j < a.pointer[i + 1]; ++j) 
		{
			x[a.cols[j]] = j;
		}
		
		for (int j = 0; j < b.n; ++j) 
		{
			s = 0.0;
			for (int q = b.pointer[j]; q < b.pointer[j + 1]; ++q) 
			{
				ind = b.cols[q];
				if (x[ind] >= 0) 
				{
					s += a.values[x[ind]] * b.values[q];
				}
			}
			if (abs(s) > eps) 
			{
				values.push_back(s);
				cols.push_back(j);
				nz++;
			}
		}
		pointer.push_back(nz);
	}

	initmtx(a.n, b.m, nz, c);

	for (int i = 0; i < nz; ++i) 
	{
		c.values[i] = values[i];
		c.cols[i] = cols[i];
	}
	c.pointer[0] = 0;
	for (int i = 0; i < c.n; ++i)
	{
		c.pointer[i + 1] = pointer[i];
	}
	clock_t finish = clock();
	time = (double)(finish - start) / CLOCKS_PER_SEC;

	return 0;
}

int mtxmultiplyparallel(crsmtx &a, crsmtx &b, crsmtx &c, double &time)
{
	if (a.n != b.n)
		return 0;

	vector<vector<double>> values(a.n);
	vector<vector<int>> cols(a.n);
	vector<int> pointer(a.n + 1);

	int nz = 0;
	clock_t start = clock();

	int i = 0;
#pragma omp parallel for
	for (i = 0; i < a.n; ++i)
	{
		vector<int> x(a.m, -1);

		for (int j = a.pointer[i]; j < a.pointer[i + 1]; ++j)
		{
			x[a.cols[j]] = j;
		}

		for (int j = 0; j < b.n; ++j)
		{
			double s = 0.0;
			for (int q = b.pointer[j]; q < b.pointer[j + 1]; ++q)
			{
				int ind = b.cols[q];
				if (x[ind] != -1)
				{
					s += a.values[x[ind]] * b.values[q];
				}
			}
			if (abs(s) > eps)
			{
				values[i].push_back(s);
				cols[i].push_back(j);
				pointer[i]++;
			}
		}
	}
	

	for (int i = 0; i < a.n; ++i) 
	{
		int tmp = pointer[i];
		pointer[i] = nz;
		nz += tmp;
	}
	pointer[a.n] = nz;

	initmtx(a.n, b.m, nz, c);
	int k = 0;
	for (int i = 0; i < a.n; ++i)
	{
		for (int j = 0; j < values[i].size(); ++j)
		{
			c.values[k] = values[i][j];
			c.cols[k] = cols[i][j];
			k++;
		}
	}
	//c.pointer[0] = 0;
	for (int i = 0; i <= c.n; ++i)
	{
		c.pointer[i] = pointer[i];
	}
	clock_t finish = clock();
	time = (double)(finish - start) / CLOCKS_PER_SEC;
	return 0;
}


void printmtx(crsmtx mtx) 
{
	int n = mtx.n;
	int m = mtx.m;
	vector<vector<double>> p(mtx.n, vector<double>(mtx.m));


	for (int i = 0; i < n; ++i)
	{
		for (int j = mtx.pointer[i]; j < mtx.pointer[i + 1]; ++j)
		{
			p[i][mtx.cols[j]] = mtx.values[j];
		}
	}
	for (int i = 0; i < n; i++) 
	{
		for (int j = 0; j < m; ++j) 
		{
			cout << p[i][j] << " ";
			
		}
		cout << endl;
	}
}

int main() 
{
	crsmtx a;
	crsmtx b;
	crsmtx bt;
	setmtxfromfile("bigA.txt", a);
	setmtxfromfile("bigB.txt", b);
	crsmtx c, cp;
	transposemtx(b, bt);
	double time1, time2;

	cout << "n = " << a.n << endl;
	cout << "m = " << a.m << endl;



	mtxmultiply(a, bt, c, time1);
	mtxmultiplyparallel(a, bt, cp, time2);


	/*for (int i = 0; i < c.notzeros; ++i) 
	{
		cout << c.values[i] << " ";
	}
	cout << endl;
	for (int i = 0; i < c.notzeros; ++i)
	{
		cout << c.cols[i] << " ";
	}
	cout << endl;

	for (int i = 0; i < cp.notzeros; ++i)
	{
		cout << cp.values[i] << " ";
	}
	cout << endl;

	for (int i = 0; i < c.notzeros; ++i)
	{
		cout << cp.cols[i] << " ";
	}
	cout << endl;*/

	//printmtx(c);
	////cout << endl;
	//printmtx(cp);


	double diff = -1;
	//comparemtx(c, cp, diff);
	cout << "diff = " << diff;
	cout << endl;
	cout << "linear time: " << time1 << endl;
	cout << "parallel time: " << time2 << endl;


	deletemtx(a);
	deletemtx(b);
	deletemtx(bt);
	deletemtx(c);
	deletemtx(cp);



	return 0;
}