#include <omp.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;
const double eps = 0.000001;
struct crsmtx 
{
	int n, m;
	int notzeros;
	double* values;
	int* cols;
	int* pointer;

	crsmtx(crsmtx& mtx) 
	{
		this->n = mtx.n;
		this->m = mtx.m;
		this->notzeros = mtx.notzeros;
		this->cols = new int[this->notzeros];
		this->values = new double[this->notzeros];
		this->pointer = new int[this->n + 1];
		for (int i = 0; i < this->notzeros; ++i)
		{
			this->cols[i] = mtx.cols[i];
			this->values[i] = mtx.values[i];
		}
		for (int i = 0; i < n + 1; ++i)
		{
			this->pointer[i] = mtx.pointer[i];
		}
	}
	crsmtx() {}
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

void setRandomMtx(int n, int m, crsmtx &mtx)
{
	double a;
	vector<double> vals;
	vector<int> cols;
	vector<int> pointer(n + 1, 0);
	//initmtx(n, m, nz, mtx);
	int cnt = 0;
	//mtx.pointer[0] = 0;
	srand(time(NULL));
	double p;
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			p = rand() % 100 + 1;
			if (p <= 5) 
			{
				vals.push_back((double)rand() / RAND_MAX * 100);
				cols.push_back(j);
				cnt++;
			}
		}
		pointer[i + 1] = cnt;
	}

	initmtx(n, m, vals.size(), mtx);

	for (int i = 0; i < vals.size(); ++i) 
	{
		mtx.values[i] = vals[i];
		mtx.cols[i] = cols[i];
	}

	for (int i = 0; i < pointer.size(); ++i)
	{
		mtx.pointer[i] = pointer[i];
	}

}

void setmtxfromfileFast(char* file, crsmtx &mtx)
{
	double a;
	int n, m, nz;
	FILE* f;
	char c;
	
	f = fopen(file, "r");
	//ifstream fin;
	//fin.open(file);
	//fin >> n >> m >> nz;
	fscanf(f, "%d %d %d", &n, &m, &nz);
	initmtx(n, m, nz, mtx);
	int cnt = 0;
	mtx.pointer[0] = 0;
	//fscanf(f, "%c", &c);
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			///fin >> a;
			fscanf(f, "%lf", &a);
			//fscanf(f, "%c", &c);

			if (abs(a) > eps)
			{
				mtx.values[cnt] = a;
				mtx.cols[cnt] = j;
				cnt++;
			}
			

		}
		//fscanf(f, "%c", &c);
		mtx.pointer[i + 1] = cnt;
	}
	//fin.close();
	fclose(f);

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

void getMtxTemplate(crsmtx a, crsmtx b, crsmtx &c) 
{
	int **templ = new int*[a.n];
	int nz2 = 0;
	for (int i = 0; i < a.n; ++i)
	{
		int *x = new int[a.m];
		templ[i] = new int[a.m];

		memset(x, 0, sizeof(int)*a.m);
		memset(templ[i], sizeof(int) * 0, a.m);

		for (int j = a.pointer[i]; j < a.pointer[i + 1]; ++j) 
		{
				x[a.cols[j]] = 1;
		}

		for (int j = 0; j < b.n; ++j) 
		{
			int nz = 0;

			for (int q = b.pointer[j]; q < b.pointer[j + 1]; ++q)
			{
				if (x[b.cols[q]] == 1) 
				{
					nz++;
					break;
				}
			}
			if (nz > 0)
			{
				templ[i][j] = 1;
				nz2++;
			}
				
		}
		delete[] x;
	}
	initmtx(a.n, a.m, nz2, c);
	
	int cnt = 0;
	
	
	for (int i = 0; i < a.n; ++i) 
	{
		for (int j = 0; j < a.m; ++j) 
		{
			if (templ[i][j] >  0)
			{
				c.cols[cnt] = j;
				c.values[cnt] = 1;
				cnt++;
			}
			
		}
		c.pointer[i + 1] = cnt;
		//delete[] templ[i];
	}

	if (cnt != nz2) cout << "fail!" << endl;

	for (int i = 0; i < a.n; ++i)
		delete templ[i];

	delete[] templ;

}

void mtxMultiplyNew(crsmtx a, crsmtx b, crsmtx &c, double& time) 
{
	clock_t start = clock();

	for (int i = 0; i < c.n; ++i) 
	{
		int *x = new int[a.m];
		memset(x, -1, sizeof(int)*a.m);

		for (int j = a.pointer[i]; j < a.pointer[i + 1]; ++j)
		{
			x[a.cols[j]] = j;
		}

		for (int j = c.pointer[i]; j < c.pointer[i + 1]; ++j) 
		{
			double sum = 0.0;

			int col = c.cols[j];
			for (int q = b.pointer[col]; q < b.pointer[col + 1]; ++q)
			{
				int ind = b.cols[q];
				if (x[ind] != -1)
					sum += a.values[x[ind]] * b.values[q];
			}
			c.values[j] = sum;
		}
		delete[] x;
	
	}

	clock_t finish = clock();
	time = (double)(finish - start) / CLOCKS_PER_SEC;

}

void mtxMultiplyNewParallel(crsmtx a, crsmtx b, crsmtx &c, double& time)
{
	clock_t start = clock();

#pragma omp parallel for
	for (int i = 0; i < c.n; ++i)
	{
		int *x = new int[a.m];
		memset(x, -1, sizeof(int)*a.m);

		for (int j = a.pointer[i]; j < a.pointer[i + 1]; ++j)
		{
			x[a.cols[j]] = j;
		}

		for (int j = c.pointer[i]; j < c.pointer[i + 1]; ++j)
		{
			double sum = 0.0;

			int col = c.cols[j];
			for (int q = b.pointer[col]; q < b.pointer[col + 1]; ++q)
			{
				int ind = b.cols[q];
				if (x[ind] != -1)
					sum += a.values[x[ind]] * b.values[q];
			}
			c.values[j] = sum;
		}
		delete[] x;

	}

	clock_t finish = clock();
	time = (double)(finish - start) / CLOCKS_PER_SEC;

}

int main() 
{
	crsmtx a;
	crsmtx b;
	crsmtx bt;
	//setmtxfromfile("Text.txt", a);
	//setmtxfromfile("Text.txt", b);

	//setmtxfromfileFast("bigA.txt", a);
	//setmtxfromfileFast("bigB.txt", b);


	setRandomMtx(10000, 10000, a);
	setRandomMtx(10000, 10000, b);



	crsmtx c;
	transposemtx(b, bt);
	double time1, time2;

	cout << "n = " << a.n << endl;
	cout << "m = " << a.m << endl;



	//mtxmultiply(a, bt, c, time1);
	//mtxmultiplyparallel(a, bt, cp, time2);


	getMtxTemplate(a, bt, c);
	//getMtxTemplate(a, bt, cp);
	
	crsmtx cp(c);

	cout << "gets template!" << endl;

	

	mtxMultiplyNew(a, bt, c, time1);
	mtxMultiplyNewParallel(a, bt, cp, time2);

	

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