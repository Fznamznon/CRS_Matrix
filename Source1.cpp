//#include <fstream>
//#include <stdlib.h>
//#include <time.h>
//#include <iostream>
//
//using namespace std;
//
//int n = 10000;
//int m = 10000;
//
//int nz = 0;
//
//int main() 
//{
//	ofstream fo;
//	fo.open("bigB.txt");
//	fo << n << " " << m << " " << nz << endl;
//	srand(time(NULL));
//	double p;
//	double val;
//	for (int i = 0; i < n; ++i) 
//	{
//		for (int j = 0; j < m; ++j) 
//		{
//			p = rand() % 100 + 1;
//			if (p <= 5) 
//			{
//				val = (double)rand() / RAND_MAX * 100;
//				nz++;
//			}
//			else
//			{
//				val = 0;
//			}
//			fo << val << " ";
//		}
//		fo << endl;
//	}
//
//	fo.close();
//
//	cout << nz << endl;
//
//	return 0;
//}