#include <iomanip>
#include <iostream>
#include <vector>
#include <map>
#include <stdio.h>

typedef unsigned long long tn;

using namespace std;

void printList(const vector<double>& numbers)
{
	constexpr int precision = 4;

	cout.setf(ios::fixed, ios::floatfield);
	cout.precision(precision);

	cout << "{ ";
	for (tn i = 0; i < numbers.size() - 1; i++)
		cout << numbers[i] << ", ";
	cout << numbers[numbers.size() - 1] << " }\n";
}
void printHistogram(const vector<double>& numbers, double start, double stopExcluding, int subdivisionsNum = 10)
{
	constexpr int precision = 2;
	const int width = max(log10(stopExcluding), start < 0 ? log10(-start) : 0) + precision + 2;

	cout.setf(ios::fixed, ios::floatfield);
	cout.precision(precision);
	
	cout << setw(2 * width + 5) << "Interval    |" << "  Frequency \n";

	int buff, l = numbers.size();
	double s, e;
	for(int i = 0; i < subdivisionsNum; i++)
	{
		buff = 0;
		s = (i * stopExcluding + (subdivisionsNum - i) * start) / subdivisionsNum;
		e = ((i + 1) * stopExcluding + (subdivisionsNum - i - 1) * start) / subdivisionsNum;

		cout << "[" << setw(width) << s << ";" << setw(width) << e << ")     ";

		for (int j = 0; j < l; j++)
			if (numbers[j] >= s && numbers[j] < e)
				buff++;
		cout << (double)buff / l << "\n";
	}
}

// const tn n, tn X0 = 2147483,
// const tn m = 2147483648, const tn c = 2147483647, const tn a = 2147483637
vector<double> method1(
	const tn n, tn X0 = 42949672,
	const tn m = 4294967291, const tn c = 4294967279, const tn a = 4294967231)
{
	vector<double> ret(n);
	for(tn i = 0; i < n; i++)
	{
		X0 = (a * X0 % m + c) % m;
		ret[i] = (double)X0 / m;
	}
	return ret;
}

vector<double> method2(
	const tn n, tn X0 = 42949672,
	const tn m = 4294967291, const tn c = 4294967279, const tn d = 4294967197, const tn a = 4294967231)
{
	vector<double> ret(n);
	for (tn i = 0; i < n; i++)
	{
		X0 = (d * (X0 * X0 % m) % m + c) % m;
		ret[i] = (double)X0 / m;
	}
	return ret;
}

int main()
{
	printList(method2(100));
	printHistogram(method2(100), 0, 1);
}