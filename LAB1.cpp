#include <iomanip>
#include <iostream>
#include <vector>
#include <map>
#include <stdio.h>

typedef long long tn;

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

//vector<vector<double>(*)()>

// METHOD1
// 31)
// const tn n, tn X0 = 2147483,
// const tn m = 2147483648, const tn c = 2147483647, const tn a = 2147483637
//
// 32)
// const tn n, tn X0 = 42949672,
// const tn m = 4294967296, const tn c = 4294967291, const tn a = 4294967157
//
// const tn m = 4294967291, const tn c = 4294967279, const tn a = 4294967231

vector<double> method1(
	const tn n, tn X0 = 42949672,
	const tn m = 4294967296, const tn c = 4294967291, const tn a = 4294967157)
{
	vector<double> ret(n);
	for(tn i = 0; i < n; i++)
	{
		X0 = (a * X0 % m + c) % m;
		ret[i] = (double)X0 / m;
	}
	return ret;
}

// METHOD2
// 32)
// const tn n, tn X0 = 42949672,
// const tn m = 4294967296, const tn c = 4294967291, const tn d = 4294967156, const tn a = 4294967157

vector<double> method2(
	const tn n, tn X0 = 42949672,
	const tn m = 4294967291, const tn c = 4294967279, const tn d = 4294967231, const tn a = 4294967197)
{
	vector<double> ret(n);
	for (tn i = 0; i < n; i++)
	{
		X0 = ( (d*(X0*X0%m) % m + a*X0%m) % m + c) % m;
		ret[i] = (double)X0 / m;
	}
	return ret;
}

vector<double> method3(
	const tn n, tn X0 = 1247437, tn X1 = 224743647, const tn m = 4294967291)
{
	vector<double> ret(n);
	tn buff;
	for (tn i = 0; i < n; i++)
	{
		buff = (X1 + X0) % m;
		X0 = X1; X1 = buff;
		ret[i] = (double)buff / m;
	}
	return ret;
}

tn inverse(tn x, tn m)
{
	tn m0 = m;
	tn a = 0, b = 1, buff;
	while (x > 1 && m > 0) {
		buff = a;
		a = b - x / m * a;
		b = buff;

		buff = m;
		m = x % m;
		x = buff;
	}
	if (b < 0)
		b += m0;
	return b;
}
// METHOD4
// 32)
// const tn n, tn X0 = 42949672,
// const tn m = 4294967296, const tn c = 4294967290, const tn a = 4294967157
vector<double> method4(
	const tn n, tn X0 = 42949672,
	const tn m = 4294967197, const tn c = 4294967291, const tn a = 4294967157)
{
	vector<double> ret(n);
	tn buff;
	for (tn i = 0; i < n; i++)
	{
		X0 = (a * inverse(X0, m) % m + c) % m;
		ret[i] = (double)X0 / m;
	}
	return ret;
}

int main()
{
	printList(method4(100));
	printHistogram(method4(10000), 0, 1);
}