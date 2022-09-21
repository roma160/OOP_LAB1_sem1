#include <iomanip>
#include <iostream>
#include <vector>
#include <map>
#include <stdio.h>

using namespace std;

typedef long long tn;
typedef vector<double>(*method_type)(tn n);

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
	const int width = max(log10(stopExcluding), start < 0 ? log10(-start) : 0) + precision + 2 + (start < 0);

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

		cout.flush();
		cout << "[" << setw(width) << s << ";" << setw(width) << e << ")     ";

		for (int j = 0; j < l; j++)
			if (numbers[j] >= s && numbers[j] < e)
				buff++;
		cout << (double)buff / l << "\n";
	}
}



vector<tn> default_method(int i, tn n);

vector<double> U(const vector<tn> S, const tn m)
{
	tn n = S.size();
	vector<double> ret(n);
	for (tn i = 0; i < n; i++)
		ret[i] = (double)S[i] / (m-1);
	return ret;
}

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
vector<tn> method1(
	const tn n, const tn X0 = 42949672,
	const tn m = 4294967296, const tn c = 4294967291, const tn a = 4294967157)
{
	vector<tn> ret(n);
	ret[0] = X0;
	for(tn i = 1; i < n; i++)
		ret[i] = (a * ret[i-1] % m + c) % m;
	return ret;
}

// METHOD2
// 32)
// const tn n, tn X0 = 42949672,
// const tn m = 4294967296, const tn c = 4294967291, const tn d = 4294967156, const tn a = 4294967157
vector<tn> method2(
	const tn n, const tn X0 = 42949672,
	const tn m = 4294967291, const tn c = 4294967279, const tn d = 4294967231, const tn a = 4294967197)
{
	vector<tn> ret(n);
	ret[0] = X0;
	for (tn i = 1; i < n; i++)
		ret[i] = ((d * (ret[i - 1] * ret[i - 1] % m) % m + a * X0 % m) % m + c) % m;
	return ret;
}

vector<tn> method3(
	const tn n, const tn X0 = 1247437, const tn X1 = 224743647, const tn m = 4294967291)
{
	vector<tn> ret(n);
	ret[0] = X0;
	ret[1] = X1;
	for (tn i = 2; i < n; i++)
		ret[i] = (ret[i - 2] + ret[i - 1]) % m;
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
vector<tn> method4(
	const tn n, const tn X0 = 42949672,
	const tn m = 4294967197, const tn c = 4294967291, const tn a = 4294967157)
{
	vector<tn> ret(n);
	ret[0] = X0;
	for (tn i = 1; i < n; i++)
		ret[i] = (a * inverse(ret[i-1], m) % m + c) % m;
	return ret;
}

vector<tn> method5(const tn n, const tn m = 4294967291)
{
	vector<tn> ret(n);
	vector<tn> X = default_method(2, n);
	vector<tn> Y = default_method(4, n);
	for (tn i = 0; i < n; i++)
		ret[i] = (X[i] + (Y[i] > X[i] ? m : 0) - Y[i]) % m;
	return ret;
}

vector<double> method6(vector<double> rnd, const double m = 0, const double sigma = 1)
{
	const int D = 12, d = 6;
	const tn n = rnd.size();
	vector<double> ret(n / 12);
	for(tn i = 1; i <= n; i++)
	{
		ret[(i - 1) / D] += rnd[i - 1];
		if (i % D == 0)
			ret[(i - 1) / D] = m + (ret[(i - 1) / D] - d) * sigma;
	}
	return ret;
}

vector<tn> default_method(const int i, const tn n)
{
	switch (i)
	{
	case 1:
		return method1(n);
	case 2:
		return method2(n);
	case 3:
		return method3(n);
	case 4:
		return method4(n);
	default:
		return vector<tn>();
	}
}

int main()
{
	vector<double> r0 = U(method5(7200), 4294967291);
	vector<double> r = method6(r0);

	printList(r);
	cout << "\na\n";
	printHistogram(r, -3, 3, 11);
	cout << "\na\n";
	printHistogram(r0, 0, 1);
	cout << "\na\n";

	//printList(method4(100));
	//printHistogram(method4(10000), 0, 1);
}