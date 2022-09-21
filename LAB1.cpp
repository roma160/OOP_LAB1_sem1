#include <iomanip>
#include <iostream>
#include <vector>
#include <map>
#include <stdio.h>

using namespace std;

typedef long long tn;
typedef vector<double>(*method_type)(tn n);


// PRINTING METHODS SECTION

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


// MOD METHODS SECTION

struct def_mod_args
{
	const tn X0, m;
	def_mod_args(const tn X0, const tn m) : X0(X0), m(m) {}
};
vector<tn> invoke_mod_method(int i, tn n, const def_mod_args* args = nullptr);

vector<double> U(const vector<tn> S, const tn m)
{
	tn n = S.size();
	vector<double> ret(n);
	for (tn i = 0; i < n; i++)
		ret[i] = (double)S[i] / (m-1);
	return ret;
}

struct method1_args : def_mod_args
{
	const tn c, a;
	method1_args(const tn X0, const tn m, const tn c, const tn a) :
		def_mod_args(X0, m), c(c), a(a) {}
};
const method1_args m1_31b{ 2147483, 2147483648, 2147483647, 2147483637 };
const method1_args m1_32b{ 42949672, 4294967296, 4294967291, 4294967157 };
const method1_args m1_p{ 42949672, 4294967291, 4294967279, 4294967231 };
vector<tn> method1(const tn n, const def_mod_args* dargs = &m1_p)
{
	const method1_args* args = (const method1_args*) dargs;
	vector<tn> ret(n);
	ret[0] = args->X0;
	for(tn i = 1; i < n; i++)
		ret[i] = (args->a * ret[i-1] % args->m + args->c) % args->m;
	return ret;
}

struct method2_args: def_mod_args
{
	const tn c, d, a;
	method2_args(const tn X0, const tn m, const tn c, const tn d, const tn a) :
		def_mod_args(X0, m), c(c), d(d), a(a) {}
};
const method2_args m2_32b{ 42949672, 4294967296, 4294967291, 4294967156, 4294967157 };
const method2_args m2_p{ 42949672, 4294967291, 4294967279, 4294967231, 4294967197 };
vector<tn> method2(const tn n, const def_mod_args* dargs = &m2_p)
{
	const method2_args* args = (const method2_args*) dargs;
	vector<tn> ret(n);
	ret[0] = args->X0;
	for (tn i = 1; i < n; i++)
		ret[i] = (
			(args->d * (ret[i - 1] * ret[i - 1] % args->m) % args->m + 
				args->a * args->X0 % args->m) % args->m + 
				args->c
			) % args->m;
	return ret;
}

struct method3_args : def_mod_args
{
	const tn X1;
	method3_args(const tn X0, const tn X1, const tn m) :
		def_mod_args(X0, m), X1(X1) {}
};
const method3_args m3_p(1247437, 224743647, 4294967291);
vector<tn> method3(const tn n, const def_mod_args* dargs = &m3_p)
{
	const method3_args* args = (const method3_args*) dargs;
	vector<tn> ret(n);
	ret[0] = args->X0;
	ret[1] = args->X1;
	for (tn i = 2; i < n; i++)
		ret[i] = (ret[i - 2] + ret[i - 1]) % args->m;
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
struct method4_args: def_mod_args
{
	const tn c, a;
	method4_args(const tn X0, const tn m, const tn c, const tn a) :
		def_mod_args(X0, m), c(c), a(a) {}
};
const method4_args m4_32b(42949672, 4294967296, 4294967290, 4294967157);
const method4_args m4_p(42949672, 4294967197, 4294967291, 4294967157);
vector<tn> method4(const tn n, const def_mod_args* dargs = &m4_p)
{
	const method4_args* args = (const method4_args*)dargs;
	vector<tn> ret(n);
	ret[0] = args->X0;
	for (tn i = 1; i < n; i++)
		ret[i] = (args->a * inverse(ret[i-1], args->m) % args->m + args->c) % args->m;
	return ret;
}

struct method5_args: def_mod_args
{
	const int method_x_i;
	const def_mod_args* x_args;

	const int method_y_i;
	const def_mod_args* y_args;

	method5_args(const tn m, const int method_x_i, const def_mod_args* x_args,
			const int method_y_i, const def_mod_args* y_args) :
		def_mod_args(-1, m), method_x_i(method_x_i), x_args(x_args), method_y_i(method_y_i), y_args(y_args) {}
};
const method5_args m5_p(4294967291, 2, nullptr, 4, nullptr);
vector<tn> method5(const tn n, const def_mod_args* dargs = &m5_p)
{
	const method5_args* args = (const method5_args*)dargs;
	vector<tn> ret(n);
	vector<tn> X = invoke_mod_method(args->method_x_i, n, args->x_args);
	vector<tn> Y = invoke_mod_method(args->method_y_i, n, args->y_args);
	for (tn i = 0; i < n; i++)
		ret[i] = (X[i] + (Y[i] > X[i] ? args->m : 0) - Y[i]) % args->m;
	return ret;
}

vector<tn> invoke_mod_method(int i, tn n, const def_mod_args* args)
{

	switch (i)
	{
	case 1:
		if (args == nullptr) return method1(n);
		return method1(n, args);
	case 2:
		if (args == nullptr) return method2(n);
		return method2(n, args);
	case 3:
		if (args == nullptr) return method3(n);
		return method3(n, args);
	case 4:
		if (args == nullptr) return method4(n);
		return method4(n, args);
	case 5:
		if (args == nullptr) return method5(n);
		return method5(n, args);
	default:
		return vector<tn>();
	}
}



// RND METHODS SECTION

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