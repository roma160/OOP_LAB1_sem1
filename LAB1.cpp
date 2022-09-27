#include <iomanip>
#include <iostream>
#include <vector>
#include <map>
#include <stdio.h>
#include <cmath>
#include <math.h>

#define _USE_MATH_DEFINES


using namespace std;

typedef long long tn;


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
	virtual ~def_mod_args() = default;
	tn X0, m;
	def_mod_args(const tn X0, const tn m) : X0(X0), m(m) {}

	virtual def_mod_args* copy() const
	{ return new def_mod_args(*this); }
};
vector<tn> invoke_mod_method(int i, tn n, const def_mod_args* args = nullptr);

vector<double> U(const vector<tn> S, const def_mod_args* args)
{
	tn n = S.size();
	vector<double> ret(n);
	for (tn i = 0; i < n; i++)
		ret[i] = (double)S[i] / (args->m-1);
	return ret;
}

struct method1_args : def_mod_args
{
	const tn c, a;
	method1_args(const tn X0, const tn m, const tn c, const tn a) :
		def_mod_args(X0, m), c(c), a(a) {}

	def_mod_args* copy() const override
	{ return new method1_args(*this); }
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

	def_mod_args* copy() const override
	{ return new method2_args(*this); }
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

	def_mod_args* copy() const override
	{ return new method3_args(*this); }
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

	def_mod_args* copy() const override
	{ return new method4_args(*this); }
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
struct def_rnd_args
{
	const int mod_i;
	const def_mod_args* mod_args;

	def_rnd_args(const int mod_i, const def_mod_args* mod_args) :
		mod_i(mod_i), mod_args(mod_args) {}
};

struct method6_args: def_rnd_args
{
	const double m, sigma;

	method6_args(const int mod_i, const def_mod_args* mod_args, const double m, const double sigma):
		def_rnd_args(mod_i, mod_args), m(m), sigma(sigma) {}
};
const method6_args m6_p(4, &m4_p, 0, 1);
vector<double> method6(const tn n, const def_rnd_args* dargs = &m6_p)
{
	const method6_args* args = (const method6_args*)dargs;
	const int D = 12, d = 6;
	const tn dn = D * n;
	vector<double> rnd = U(
		invoke_mod_method(args->mod_i, dn, args->mod_args),
		args->mod_args
	);

	vector<double> ret(n);
	for(tn i = 1; i <= dn; i++)
	{
		ret[(i - 1) / D] += rnd[i - 1];
		if (i % D == 0)
			ret[(i - 1) / D] = args->m + (ret[(i - 1) / D] - d) * args->sigma;
	}
	return ret;
}

double sign(double n)
{
	if (n < 0) return -1;
	if (n > 0) return 1;
	return 0;
}
const def_rnd_args m7_p(4, &m4_p);
vector<double> method7(const tn n, const def_rnd_args* dargs = &m7_p)
{
	const tn N = n / 2 + n % 2;
	const method6_args* args = (const method6_args*)dargs;
	vector<double> rnd = U(
		invoke_mod_method(args->mod_i, 2*N, args->mod_args),
		args->mod_args
	);

	double S;
	for(tn i = 0; i < N; i++)
	{
		double& a = rnd[2 * i], & b = rnd[2 * i + 1];
		do
		{
			a = 2 * a - sign(a);
			b = 2 * b - sign(b);
		} while ((S = a*a + b*b) >= 1);

		S = sqrt(-2 * log(S) / S);
		a *= S;
		b *= S;
	}
	return rnd;
}

struct method8_args : def_rnd_args
{
	const int mod2_i;
	const def_mod_args* mod2_args;

	method8_args(
		const int mod1_i, const def_mod_args* mod1_args,
		const int mod2_i, const def_mod_args* mod2_args) :
		def_rnd_args(mod1_i, mod1_args), mod2_i(mod2_i), mod2_args(mod2_args) {}
};
const method8_args m8_p(1, &m1_p, 4, &m4_p);
vector<double> method8(const tn n, const def_rnd_args* dargs = &m8_p)
{
	const double m = sqrt(8 / exp(1.0));
	const double e4_4 = 4 * exp(0.25);
	const double e4_135 = 4 * exp(-1.35);
	const tn buffer_n = 2 * n;

	method8_args args = *(method8_args*)dargs;
	def_mod_args* args1 = args.mod_args->copy(), * args2 = args.mod2_args->copy();
	args.mod_args = args1;
	args.mod2_args = args2;

	bool restart;
	vector<double> ret(n);
	vector<double> u(buffer_n), v(buffer_n);
	
	for (tn i = 0; i < n; i++)
	{
		double& x = ret[i];
		restart = true;
		for (tn j = buffer_n; restart; j++)
		{
			restart = false;
			if (j >= buffer_n)
			{
				args1->X0 = u[buffer_n - 1] * args1->m;
				args2->X0 = v[buffer_n - 1] * args2->m;

				u = U(
					invoke_mod_method(args.mod_i, buffer_n, args.mod_args),
					args.mod_args
				);
				v = U(
					invoke_mod_method(args.mod2_i, buffer_n, args.mod2_args),
					args.mod2_args
				);

				j = 0;
			}

			if (abs(u[j]) <= 1e-20)
			{
				restart = true;
				continue;
			}
			x = (v[j] - 0.5) * m / u[j];

			if (x * x <= 5 - e4_4 * u[j]) continue;
			if (x * x >= e4_135 / u[j] + 1.4 || x * x > -4 * log(u[j]))
			{
				restart = true;
				continue;
			}
		}
	}

	delete args1;
	delete args2;
	return ret;
}

struct method9_args: def_rnd_args
{
	const double mu;

	method9_args(const int mod_i, const def_mod_args* mod_args, const double mu) :
		def_rnd_args(mod_i, mod_args), mu(mu) {}
};
const method9_args m9_p(4, &m4_p, 1);
vector<double> method9(const tn n, const def_rnd_args* dargs = &m9_p)
{
	const method9_args* args = (const method9_args*)dargs;
	vector<double> rnd = U(
		invoke_mod_method(args->mod_i, n, args->mod_args),
		args->mod_args
	);

	const double m = -args->mu;

	vector<double> ret(n);
	for (tn i = 0; i < n; i++) ret[i] = m * log(rnd[i]);
	return ret;
}

struct method10_args : def_rnd_args
{
	const double a;

	method10_args(const int mod_i, const def_mod_args* mod_args, const double a) :
		def_rnd_args(mod_i, mod_args), a(a) {}
};
const method10_args m10_p(4, &m4_p, 1);
vector<double> method10(const tn n, const def_rnd_args* dargs = &m10_p)
{
	method10_args args = *(method10_args*)dargs;
	def_mod_args* args1 = args.mod_args->copy();
	args.mod_args = args1;

	const double pi = atan(1) * 4;
	const double& a = args.a;
	const double sq_a = sqrt(2 * a - 1);
	const tn buffer_n = 2 * n;

	bool restart;
	vector<double> ret(n);
	vector<double> rnd(buffer_n);

	for (tn i = 0; i < n; i++)
	{
		double& x = ret[i];
		restart = true;
		for (tn j = buffer_n; restart; j += 2)
		{
			restart = false;
			if (j >= buffer_n - 1)
			{
				args1->X0 = rnd[buffer_n - 1] * args1->m;
				rnd = U(
					invoke_mod_method(args.mod_i, buffer_n, args.mod_args),
					args.mod_args
				);

				j = 0;
			}

			double Y = tan(pi * rnd[j]);
			double X = sqrt(2 * a - 1) * Y + a - 1;
			if(X <= 0 || rnd[j+1] > exp((a-1)*log(X/(a-1)) - sq_a*Y) * (1 + Y * Y))
			{
				restart = true;
				continue;
			}
			ret[i] = X;
		}
	}
	return ret;
}

int main()
{
	vector<double> r = method10(1000);

	printList(r);

	//printList(method4(100));
	printHistogram(r, -2, 2);
}