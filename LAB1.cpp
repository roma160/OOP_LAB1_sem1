#include <iomanip>
#include <iostream>
#include <vector>
#include <map>
#include <stdio.h>

using namespace std;

void printHistogram(const vector<double>& numbers, double start, double stopExcluding, int subdivisionsNum = 10)
{
	constexpr int precision = 2;
	const int width = max(log10(stopExcluding), start < 0 ? log10(-start) : 0) + precision + 2;

	cout.setf(ios::fixed, ios::floatfield);
	cout.precision(precision);
	
	cout << setw(2 * width + 5) << "Interval    |" << "  Frequency \n";

	int buff;
	double s, e;
	for(int i = 0; i < subdivisionsNum; i++)
	{
		buff = 0;
		s = (i * stopExcluding + (subdivisionsNum - i) * start) / subdivisionsNum;
		e = ((i + 1) * stopExcluding + (subdivisionsNum - i - 1) * start) / subdivisionsNum;

		cout << "[" << setw(width) << s << ";" << setw(width) << e << ")     ";

		for (int j = 0; j < numbers.size(); j++)
			if (numbers[j] >= s && numbers[j] < e)
				buff++;
		cout << (double)buff / subdivisionsNum << "\n";
	}
}

int main()
{
	printHistogram(vector<double>{0, 1, 1, 1, 2, 5, 5}, 0, 5, 10);
}