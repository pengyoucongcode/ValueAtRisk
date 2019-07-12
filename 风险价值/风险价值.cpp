// 风险价值.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include "pch.h"
#include <iostream>
#include <math.h>
#include <errno.h>
using namespace std;

const double PI = 3.1415926;
double Normal(double z)
{
	double temp;
	temp = exp((-1) * z * z / 2) / sqrt(2 * PI);
	return temp;

}
/***************************************************************/
/* 返回标准正态分布的累积函数，该分布的平均值为 0，标准偏差为 1。                           */
/***************************************************************/
double NormSDist(const double z)
{
	// this guards against overflow
	if (z > 6) return 1;
	if (z < -6) return 0;

	static const double gamma = 0.231641900,
		a1 = 0.319381530,
		a2 = -0.356563782,
		a3 = 1.781477973,
		a4 = -1.821255978,
		a5 = 1.330274429;

	double k = 1.0 / (1 + fabs(z) * gamma);
	double n = k * (a1 + k * (a2 + k * (a3 + k * (a4 + k * a5))));
	n = 1 - Normal(z) * n;
	if (z < 0)
		return 1.0 - n;

	return n;
}


/***************************************************************/
/* 返回标准正态分布累积函数的逆函数。该分布的平均值为 0，标准偏差为 1。                    */
/***************************************************************/
double normsinv(const double p)
{
	static const double LOW = 0.02425;
	static const double HIGH = 0.97575;

	/* Coefficients in rational approximations. */
	static const double a[] =
	{
		-3.969683028665376e+01,
			2.209460984245205e+02,
			-2.759285104469687e+02,
			1.383577518672690e+02,
			-3.066479806614716e+01,
			2.506628277459239e+00
	};

	static const double b[] =
	{
		-5.447609879822406e+01,
			1.615858368580409e+02,
			-1.556989798598866e+02,
			6.680131188771972e+01,
			-1.328068155288572e+01
	};

	static const double c[] =
	{
		-7.784894002430293e-03,
			-3.223964580411365e-01,
			-2.400758277161838e+00,
			-2.549732539343734e+00,
			4.374664141464968e+00,
			2.938163982698783e+00
	};

	static const double d[] =
	{
		7.784695709041462e-03,
			3.224671290700398e-01,
			2.445134137142996e+00,
			3.754408661907416e+00
	};

	double q, r;

	errno = 0;

	if (p < 0 || p > 1)
	{
		errno = EDOM;
		return 0.0;
	}
	else if (p == 0)
	{
		errno = ERANGE;
		return -HUGE_VAL /* minus "infinity" */;
	}
	else if (p == 1)
	{
		errno = ERANGE;
		return HUGE_VAL /* "infinity" */;
	}
	else if (p < LOW)
	{
		/* Rational approximation for lower region */
		q = sqrt(-2 * log(p));
		return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
			((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
	}
	else if (p > HIGH)
	{
		/* Rational approximation for upper region */
		q = sqrt(-2 * log(1 - p));
		return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
			((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
	}
	else
	{
		/* Rational approximation for central region */
		q = p - 0.5;
		r = q * q;
		return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
			(((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
	}
}

//计算投资组合期望收益率
double expectedRateOfReturn(double num1[], double num2[],int len)
{
	double result = 0.0;
	for (int i = 0; i < len; i++)
	{
		result += (num1[i] * num2[i]);
	}
	return result;
}

int main()
{
	//double count, bzc, jz, a; // 依次为总资本、标准差、收益率均值，置信水平
	//cout << "请依次输出总资本（万元）、标准差、收益率均值，置信水平：";
	//cin >> count >> bzc >> jz >> a;
	//double Z = normsinv(1 - a);
	//cout << "Z值：" << Z << endl;
	//double VaR = count * (Z * bzc + jz);
	//cout << "风险价值为(万元）：" << VaR << endl;
}

