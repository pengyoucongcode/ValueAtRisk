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

//计算投资组合期望收益率,参数依次为期望收益率数组、投资比例数组、股票支数
double expectedRateOfReturn(double num1[], double num2[],int len) 
{
	double result = 0.0;
	for (int i = 0; i < len; i++)
	{
		result += (num1[i] * num2[i]);
	}
	return result;
}

//计算协方差矩阵,参数依次为协方差矩阵、相关系数矩阵、标准差数组、维数
void covarianceMatrix(double** num, double** num1, double num2[], int len) 
{
	for (int i = 0; i < len; i++)
	{
		for (int j = 0; j < len; j++)
		{
			num[i][j] = num1[i][j] * num2[i] * num2[j];
		}
	}
}

// 计算组合标准差(分散标准差）,参数依次为标准差、投资比例、相关系数矩阵、维度
double StandardDeviation(double num1[], double num2[], double** num3,int len)
{
	double* temp1, *temp2;
	temp1 = new double[len];
	temp2 = new double[len];
	for (int i = 0; i < len; i++)
	{
		double t = num1[i] * num2[i];
		temp1[i] = pow(t, 2);
	}
	int k = 0;
	for (int i = 0; i < len-1; i++)
	{
		for (int j = i + 1; j < len; j++)
		{
			temp2[k++] = 2 * (num1[i] * num2[i]) * (num1[j] * num2[j]) * num3[i][j];
		}
	}
	double sum = 0.0;
	for (int i = 0; i < len; i++)
		sum += (temp1[i] + temp2[i]);
	double result = sqrt(sum);
	return result;
}

// 非分散标准差计算，参数依次为标准差数组、投资比例数组、长度
double nonDispersiveStandardDeviation(double num1[], double num2[],int len)
{
	double sum = 0.0;
	for (int i = 0; i < len; i++)
		sum += (num1[i] * num2[i]);
	return sum;

}

// 风险价值计算，参数依次为投资总额、Z值（标准正态分布的抽样分位数）、标准差、期望收益率
double VaR(double W, double Z, double R, double Q)
{
	double result = W*(Z*R+Q);
	return result;
}

// 非分散风险价值，参数依次为总资本、Z值（标准正态分布的抽样分位数）、投资组合非分散标准差
double VaR_Undx(double W, double Z, double R)
{
	double result = W * abs(Z) * R;
	return result;
}

// 分散风险价值计算，参数依次为总资本、Z值（标准正态分布的抽样分位数）、投资组合分散标准差
double VaR_dx(double W, double Z, double R)
{
	double result = W * abs(Z) * R;
	return result;
}
int main()
{
	int len; // 股票支 数
	double count; // 投资总资本
	double a; // 置信度
	double *rateOfReturn; // 期望收益率数组 
	double* standardDeviation; // 标准差数组
	double* proportion; // 投资比例数组
	double** correlationCoefficient; // 相关系数矩阵
	double** covariance; // 协方差矩阵
	cout << "请输入总资本：";
	cin >> count;
	cout << "请输入一共有多少支股票：";
	cin >> len;
	correlationCoefficient = new double* [len];
	covariance = new double* [len];
	for (int i = 0; i < len; i++)
	{
		correlationCoefficient[i] = new double[len];
		covariance[i] = new double[len];
	}
	rateOfReturn = new double[len];
	standardDeviation = new double[len];
	proportion = new double[len];
	cout << "请依次输入每支股票的收益率、标准差、投资比例：" << endl;
	for (int j = 0; j < len; j++)
		cin >> rateOfReturn[j] >> standardDeviation[j] >> proportion[j];
	cout << "请输入股票相关系数矩阵（规则：从左到右，从上到下）：" << endl;
	for (int x = 0; x < len; x++)
	{
		for (int y = 0; y < len; y++)
			cin >> correlationCoefficient[x][y];
	}
	double combinedRateOfReturn = expectedRateOfReturn(rateOfReturn, proportion, len);//组合期望收益率
	covarianceMatrix(covariance, correlationCoefficient, standardDeviation, len);
	cout << "该投资组合的期望收益率为：" << combinedRateOfReturn << endl;
	cout << "协方差矩阵：" << endl;
	for (int i = 0; i < len; i++)
	{
		for (int j = 0; j < len; j++) 
			cout << covariance[i][j] << " ";
		cout << endl;
	}
	double combinedStandardDeviation = StandardDeviation(standardDeviation,proportion, correlationCoefficient, len); // 分散标准差
	cout << "该组合的标准差(分散标准差）：" << combinedStandardDeviation << endl;
	double non_dispersiveStandardDeviation = nonDispersiveStandardDeviation(standardDeviation, proportion, len); // 非分散标准差
	cout << "该组合的非分散标准差：" << non_dispersiveStandardDeviation << endl;
	cout << "请输入置信度：";
	cin >>a;
	double Z = normsinv(1 - a);
	double combinedRiskValue = VaR(count, Z, combinedStandardDeviation, combinedRateOfReturn); //组合风险价值
	cout << "组合风险价值：" << combinedRiskValue << endl;
	double diversified_risk_value = VaR_dx(count, Z, combinedStandardDeviation);
	double non_distributed_risk_value = VaR_Undx(count, Z, non_dispersiveStandardDeviation);
	cout << "该组合的分散风险价值：" << diversified_risk_value << " 非分散风险价值：" << non_distributed_risk_value << endl;
	return 0;
}