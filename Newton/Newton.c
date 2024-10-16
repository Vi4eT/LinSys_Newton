#include "Misc.h"

static double sub(double q, int n)	//substitute q(q-1)...(q-n+1)
{
	double res = q;
	for (int i = 1; i < n; i++)
		res *= (q - i);
	return res;
}

static double fact(int n)
{
	double d = 1;
	for (int i = 2; i <= n; i++)
		d *= i;
	return d;
}

static int newton(int n, double* x, double* y, double* c)
{
	if (n < 2)
		return 0;
	double** dif = (double**)malloc(n * sizeof(double*));
	if (!dif)
		return 0;
	for (int i = 0; i < n; i++)
	{
		dif[i] = (double*)malloc(n * sizeof(double));			//C6385 if n - i
		if (!dif[i])
		{
			for (int j = 0; j < i; j++)
				free(dif[j]);
			free(dif);
			return 0;
		}
	}
	for (int i = 0; i < n; i++)
		dif[i][0] = y[i];
	for (int i = 1; i < n; i++)
		for (int j = 0; j < n - i; j++)
			if (j + 1 < n)										//C6385 if absent
				dif[j][i] = dif[j + 1][i - 1] - dif[j][i - 1];	//forward finite differences table
	for (int i = 0; i < n; i++)
		c[i] = dif[0][i] / (fact(i) * pow(x[1] - x[0], i));
	for (int i = 0; i < n; i++)
		free(dif[i]);
	free(dif);
	return 1;
}

static double interpolate(int n, double* x, double* c, double xi)
{
	double res = c[0], q = (xi - x[0]) / (x[1] - x[0]);
	for (int i = 1; i < n; i++)
		res += sub(q, i) * c[i];
	if (fabs(res) < 1e-13)	//avoids machine epsilon
		res = 0;
	return res;
}

int main(int argc, char const* argv[])
{
	FILE* in = stdin;
	if (setupfiles(argc, argv, &in))
		return 1;
	int n = 0;
	scan("Number of data points: ", "%d", &n, valn);
	double* x = (double*)malloc(n * sizeof(double));
	double* y = (double*)malloc(n * sizeof(double));
	double* c = (double*)malloc(n * sizeof(double));
	if (!x || !y || !c)
		return stdexit("Memory allocation error.\n", x, y, c, in);

	system("cls");
	put("Input X values:\n");
	if (!getarr(in, x, n))
		return stdexit("", x, y, c, in);
	put("\nInput Y values:\n");
	if (!getarr(in, y, n))
		return stdexit("", x, y, c, in);

	print("\n\nNumber of data points = %d\n\nNewton polynomial coefficients:\n", n);
	if (newton(n, x, y, c))
		for (int i = 0; i < n; i++)
			print("c[%d] = %lg\n", i, c[i]);
	else
		return stdexit("Newton function error.\n", x, y, c, in);

	put("\nInterpolation:\nX               Y\n");
	for (double xi = x[0], dx = 0.1; xi <= x[n - 1] + dx / 10.0; xi += dx) //dx/10.0 avoids machine epsilon
		print("%wg %wg\n", xi, interpolate(n, x, c, xi));

	printf("\nLog saved in \"output.txt\".\n");
	return stdexit("", x, y, c, in);
}