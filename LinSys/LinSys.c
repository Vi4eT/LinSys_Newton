#include "Misc.h"

static int dominance(const double* a, int n)
{
	double sum = 0;
	for (int i = 0; i < n; i++, sum = 0)
	{
		for (int j = 0; j < n; j++)
			if (i != j)
				sum += fabs(a[i * n + j]);
		if (fabs(a[i * n + i]) <= sum)		//strict diagonal dominance
			return 0;
	}
	return 1;
}

static int verify(int n, double e, const double* a, const double* b, double* x)
{
	double line = 0;
	for (int i = 0; i < n; i++, line = 0)
	{
		for (int j = 0; j < n; j++)
			line += a[i * n + j] * round(x[j] / e) * e;
		print("%lg == %lg\n", line, b[i]);
		if (fabs(line - b[i]) > e)
			return 0;
	}
	return 1;
}

static int slau(int n, double e, const double* a, const double* b, double* x)	//name is mandatory
{
	double* A = (double*)malloc(n * n * sizeof(double));
	double* B = (double*)malloc(n * sizeof(double));
	double* x0 = (double*)malloc(n * sizeof(double));
	if (!A || !B || !x0)
	{
		free(A);
		free(B);
		free(x0);
		return 0;
	}
	for (int i = 0; i < n; i++)		//iterative matrix
	{
		double tmp = a[i * n + i];
		for (int j = 0; j < n; j++)
			if (i == j)
				A[i * n + j] = 0;
			else
				A[i * n + j] = -a[i * n + j] / tmp;
		B[i] = x0[i] = b[i] / tmp;
	}
	for (int iter = 1; ; iter++)
	{
		print("%i: ", iter);
		double minacc = 0;
		for (int i = 0; i < n; i++)
		{
			x[i] = B[i];
			for (int j = 0; j < n; j++)
				x[i] += A[i * n + j] * x0[j];
			double acc = fabs(x[i] - x0[i]);
			if (minacc < acc)
				minacc = acc;
			print("%pg ", (int)ceil(fabs(log10(e))) + 1, x[i]); //precision, e > 0
			x0[i] = x[i];
		}
		put("\n");
		if (minacc < e)
			break;
	}
	free(A);
	free(B);
	free(x0);
	return 1;
}

int main(int argc, char const* argv[])
{
	FILE* in = stdin;
	if (setupfiles(argc, argv, &in))
		return 1;
	int n = 0;
	double e = 0;
	scan("Matrix row count: ", "%d", &n, valn);
	scan("Accuracy: ", "%lg", &e, valacc);
	double* a = (double*)malloc(n * n * sizeof(double));
	double* b = (double*)malloc(n * sizeof(double));
	double* x = (double*)malloc(n * sizeof(double));
	if (!a || !b || !x)
		return stdexit("Memory allocation error.\n", a, b, x, in);

	system("cls");
	put("Input matrix:\n");
	for (int i = 0; i < n; i++)
	{
		if (!getarr(in, &a[i * n], n))
			return stdexit("", a, b, x, in);
		put("\n");
	}
	put("\nInput vector:\n");
	if (!getarr(in, b, n))
		return stdexit("", a, b, x, in);
	print("\n\nMatrix size == %dx%d, accuracy == %lg\n", n, n, e);

	if (!dominance(a, n))
		return !stdexit("Matrix is not diagonally dominant.\n", a, b, x, in);
	put("\nIterations:\n");
	if (!slau(n, e, a, b, x))
		return stdexit("Memory allocation error.\n", a, b, x, in);

	put("\nVerification:\n");
	if (!verify(n, e, a, b, x))
		put("Verification failed!\n");
	else
		put("Solution verified.\n");

	printf("\nLog saved in \"output.txt\".\n");
	return stdexit("", a, b, x, in);
}