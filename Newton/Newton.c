#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>

FILE* out;

static void put(const char* format)
{
	fputs(format, stdout);
	fputs(format, out);
}

static void print(char* format, ...)
{
	int a = 0;
	double d = 0;
	char* c = format;
	va_list args;
	va_start(args, format);
	for (; *c; c++)
	{
		if (*c != '%')
		{
			printf("%c", *c);
			fprintf(out, "%c", *c);
			continue;
		}
		switch (*++c)
		{
			case 'd':
				a = va_arg(args, int);
				printf("%d", a);
				fprintf(out, "%d", a);
				break;
			case 'l':
				c++;
				d = va_arg(args, double);
				printf("%-7lg", d);
				fprintf(out, "%-7lg", d);
				break;
			case 'w':
				c++;
				d = va_arg(args, double);
				printf("%-15lg", d);
				fprintf(out, "%-15lg", d);
				break;
			default:
				printf("%c", *c);
				fprintf(out, "%c", *c);
		}
	}
	va_end(args);
}

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
		dif[i] = (double*)malloc((n - i) * sizeof(double));
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
			if (j + 1 < n)
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
	double res = c[0];
	for (int i = 1; i < n; i++)
		res += sub((xi - x[0]) / (x[1] - x[0]), i) * c[i];
	if (fabs(res) < 1e-13)	//avoids machine epsilon
		res = 0;
	return res;
}

int main(int argc, char const* argv[])
{
	FILE* in = stdin;
	int n = 0;
	if (argc > 1 && fopen_s(&in, argv[1], "r"))
	{
		printf("Cannot open file.\n");
		return 1;
	}
	if (fopen_s(&out, "output.txt", "w"))
	{
		printf("Cannot open file.\n");
		fclose(in);
		return 1;
	}
	if (argc > 2)
	{
		printf("Too many input arguments.\n");
		fclose(in);
		fclose(out);
		return 0;
	}
	while (1)
	{
		printf("Number of data points: ");
		fscanf_s(stdin, "%d", &n);
		if (n <= 1)
		{
			while (getchar() != '\n'); //skip junk
			printf("Wrong number.\n");
		}
		else
			break;
	}
	double* x = (double*)malloc(n * sizeof(double));
	double* y = (double*)malloc(n * sizeof(double));
	double* c = (double*)malloc(n * sizeof(double));
	if (!x || !y || !c)
	{
		put("Memory allocation error.\n");
		free(x);
		free(y);
		free(c);
		fclose(in);
		fclose(out);
		return 1;
	}
	system("cls");
	put("Input X values:\n");
	for (int i = 0; i < n; i++)
		switch (fscanf_s(in, "%lg", &x[i]))
		{
			case 0:
				fgetc(in);
				i--;
				break;
			case EOF:
				put("Not enough values.\n");
				goto end;
			default:
				print("%lg ", x[i]);
		}
	put("\nInput Y values:\n");
	for (int i = 0; i < n; i++)
		switch (fscanf_s(in, "%lg", &y[i]))
		{
			case 0:
				fgetc(in);
				i--;
				break;
			case EOF:
				put("Not enough values.\n");
				goto end;
			default:
				print("%lg ", y[i]);
		}
	print("\n\nNumber of data points = %d\n\nNewton polynomial coefficients:\n", n);
	if (newton(n, x, y, c))
		for (int i = 0; i < n; i++)
			print("c[%d] = %lg\n", i, c[i]);
	else
	{
		put("Newton function error.\n");
		free(x);
		free(y);
		free(c);
		fclose(in);
		fclose(out);
		return 1;
	}
	put("\nInterpolation:\nX               Y\n");
	for (double xi = x[0], dx = 0.1; xi <= x[n - 1] + dx / 10.0; xi += dx) //dx/10.0 avoids machine epsilon
		print("%wg %wg\n", xi, interpolate(n, x, c, xi));
end:
	printf("\nLog saved in \"output.txt\".\n");
	free(x);
	free(y);
	free(c);
	fclose(in);
	fclose(out);
	return 0;
}