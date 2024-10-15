#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
FILE* out;
void put(const char* format)
{
	fputs(format, stdout);
	fputs(format, out);
}
void print(char* format, ...)
{
	int a = 0, e = 0;
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
			case 'i':
				a = va_arg(args, int);
				printf("%2d", a);
				fprintf(out, "%2d", a);
				break;
			case 'l':
				c++;
				d = va_arg(args, double);
				printf("%lg", d);
				fprintf(out, "%lg", d);
				break;
			case 'v':
				c++;
				d = va_arg(args, double);
				printf("%-9lg", d);
				fprintf(out, "%-9lg", d);
				break;
			case 'w':
				c++;
				e = va_arg(args, int);
				d = va_arg(args, double);
				printf("%-15.*lg", e, d);
				fprintf(out, "%-15.*lg", e, d);
				break;
			default:
				printf("%c", *c);
				fprintf(out, "%c", *c);
		}
	}
	va_end(args);
}

static int count(double d)
{
	int count = 0;
	double threshold = 5e-14;
	while (fabs(d - round(d)) > threshold)
	{
		d *= 10.0;
		threshold *= 10.0;
		count++;
	}
	return count;
}

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

static int check(int n, double e, const double* a, const double* b, double* x)
{
	put("\nCheck:\n");
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

static int slau(int n, double e, const double* a, const double* b, double* x)
{
	if (!dominance(a, n))
	{
		put("Matrix is not diagonally dominant.\n");
		return 0;
	}
	double* A = (double*)malloc(n * n * sizeof(double));
	double* B = (double*)malloc(n * sizeof(double));
	double* x0 = (double*)malloc(n * sizeof(double));
	double* acc = (double*)malloc(n * sizeof(double)); //del
	if (!A || !B || !x0 || !acc)
	{
		put("Memory allocation error.\n");
		free(A);
		free(B);
		free(x0);
		free(acc);
		return 0;
	}
	double tmp = 0;
	for (int i = 0; i < n; i++)	//iterative matrix
	{
		tmp = a[i * n + i];
		for (int j = 0; j < n; j++)
			if (i == j)
				A[i * n + j] = 0;
			else
				A[i * n + j] = -a[i * n + j] / tmp;
		B[i] = x0[i] = b[i] / tmp;
	}
	put("Iterations:\n");
	for (int iter = 1, d = count(e); ; iter++) //d???
	{
		print("%i: ", iter);
		for (int i = 0; i < n; i++)
		{
			x[i] = B[i];
			for (int j = 0; j < n; j++)
				x[i] += A[i * n + j] * x0[j];
			print("%wg ", d, x[i]);
		}
		put("\n");
		tmp = 0;
		for (int i = 0; i < n; i++)
		{
			acc[i] = fabs(x[i] - x0[i]);
			if (tmp < acc[i])
				tmp = acc[i];
		}
		if (tmp < e)
			break;
		for (int i = 0; i < n; i++)
			x0[i] = x[i];
	}
	free(A);
	free(B);
	free(x0);
	free(acc);
	if (!check(n, e, a, b, x))
	{
		put("Solution check failed!\n");
		return 0;
	}
	else
	{
		put("Solution checked.\n");
		return 1;
	}
}

int main(int argc, char const* argv[])
{
	FILE* in = stdin;
	int n = 0;
	double e = 0;
	if (argc > 1 && fopen_s(&in, argv[1], "r"))
	{
		printf("Cannot open input file.\n");
		return 1;
	}
	if (fopen_s(&out, "output.txt", "w"))
	{
		printf("Cannot open output file.\n");
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
		printf("Matrix size: ");
		fscanf_s(stdin, "%d", &n);
		if (n <= 1)
		{
			while (getchar() != '\n');		//skip junk
			printf("Wrong size.\n");
		}
		else
			break;
	}
	while (getchar() != '\n');				//skip junk
	while (1)
	{
		printf("Accuracy: ");
		fscanf_s(stdin, "%lg", &e);
		if (e >= 1 || e <= 1e-13)			//avoid machine zero
		{
			while (getchar() != '\n');		//skip junk
			printf("Wrong accuracy.\n");
		}
		else
			break;
	}
	double* a = (double*)malloc(n * n * sizeof(double));
	double* b = (double*)malloc(n * sizeof(double));
	double* x = (double*)malloc(n * sizeof(double));
	if (!a || !b || !x)
	{
		put("Memory allocation error.\n");
		free(a);
		free(b);
		free(x);
		fclose(in);
		fclose(out);
		return 1;
	}
	printf("\n");
	put("Input matrix:\n");
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			switch (fscanf_s(in, "%lg", &a[i * n + j]))
			{
				case 0:
					fgetc(in);
					j--;
					break;
				case EOF:
					put("Not enough coefficients.\n");
					goto end;
				default:
					print("%vg ", a[i * n + j]);
			}
		put("\n");
	}
	put("\nInput vector:\n");
	for (int i = 0; i < n; i++)
		switch (fscanf_s(in, "%lg", &b[i]))
		{
			case 0:
				fgetc(in);
				i--;
				break;
			case EOF:
				put("Not enough coefficients.\n");
				goto end;
			default:
				print("%vg ", b[i]);
		}
	print("\n\nMatrix size == %dx%d, accuracy == %lg.\n", n, n, e);
	slau(n, e, a, b, x);
end:
	printf("Log saved in \"output.txt\".\n");
	free(a);
	free(b);
	free(x);
	fclose(in);
	fclose(out);
	return 0;
}