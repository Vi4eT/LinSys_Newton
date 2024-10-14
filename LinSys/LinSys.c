#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <stdarg.h>
FILE *out;
void put(const char* format)
{
	fputs(format, stdout);
	fputs(format, out);
}
int count(double d)
{
	int count = 0;
	double threshold = 5e-14;
	while(fabs(d - round(d)) > threshold)
	{
		d *= 10.0;
		threshold *= 10.0;
		count++;
	}
	return count;
}
void print(char* format, ...)
{
	int a = 0, e = 0;
	double d = 0;
	char *c = format;
	va_list args;
	va_start(args, format);
	for(; *c; c++)
	{
		if(*c != '%')
		{
			printf("%c", *c);
			fprintf(out, "%c", *c);
			continue;
		}
		switch(*++c)
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
int dominance(const double *a, int n)
{
	int i = 0, j = 0;
	double sum = 0;
	for(i = 0; i < n; sum = 0, i++)
	{
		for(j = 0; j < n; j++)
			if(i != j)
				sum += fabs(a[i * n + j]);
		if(fabs(a[i * n + i]) <= sum)//strict diagonal dominance
			return 0;
	}
	return 1;
}
int check(int n, double e, const double *a, const double *b, double *x)
{
	int i = 0, j = 0;
	double line = 0;
	put("\nCheck:\n");
	for(i = 0; i < n; line = 0, i++)
	{
		for(j = 0; j < n; j++)
			line += a[i * n + j] * round(x[j] / e) * e;
		if(fabs(line - b[i]) > e)
		{
			print("%lg != %lg\n", line, b[i]);
			return 0;
		}
		print("%lg == %lg\n", line, b[i]);
	}
	return 1;
}
int slau(int n, double e, const double *a, const double *b, double *x)
{
	int i = 0, j = 0, k = 0, d = count(e);
	double *A = NULL, *B = NULL, *x0 = NULL, *acc = NULL, tmp = 0;
	if(dominance(a, n))
	{
		A = (double*)malloc(n * n * sizeof(double));
		B = (double*)malloc(n * sizeof(double));
		x0 = (double*)malloc(n * sizeof(double));
		acc = (double*)malloc(n * sizeof(double));
		if(!A || !B || !x0 || !acc)
		{
			put("Memory allocation error.\n");
			free(A);
			free(B);
			free(x0);
			free(acc);
			return 0;
		}
		for(i = 0; i < n; i++)//iterative matrix
		{
			tmp = a[i * n + i];
			for(j = 0; j < n; j++)
				if(i == j)
					A[i * n + j] = 0;
				else
					A[i * n + j] = -a[i * n + j] / tmp;
			B[i] = x0[i] = b[i] / tmp;
		}
		/*iterative matrix norm
		for(i = 0; i < n; i++)
		{
			sum[i] = 0;
			for(j = 0; j < n; j++)
			{
				sum[i] += fabs(A[i * n + j]);
				printf("%.4f ", A[i * n + j]);
			}
			if(max < sum[i])
				max = sum[i];
			printf("\n");
		}
		if(max >= 1)
		{
			printf("Iterative method is not convergent.\n");
			free(A);
			free(B);
			free(x0);
			return 0;
		}
		printf("\n");*/
		put("Iterations:\n");
		for(;;)//iterative method
		{
			print("%i: ", ++k);
			for(i = 0; i < n; i++)
			{
				x[i] = B[i];
				for(j = 0; j < n; j++)
					x[i] += A[i * n + j] * x0[j];
				print("%wg ", d, x[i]);
			}
			put("\n");
			tmp = 0;
			for(i = 0; i < n; i++)
			{
				acc[i] = fabs(x[i] - x0[i]);
				if(tmp < acc[i])
					tmp = acc[i];
			}
			if(tmp < e)
				break;
			for(i = 0; i < n; i++)
				x0[i] = x[i];
		}
		free(A);
		free(B);
		free(x0);
		free(acc);
		if(!check(n, e, a, b, x))
		{
			put("Solution check failed.\n");
			return 0;
		}
		else
		{
			put("Solution checked.\n");
			return 1;
		}
	}
	else
	{
		put("Matrix is not diagonally dominant.\n");
		return 0;
		/*row switching
		A = (double*)malloc(n * n * sizeof(double));
		queue = (int*)malloc(n * sizeof(int));
		if(A && queue)
		{
			for(i = 0; i < n; i++)
				queue[i] = -1;
			for(i = 0; i < n; i++)
			{
				tmp = 0;
				k = 0;
				for(j = 0; j < n; j++)
				{
					if(tmp <= fabs(a[i * n + j]))
					{
						tmp = fabs(a[i * n + j]);
						k = j;
					}
				}
				if(a[i * n + i] == tmp)
					queue[i] = i;
				else
					if(queue[k] == -1)
						queue[k] = i;
					else
						for(l = 0; l < n; l++)
							if(queue[l] == -1)
							{
								queue[l] = i;
								break;
							}
			}
			for(i = 0; i < n; i++)
			{
				for(j = 0; j < n; j++)
				{
					A[i * n + j] = a[queue[i] * n + j];
					printf("%.4f ", A[i * n + j]);
				}
				printf("\n");
			}
			slau(n, e, A, b, x);
			free(A);
			free(queue);
		}
		else
		{
			printf("Memory allocation error.\n");
			return 0;
		}*/
	}
}
int main(int argc, char const* argv[])
{
	FILE *in = stdin;
	int n = 0, i = 0, j = 0;
	double *a = NULL, *b = NULL, *x = NULL, e = 0;
	if(argc > 1 && fopen_s(&in, argv[1], "r"))
	{
		printf("Cannot open file.\n");
		return 1;
	}
	if(fopen_s(&out, "output.txt", "w"))
	{
		printf("Cannot open file.\n");
		fclose(in);
		return 1;
	}
	if(argc > 2)
	{
		printf("Too many input arguments.\n");
		fclose(in);
		fclose(out);
		return 0;
	}
	for(;;)
	{
		printf("Matrix size: ");
		fscanf_s(stdin, "%d", &n);
		if(n <= 1)
		{
			while(getchar() != '\n');//skip trash
			printf("Wrong size.\n");
		}
		else
			break;
	}
	while(getchar() != '\n');//skip trash
	for(;;)
	{
		printf("Accuracy: ");
		fscanf_s(stdin, "%lg", &e);
		if(e >= 1 || e <= 1e-13)//avoid machine zero
		{
			while(getchar() != '\n');//skip trash
			printf("Wrong accuracy.\n");
		}
		else
			break;
	}
	a = (double*)malloc(n * n * sizeof(double));
	b = (double*)malloc(n * sizeof(double));
	x = (double*)malloc(n * sizeof(double));
	if(!a || !b || !x)
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
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
			switch(fscanf_s(in, "%lg", &a[i * n + j]))
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
	for(i = 0; i < n; i++)
		switch(fscanf_s(in, "%lg", &b[i]))
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