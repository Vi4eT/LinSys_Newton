#include "Misc.h"

int setupfiles(int argc, char const* argv[], FILE** in)
{
	if (argc > 1 && fopen_s(in, argv[1], "r")) 
	{
		printf("Cannot open input file.\n");
		return 1;
	}
	if (fopen_s(&out, "output.txt", "w"))
	{
		printf("Cannot open output file.\n");
		fclose(*in);
		return 1;
	}
	if (argc > 2) 
	{
		printf("Too many input arguments.\n");
		fclose(*in);
		fclose(out);
		return 1;
	}
	return 0;
}

void put(const char* buffer)
{
	fputs(buffer, stdout);
	fputs(buffer, out);
}

void print(char* format, ...)
{
	int a = 0, precision = 0;
	double d = 0;
	va_list args;
	va_start(args, format);
	for (char* c = format; *c; c++)
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
				printf("%-7lg", d);
				fprintf(out, "%-7lg", d);
				break;
			case 'w':
				c++;
				d = va_arg(args, double);
				printf("%-15lg", d);
				fprintf(out, "%-15lg", d);
				break;
			case 'p':
				c++;
				precision = va_arg(args, int);
				d = va_arg(args, double);
				printf("%-20.*lg", precision, d);
				fprintf(out, "%-20.*lg", precision, d);
				break;
			default:
				printf("%c", *c);
				fprintf(out, "%c", *c);
		}
	}
	va_end(args);
}

int valn(void* n)
{
	return *(int*)n <= 2;
}

int valacc(void* e)
{
	return *(double*)e > 0.1 || *(double*)e < 1e-13;	//avoid machine zero
}

void scan(const char* prompt, const char* format, void* x, int(*val)(void*))
{
	while (1)
	{
		printf("%s", prompt);
		fscanf_s(stdin, format, x);
		if (val(x))
		{
			while (getchar() != '\n');	//skip junk
			printf("Error!\n");
		}
		else
			break;
	}
	while (getchar() != '\n');
}

int getarr(FILE* in, double* arr, int n)
{
	for (int i = 0; i < n; i++)
		switch (fscanf_s(in, "%lg", &arr[i]))
		{
			case 0:
				fgetc(in);
				i--;
				break;
			case EOF:
				put("\nNot enough values.\n");
				return 0;
			default:
				print("%lg ", arr[i]);
		}
	return 1;
}

int stdexit(const char* buffer, double* x, double* y, double* c, FILE* in)
{
	put(buffer);
	free(x);
	free(y);
	free(c);
	fclose(in);
	fclose(out);
	if (!*buffer)
		return 0;
	return 1;
}