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
void print(char* format, ...)
{
  int a = 0;
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
double sub(double q, int n)//substitute q(q-1)...(q-n+1)
{
  int i = 1;
  double res = q;
  for(; i < n; i++)
    res *= (q - i);
  return res;
}
double fact(int n)//factorial
{
  int i = 2;
  double d = 1.0;
  for(; i <= n; i++)
    d *= i;
  return d;
}
int newton(int n, double *x, double *y, double *c)
{
  int i = 0, j = 0;
  double **dif = NULL, dx = 0.1, xi = x[0], q = 0;
  dif = (double**)malloc(n * sizeof(double*));
  if(!dif)
  {
    put("Memory allocation error.\n");
    return 0;
  }
  for(i = 0; i < n; i++)
  {
    dif[i] = NULL;
    dif[i] = (double*)malloc((n - i) * sizeof(double));
    if(!dif[i])
    {
      put("Memory allocation error.\n");
      for(j = 0; j < i; j++)
        free(dif[i]);
      free(dif);
      return 0;
    }
  }
  for(i = 0; i < n; i++)//forward finite differences table
    dif[i][0] = y[i];
  for(i = 1; i < n; i++)
    for(j = 0; j < n - i; j++)
      dif[j][i] = dif[j + 1][i - 1] - dif[j][i - 1];
  /*c[0] = dif[0][0];
  printf("%lf\n", c[0]);
  xi = (x[0] + x[1]) / 2.0;
  q = (xi - x[0]) / (x[1] - x[0]);
  for(i = 1; i < n; i++)
  {
    c[i] = (sub(q, i) * dif[0][i]) / fact(i);
    c[i] /= pow(xi, i);
    printf("%lf\n", c[i]);
  }*/
  put("Interpolation:\nX               Y\n");
  for(j = 0; xi <= x[n - 1] + dx / 10.0; xi += dx, j++)//dx / 10.0 avoids machine epsilon
  {
    c[j] = dif[0][0];
    q = (xi - x[0]) / (x[1] - x[0]);
    for(i = 1; i < n; i++)
      c[j] += (sub(q, i) * dif[0][i]) / fact(i);
    if(fabs(c[j]) < 1e-13)
      c[j] = 0;
    print("%wg %wg\n", xi, c[j]);
  }
  for(i = 0; i < n; i++)
    free(dif[i]);
  free(dif);
  return 1;
}
int main(int argc, char const* argv[])
{
  FILE *in = stdin;
  int n = 0, i = 0;
  double *x = NULL, *y = NULL, *c = NULL;
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
    printf("Number of data points: ");
    fscanf_s(stdin, "%d", &n);
    if(n <= 1)
    {
      while(getchar() != '\n');//skip trash
      printf("Wrong number.\n");//hotline, baby
    }
    else
      break;
  }
  x = (double*)malloc(n * sizeof(double));
  y = (double*)malloc(n * sizeof(double));
  c = (double*)malloc(n * 10 * sizeof(double));
  if(!x || !y || !c)
  {
    put("Memory allocation error.\n");
    free(x);
    free(y);
    free(c);
    fclose(in);
    fclose(out);
    return 1;
  }
  printf("\n");
  put("Input X values:\n");
  for(i = 0; i < n; i++)
    switch(fscanf_s(in, "%lg", &x[i]))
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
  for(i = 0; i < n; i++)
    switch(fscanf_s(in, "%lg", &y[i]))
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
  print("\n\nNumber of data points == %d.\n", n);
  newton(n, x, y, c);
end:
  printf("Log saved in \"output.txt\".\n");
  free(x);
  free(y);
  free(c);
  fclose(in);
  fclose(out);
	return 0;
}