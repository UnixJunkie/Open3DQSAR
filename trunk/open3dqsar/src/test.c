#include <stdio.h>


int main(int argc, char **argv)
{
  long one = 1;
  
  printf ("%d\n", (int)(!(*((char *)(&one)))));
  return 0;
}
