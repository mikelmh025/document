/*
*   standardInput.c
*   Reads three integers from stdin and prints them to stdout
*/
#include<stdio.h>
#include<stdlib.h>
int main(){
   int n, i; int x[3];
   printf("Enter three integers separated by spaces, then press return: ");
   n = scanf(" %d %d %d", &x[0], &x[1], &x[2]);
   printf("%d numbers were successfully read: ", n);
   for(i=0; i<n; i++) printf("%d ", x[i]);
   printf("\n");
   return EXIT_SUCCESS;
}
