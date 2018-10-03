int main () {
   /*This function is for plotting*/
   /* local variable definition */
   int b = 2;
   double a[] ={ 0.1, 2};
   double ret; 
   /* calling a function to get max value */
   int i, j; 
   FILE *fp;
   fp = fopen("numbers.txt", "w");

   for (j = 0; j<100; j++){
       for ( i = 0; i<100; i++){
           a[0] =  i*M_PI/100;
           a[1] = j*0.02;
           fprintf(fp,"%4.9f %.9f %2.4f\n ", a[0],a[1],f(b,a));
       }
   }
   fclose(fp);
   return 0;
}

