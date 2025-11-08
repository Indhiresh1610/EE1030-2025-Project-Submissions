#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define STB_IMAGE_IMPLEMENTATION
#include"../c_libs/stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../c_libs/stb_image_write.h"
double* imageload(const char* filename, int* m, int* n) {
    int w, h, c;
    unsigned char *data = stbi_load(filename, &w, &h, &c, 1);
    
    if (data == NULL) {
        printf("Error in loading image '%s'.\n", filename);
        return NULL;
    }
    double* matrix = (double*)malloc(w * h * sizeof(double));
    if (matrix == NULL) {
        printf("Error: Failed to allocate memory for matrix.\n");
        stbi_image_free(data);
        return NULL;
    }

    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {

            unsigned char pixelvalue = data[i * w + j];
            
            matrix[i * w + j] = (double)pixelvalue;
        }
    }

    stbi_image_free(data);
    *m = h;
    *n = w;
    return matrix;
}

void pngimage(const char* filename, int m, int n, double* matrix) {
    
    unsigned char* outpixels = (unsigned char*)malloc(m * n * sizeof(unsigned char));
    if (outpixels == NULL) {
        printf("Error in  allocating memory\n");
        return;
    }
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double pixelvalue = matrix[i * n + j];
            if (pixelvalue < 0.0){
                 pixelvalue = 0.0;
            }
            if (pixelvalue > 255.0){
                 pixelvalue = 255.0;
            }
            
            outpixels[i * n + j] = (unsigned char)round(pixelvalue);
        }
    }
    int result = stbi_write_png(filename, n, m, 1, outpixels,n);
    
    if (result == 0) {
        printf("Error in creating image '%s'.\n", filename);
    }

    free(outpixels);
}
void jpgimage(const char* filename, int m, int n, double* matrix) {

    unsigned char* outpixels = (unsigned char*)malloc(m * n * sizeof(unsigned char));
    if (outpixels == NULL) {
        printf("Error in  allocating memory\n");
        return;
    }


    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double pixelvalue = matrix[i * n + j];
            
            if (pixelvalue < 0.0){
             pixelvalue = 0.0;
        }
            if (pixelvalue > 255.0){

             pixelvalue = 255.0;
            }
            
            outpixels[i * n + j] = (unsigned char)round(pixelvalue);
        }
    }

    int quality = 30; 

    int result = stbi_write_jpg(filename, n, m, 1, outpixels, quality);
    
    if (result == 0) {
        printf("Error: Failed to write JPG image to '%s'.\n", filename);
    }

    free(outpixels);
}
double norm(int n,double v[n]){
    double sum=0.0;
    for(int i=0;i<n;i++){
        sum+=v[i]*v[i];
    }
    double a=sqrt(sum);
    return a;
}
void unit(int n,double v[n]){
    double mag=norm(n,v);
    if (mag < 1e-9) {
        return; 
    }
    for(int i=0;i<n;i++){
        v[i]=v[i]/mag;
    }
}
void matvec(int m,int n,double* a,double v[n],double b[m]){
    for(int i=0;i<m;i++){
        double sum=0.0;
        for(int j=0;j<n;j++){
            sum+=a[i*n +j]*v[j];
        }
        b[i]=sum;
    }
}
void trans(int m,int n, double* a,double* aT){
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            aT[j*m+i]=a[i*n+j];
        }
    }
}
void converge(int m, int n, double* a, double* aT, double v[n],double* v_old) {
    
    unit(n, v);       
    if (v_old == NULL) {
        printf("Error: Failed to allocate memory in converge.\n");
        return;
    }

    double tolerance = 1e-4; 
    int max_iterations = 25;

    for (int i = 0; i < max_iterations; i++) {
        
        for(int j=0; j<n; j++){
            v_old[j] = v[j];
        }


        double b[m];
        double c[n];
        matvec(m, n, a, v, b);
        matvec(n, m, aT, b, c);


        for (int j = 0; j < n; j++) {
            v[j] = c[j];
        }
        unit(n, v); 

        
        double diff_norm = 0.0;
        double flip_diff_norm = 0.0;
        for (int j = 0; j < n; j++) {
            diff_norm += (v[j] - v_old[j]) * (v[j] - v_old[j]);
            flip_diff_norm += (v[j] + v_old[j]) * (v[j] + v_old[j]);
        }
        
        if (diff_norm < tolerance * tolerance || flip_diff_norm < tolerance * tolerance) {
            
            break;
        }


    }

}
void rankonecomponent(int m, int n, double sigma, double u[m], double v[n], double* k) {
    
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            k[i*n +j] = sigma * u[i] * v[j];
        }
    }
}
void matrixsubtract(int m, int n, double *a, double* b) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            a[i*n +j] = a[i*n +j] - b[i*n+ j];
        }
    }
}
double frobeniusnorm(int m, int n, double* matrix) {
    double sum = 0.0;
    
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            sum += matrix[i*n+j] * matrix[i*n +j];
        }
    }
    
    return sqrt(sum);
}


int main(){

int m,n;
int k=20;
double* a = imageload("../../figs/globe/globe.jpg", &m, &n);
    if (a == NULL) {
        printf("Failed to load image.\n");
        return 1;
    }
    int max_k = (m < n) ? m : n;
    if (k > max_k) {
        printf("Error: k (%d) cannot be larger than the smallest image dimension (%d).\n", k, max_k);
        printf("Setting k to %d.\n", max_k);
        k = max_k;
    }
double sigma[k];
double* a2 = (double*)malloc(m * n * sizeof(double));
double* a3 = (double*)malloc(m * n * sizeof(double));
double* t = (double*)malloc(m * n * sizeof(double));
double* aT = (double*)malloc(n * m * sizeof(double));
double* v_old = (double*)malloc(n * sizeof(double));
double* v1 = (double*)malloc(n * k * sizeof(double));
double* u = (double*)malloc(m * k * sizeof(double));
    
if (a == NULL || a2 == NULL || a3 == NULL || t == NULL || aT == NULL || v1 == NULL || u == NULL) {
        printf("Error: Failed to allocate memory.\n");
        return 1;
    }
    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            a2[i * n + j] = a[i * n + j];
            t[i * n + j] = 0.0;
        }
    }
for(int i=0;i<k;i++){
    double v[n];
   for (int j = 0; j < n; j++) {
        v[j] = (double)rand() / (double)RAND_MAX;
    }
    trans(m, n, a, aT);
    converge(m,n,a,aT,v,v_old);

    for(int j=0;j<n;j++){
        v1[j*k +i]=v[j];
    }
    double temp[m];
    matvec(m,n,a,v,temp);
    sigma[i]=norm(m,temp);
    if (sigma[i] > 1e-9){
    for(int j=0;j<m;j++){
       temp[j]=temp[j]/sigma[i];
    }
    }
    for(int j=0;j<m;j++){
       u[j*k+i]=temp[j];
    }
   double* a1 = (double*)malloc(m * n * sizeof(double));
    rankonecomponent(m,n,sigma[i],temp,v,a1);
    matrixsubtract(m,n,a,a1);
    for(int p=0;p<m;p++){
        for(int q=0;q<n;q++){
            t[p*n +q]+=a1[p*n +q];
        }
    }
    free(a1);
}
    
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            a3[i*n +j]=a2[i*n+j]-t[i*n +j];
        }
    }
    double errornorm=frobeniusnorm(m,n,a3);
    double originalnorm=frobeniusnorm(m,n,a2);
    double percentageerror = 0.0;
    if (originalnorm > 1e-9) {
        percentageerror = (errornorm / originalnorm) * 100.0;
    }

    printf("absolutr error :%lf\n",errornorm);
    printf("percentage error:%lf",percentageerror);
    pngimage("../../figs/globe/globe_k20.png",m,n,t);
    jpgimage("../../figs/globe/globe_k20.jpg",m,n,t);
    free(a);
    free(a2);
    free(v_old);
    free(a3);
    free(t);
    free(aT);
    free(v1);
    free(u);
    return 0;
}
