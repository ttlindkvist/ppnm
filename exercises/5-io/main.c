#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<assert.h>

void printline(double x){
        printf("| %g\t| %g\t| %g\t|\n", x, cos(x), sin(x));
}
void printheading(){
    printf("-------------------------------------\n");
    printf("| x\t| cos(x)\t| sin(x)\t|\n");
    printf("-------------------------------------\n");
}

void readcmd(int argc, char **argv){
    printf("\nREADING WITH CMD\n");
    printheading();
    for(int i = 2; i<argc; i++){
        double x = atof(argv[i]);
        printline(x);
    }
    printf("-------------------------------------\n");
}
void readstdio(){
    printf("\nREADING WITH STDIN\n");
    printheading();
    double x;
    while(scanf("%lg", &x) > 0){
        printline(x);
    }
    printf("-------------------------------------\n");

}
void readfile(char **argv){
    
}

int main(int argc, char **argv){
    assert(argc >= 2);
    if(strcmp(argv[1], "CMD") == 0){
        readcmd(argc, argv);
    } else if(strcmp(argv[1], "STDIN") == 0){
        readstdio();
    } else if(strcmp(argv[1], "FILE") == 0){
        readfile(argv);
    } else{
        fprintf(stderr, "ERROR: only supports reading with argv[1] being \"CMD\", \"STDIO\" or \"FILE\"\n");
    }
    return 0;
}