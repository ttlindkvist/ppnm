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
    printf("\nREADING WITH CMD - WRITING TO STDOUT\n");
    printheading();
    for(int i = 2; i<argc; i++){
        double x = atof(argv[i]);
        printline(x);
    }
    printf("-------------------------------------\n");
}
void readstdio(){
    printf("\nREADING WITH STDIN - WRITING TO STDOUT\n");
    printheading();
    double x;
    while(scanf("%lg", &x) > 0){
        printline(x);
    }
    printf("-------------------------------------\n");

}
void readfile(int argc, char **argv){
    assert(argc >= 4);
    FILE *in_stream = fopen(argv[2], "r");
    FILE *out_stream = fopen(argv[3], "a");

    fprintf(out_stream, "\nREADING FROM FILE - WRITING TO FILE\n");
    fprintf(out_stream, "-------------------------------------\n");
    fprintf(out_stream, "| x\t| cos(x)\t| sin(x)\t|\n");
    fprintf(out_stream, "-------------------------------------\n");

    double x;
    while(fscanf(in_stream, "%lg", &x) > 0){
        fprintf(out_stream, "| %g\t| %g\t| %g\t|\n", x, cos(x), sin(x));
    }
    fprintf(out_stream, "-------------------------------------\n");

    fclose(in_stream);
    fclose(out_stream);
}

int main(int argc, char **argv){
    assert(argc >= 2);
    if(strcmp(argv[1], "CMD") == 0){
        readcmd(argc, argv);
    } else if(strcmp(argv[1], "STDIN") == 0){
        readstdio();
    } else if(strcmp(argv[1], "FILE") == 0){
        readfile(argc, argv);
    } else{
        fprintf(stderr, "ERROR: only supports reading with argv[1] being \"CMD\", \"STDIN\" or \"FILE\"\n");
        exit(EXIT_FAILURE);
    }
    exit(EXIT_SUCCESS);
}