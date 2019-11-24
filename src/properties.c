#include "properties.h"
void handle_error_with_exit(char*error_string){//uccide il processo dopo essersi accorto di un errore
    if(error_string==NULL){
	printf("error string is NULL\n");
	perror("");
        exit(EXIT_FAILURE);
    }
    printf("%s",error_string);
    //verifica che la stringa contenga il carattere newline
    //fprintf(file_log,"%s",error_string);
    exit(EXIT_FAILURE);
}