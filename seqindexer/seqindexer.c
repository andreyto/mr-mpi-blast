#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#define MAXBUFSZ 100000
#define MAXNQUERY 1000000
//#define TEST 1
#define ULONG_LONG_T unsigned long long

int main(int argc, char* argv[])
{
    if (argc < 2) {
        printf("Usage: lineindexer infile outfile\n\n");
        return 1;
    }
    
    char *inFileName = argv[1];
    FILE* p_inFile = fopen(inFileName, "r");
     
    ULONG_LONG_T seqBeginPos[MAXNQUERY];
    ULONG_LONG_T seqEndPos[MAXNQUERY];
 
    
    int nSeq = 0;
    char buff[MAXBUFSZ];
    
    char* outFileName = argv[2];
    FILE* p_outFile = fopen(outFileName, "w");
    
    /// produce an index of the file
    ULONG_LONG_T pos = ftell(p_inFile);
    fseek(p_inFile, 0, SEEK_END);
    ULONG_LONG_T fSize = ftell(p_inFile);
    rewind(p_inFile);
    printf("file size = %ld\n", fSize);
    
    /// Save index to file
    while ( fgets( buff, MAXBUFSZ, p_inFile ) != NULL ) {
        
        if (buff[0] == '>') {
            #ifdef TEST
            printf("%s\n", buff);
            #endif
            seqBeginPos[nSeq] = pos;
            nSeq++;
        }
        pos = ftell(p_inFile);       
        if (nSeq) seqEndPos[nSeq-1] = pos;
    }
    seqEndPos[nSeq-1] = fSize;
    printf("num total seqs = %d\n", nSeq);
    
    /// save
    int j;
    for (j = 0; j < nSeq; j++) {
        fprintf(p_outFile, "%ld,%ld\n", seqBeginPos[j], seqEndPos[j]);
    }
    
    ///
    /// Test
    ///
    #ifdef TEST
    
    for (j = 0; j < nSeq; j++) {
        printf("%ld %ld\n", seqBeginPos[j], seqEndPos[j]);
    }

    /// use the index to select a random line in a file
    int k = rand() % nSeq;
    fseek( p_inFile, seqBeginPos[k], SEEK_SET );
    fgets( buff, MAXBUFSZ, p_inFile );
    printf( "Random line from file %s", buff );

    fseek( p_inFile, seqBeginPos[0], SEEK_SET );
    fgets( buff, MAXBUFSZ, p_inFile );
    printf( "0 seq from file %s", buff );
    fseek( p_inFile, seqBeginPos[1], SEEK_SET );
    fgets( buff, MAXBUFSZ, p_inFile );
    printf( "1 seq from file %s", buff );
    
    fseek( p_inFile, seqBeginPos[2], SEEK_SET );
    while (ftell(p_inFile) < seqBeginPos[3]) {
        fgets( buff, MAXBUFSZ, p_inFile );
        printf( "2 seq from file %s", buff );
    }
    
    fseek( p_inFile, seqBeginPos[3], SEEK_SET );
    while (ftell(p_inFile) < seqBeginPos[4]) {
        fgets( buff, MAXBUFSZ, p_inFile );
        printf( "2 seq from file %s", buff );
    }
    
    fseek( p_inFile, seqBeginPos[4], SEEK_SET );
    while (ftell(p_inFile) < seqBeginPos[5]) {
        fgets( buff, MAXBUFSZ, p_inFile );
        printf( "2 seq from file %s", buff );
    }
    //int i;
    //ULONG_LONG_T sz = seqEndPos[i]-seqBeginPos[i]+3;
    //for (i = 0; i < nSeq; i++) {
        //char* buff2 = (char*) malloc(sizeof(char)*sz);
        //fseek( p_inFile, seqBeginPos[i], SEEK_SET );
        //ULONG_LONG_T res = fread(buff2, sizeof(char), sz+1, p_inFile);
        //if (res != sz+1) {
            //fputs ("Reading error 1",stderr); exit (3);
        //}
        //else {
            //printf("fread = %s\n", buff2);
        //}
        //free(buff2);
    //}
    #endif
    
    fclose(p_inFile);
    fclose(p_outFile);
 
    
    return 1;
}
