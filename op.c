// 20230824,katayama
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main() {
    FILE *file = fopen("trj0.g96", "r");
    FILE *fileID   = fopen("trj.bin", "wb");
    if (file == NULL || fileID == NULL) {
        perror("Unable to open the file");
        return 1;
    }

    char line[64];
    double data1, data2, data3;
    const char *start_string = "POSITIONRED";
    const char *end_string = "END";

    while (fgets(line, sizeof(line), file) != NULL) {
        if (strstr(line, start_string) != NULL) {
            while (fgets(line, sizeof(line), file) != NULL) {
                if (strstr(line, end_string) != NULL) {
                    break;
                } else {
                    sscanf(line, "%lf %lf %lf", &data1, &data2, &data3);
                    fwrite( &data1, sizeof(double), 1, fileID );
                    fwrite( &data2, sizeof(double), 1, fileID );
                    fwrite( &data3, sizeof(double), 1, fileID );
                }
            }
        }
    }

    fclose(file);
    fclose(fileID);
    return 0;
}
