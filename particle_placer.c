#include <stdio.h>

// simple cubic particle placer
int main(int argc, char* argv[])
{
	if(argc != 2)
        printf("usage: scplacer.exe N, where N is the amount of cubes per dimension\n");
    int nPerDim = 0;
    sscanf(argv[1], "%d", &nPerDim);
    if(nPerDim == 0)
        printf("usage: scplacer.exe N, where N is the amount of cubes per dimension\n");

    FILE* file;
	char buffer[16] = "sc%d.txt", filename[16] = "";
    sprintf(filename, buffer, nPerDim); // e.g. "sc10.txt"

	file = fopen(filename, "w");
    
	// print the header that is needed for the visualization program
	fprintf(file, "%d\n0.000000 %lf\n0.000000 %lf\n0.000000 %lf\n",
        nPerDim * nPerDim * nPerDim, (double) nPerDim, (double) nPerDim, (double) nPerDim);

	// now make the loop to place particles in a cubic lattice (with offset)
    for (int x = 0; x < nPerDim; x++) {
        for (int y = 0; y < nPerDim; y++) {
            for (int z = 0; z < nPerDim; z++) {
				fprintf(file, "%lf %lf %lf\n", x + 0.5, y + 0.5, z + 0.5);
            }
        }
    }

	fclose(file);
    return 0;
}
