/*
    AC - OpenMP -- SERIE
    fun_s.c
     rutinas que se utilizan en el modulo gengrupos_s.c 
****************************************************************************/
#include <math.h>
#include <float.h> // DBL_MAX
#include <stdlib.h>
#include <string.h>

#include "defineg.h"
           // definiciones

/**************************************************************************************
   1 - Funcion para calcular la distancia genetica entre dos elementos (distancia euclidea)
       Entrada:  2 elementos con NCAR caracteristicas (por referencia)
       Salida:  distancia (double)
**************************************************************************************/
double gendist (float *elem1, float *elem2)
{
	// PARA COMPLETAR
	int i;
	double dist = 0, res = 0;

	for(i = 0; i < NCAR; i++){
		res += pow((elem1[i] - elem2[i]), 2);
	}
	dist = sqrt(res);
	return dist;
}

/****************************************************************************************
   2 - Funcion para calcular el grupo (cluster) mas cercano (centroide mas cercano)
   Entrada:  nelem  numero de elementos, int
             elem   elementos, una matriz de tamanno MAXE x NCAR, por referencia
             cent   centroides, una matriz de tamanno NGRUPOS x NCAR, por referencia
   Salida:   popul  grupo mas cercano a cada elemento, vector de tamanno MAXE, por referencia
*****************************************************************************************/
void grupo_cercano (int nelem, float elem[][NCAR], float cent[][NCAR], int *popul)
{
	// PARA COMPLETAR
	// popul: grupo mas cercano a cada elemento
	int i, j;
	double min_dist = DBL_MAX, dist = 0;
	
	
    for(i = 0; i < nelem; i++) {
		for (j = 0; j < ngrupos; j++) {
			dist = gendist(elem[i], cent[j]);
			if (dist < min_dist) {
				min_dist = dist;
				popul[i] = j;
			}
		}
		min_dist = DBL_MAX;
	}	
}

/****************************************************************************************
   3 - Funcion para calcular la calidad de la particion de clusteres.
       Ratio entre a y b. El termino a corresponde a la distancia intra-cluster.
       El termino b corresponde a la distancia inter-cluster.
   Entrada:  elem     elementos, una matriz de tamanno MAXE x NCAR, por referencia
             listag   vector de NGRUPOS structs (informacion de grupos generados), por ref.
             cent     centroides, una matriz de tamanno NGRUPOS x NCAR, por referencia
   Salida:   valor del CVI (double): calidad/ bondad de la particion de clusters
*****************************************************************************************/
double silhouette_simple(float elem[][NCAR], struct lista_grupos *listag, float cent[][NCAR], float a[]){
    float b[ngrupos], s[ngrupos]; 
    
    // PARA COMPLETAR

	// aproximar a[i] de cada cluster: calcular la densidad de los grupos
    //		media de las distancia entre todos los elementos del grupo
    //   	si el numero de elementos del grupo es 0 o 1, densidad = 0

	//Calculamos la densidad media
	int numElementos=0;
	double distancia=0.0, cantidad=0.0, res = 0.0, calidad=0.0;
	int i,j,z,x;
    for (i=0; i<ngrupos; i++){
		numElementos = listag[i].nelemg;
		if(numElementos<2){
			a[i]= 0.0;
		}else{
			cantidad=0.0;
			distancia=0.0;
			for(j=0; j<numElementos; j++){
				 for(z= j+1; z<numElementos; z++){
					 distancia += gendist(&elem[listag[i].elemg[j]][0], &elem[listag[i].elemg[z]][0]); //distancia += gendist(elem[j], elem[z]); "de esta manera no pasas una direccion pasas un valor "
					 cantidad += 1.0;
				 }
			}
		a[i]= (float)(distancia/cantidad);
		}
		
		// Aproximar b[i] de cada cluster
		for (x=0; x <ngrupos; x++){
			if (i!=x){
				distancia += gendist(&cent[i][0], &cent[x][0]);
				cantidad += 1.0;
			}	
		}
		b[i]=(float)(distancia/cantidad);


		// calcular el ratio s[i] de cada cluster
		res= b[i]-a[i];
		s[i]= (float)(res/fmax(a[i],b[i]));
		calidad += s[i];
	}

	// promedio y devolver
    return calidad/(double)(ngrupos);

	}
/* Esto lo tenemos que documentar
    // aproximar b[i] de cada cluster
	for (i=0; i < ngrupos; i++){
		cantidad=0.0;
		distancia=0.0;
		for (j = 0; j <ngrupos; j++){
			if (i!=j){
				distancia += gendist(&cent[i][0], &cent[j][0]);
				cantidad += 1.0;
			}
			
		}
		b[i]=(float)(distancia/cantidad);
	}
	
	
	// calcular el ratio s[i] de cada cluster
	for (i = 0; i < ngrupos; i++){
		res= b[i]-a[i];
		s[i]= (float)(res/fmax(a[i],b[i]));
		calidad += s[i];
	}
	

	// promedio y devolver
    return calidad/(double)(ngrupos);
}
*/

/********************************************************************************************
   4 - Funcion para relizar el analisis de enfermedades
   Entrada:  listag   vector de NGRUPOS structs (informacion de grupos generados), por ref.
             enf      enfermedades, una matriz de tamaño MAXE x TENF, por referencia
   Salida:   prob_enf vector de TENF structs (informacion del análisis realizado), por ref.
*****************************************************************************************/

void analisis_enfermedades (struct lista_grupos *listag, float enf[][TENF], struct analisis *prob_enf)
{
	// PARA COMPLETAR
	// Realizar el análisis de enfermedades en los grupos:
	//		mediana máxima y el grupo en el que se da este máximo (para cada enfermedad)
	//		mediana mínima y su grupo en el que se da este mínimo (para cada enfermedad)
    float mediana, auxiliar;
    int i, j, k, z;
    int numeroElementos, minimo;

    for (i = 0; i < TENF; i++) {
        prob_enf[i].mmax = 0.0;
        prob_enf[i].mmin = 2.0;
        for (j = 0; j < ngrupos; j++) {
            numeroElementos = listag[j].nelemg;
            if (numeroElementos!=0) {
                float* array = malloc(numeroElementos * sizeof(float));
                for (k = 0; k < numeroElementos; k++) {
                    array[k]= enf[listag[j].elemg[k]][i];
                }
                for (k = 0; k < numeroElementos; k++) { //Usamos el método de ordenación selección: 
                    minimo = k;
                    for (z = k+1; z < numeroElementos; z++) {
                        if (array[z] < array[minimo]) {
                            minimo = z;
                        }
                    }
                    auxiliar = array[k];
                    array[k] = array[minimo];
                    array[minimo] = auxiliar;
                }
                mediana = array[numeroElementos/2];   //Sacamos la mediana
                if (mediana > prob_enf[i].mmax) {
                    prob_enf[i].mmax = mediana;
                    prob_enf[i].gmax = j;
                }
                if (mediana < prob_enf[i].mmin) {
                    prob_enf[i].mmin = mediana;
                    prob_enf[i].gmin = j;
                }
            }
        }
    }  
}


/***************************************************************************************************
   OTRAS FUNCIONES DE LA APLICACION
****************************************************************************************************/

void inicializar_centroides(float cent[][NCAR]){
	int i, j;
	srand (147);
	for (i=0; i<ngrupos; i++)
		for (j=0; j<NCAR/2; j++){
			cent[i][j] = (rand() % 10000) / 100.0;
			cent[i][j+(NCAR/2)] = cent[i][j];
		}
}

int nuevos_centroides(float elem[][NCAR], float cent[][NCAR], int popul[], int nelem){
	int i, j, fin;
	double discent;
	double additions[ngrupos][NCAR+1];
	float newcent[ngrupos][NCAR];

	for (i=0; i<ngrupos; i++)
		for (j=0; j<NCAR+1; j++)
			additions[i][j] = 0.0;

	// acumular los valores de cada caracteristica (100); numero de elementos al final
	for (i=0; i<nelem; i++){
		for (j=0; j<NCAR; j++) additions[popul[i]][j] += elem[i][j];
		additions[popul[i]][NCAR]++;
	}

	// calcular los nuevos centroides y decidir si el proceso ha finalizado o no (en funcion de DELTA)
	fin = 1;
	for (i=0; i<ngrupos; i++){
		if (additions[i][NCAR] > 0) { // ese grupo (cluster) no esta vacio
			// media de cada caracteristica
			for (j=0; j<NCAR; j++)
				newcent[i][j] = (float)(additions[i][j] / additions[i][NCAR]);

			// decidir si el proceso ha finalizado
			discent = gendist (&newcent[i][0], &cent[i][0]);
			if (discent > DELTA1)
				fin = 0;  // en alguna centroide hay cambios; continuar

			// copiar los nuevos centroides
			for (j=0; j<NCAR; j++)
				cent[i][j] = newcent[i][j];
		}
	}
	return fin;
}

