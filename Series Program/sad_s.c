/*
    AC - OpenMP -- SERIE
    fun_s.c
     rutinas que se utilizan en el modulo gengrupos_s.c 
****************************************************************************/
#include <math.h>
#include <float.h> // DBL_MAX
#include <stdlib.h>

#include "defineg.h"           // definiciones
#include <string.h>



/**************************************************************************************
   1 - Funcion para calcular la distancia genetica entre dos elementos (distancia euclidea)
       Entrada:  2 elementos con NCAR caracteristicas (por referencia)
       Salida:  distancia (double)
**************************************************************************************/
double gendist (float *elem1, float *elem2)
{
	// PARA COMPLETAR
	// calcular la distancia euclidea entre dos vectores
	double dist=0.0, rest=0.0;
	for(int i=0;i<40;i++){
		rest+=pow(elem2[i]-elem1[i],2);
	}
	dist=sqrt(rest);
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
	float new_min_dist=0.0, min_dist=DBL_MAX;
	int grup;
	//Iteramos cada elemento
	for(int i=0;i<nelem;i++){
		//Iteramos cada grupo/cluster
		for(int j=0;j<ngrupos;j++){
			//Obtenemos la distancia del elemento 'i' repecto al grupo 'j'
			new_min_dist=gendist(&elem[i][0],&cent[j][0]);
			//Solo para la primera pasada
			if(min_dist> new_min_dist){
				min_dist=new_min_dist;
				grup=j;
			}
		}
		min_dist=DBL_MAX;
		//Añadimos al grupo al que pertenece
		popul[i]=grup;
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
	double dist=0.0, cont=0.0, resta, calidad=0.0;
	int i,j,x;
	// aproximar a[i] de cada cluster: calcular la densidad de los grupos
    //		media de las distancia entre todos los elementos del grupo
    //   	si el numero de elementos del grupo es 0 o 1, densidad = 0
	
	//Iteramos sobre los clusteres
	for(i=0; i<ngrupos; i++){
		//Contador de el numero de distancias calculadas
		cont=0.0;
		dist=0.0;
		//Si el numero de elementos del cluster es mayor que 1 hacemos las operaciones necesarias
		if(listag[i].nelemg>1){
			//Para no repetir distancias con elementos ya usados
			for(j=0; j<listag[i].nelemg; j++){
				//Comparamos con un elemento el resto de elementos
				for(x=j+1; x<listag[i].nelemg; x++){
					dist+= gendist(&elem[listag[i].elemg[j]][0], &elem[listag[i].elemg[x]][0]);
					cont=cont+1.0;
				}
			}
			//Le asignamos la densidad al cluster
			a[i]=(float)(dist/cont); //        ---------------->CUIDADO AQUI<------------------
		}else{
			a[i]=0.0;
		}
	}
    
    // aproximar b[i] de cada cluster
	for(i=0; i<ngrupos; i++){
		cont=0.0;
		dist=0.0;
		for(j=0; j<ngrupos; j++){
			if(i!=j){
				dist+= gendist(&cent[i][0], &cent[j][0]);
				cont= cont +1.0;
			}
		}
		b[i]=(float)(dist/cont); //        ---------------->CUIDADO AQUI<------------------
	}
	
	// calcular el ratio s[i] de cada cluster
	for(i=0; i<ngrupos; i++){
		resta= b[i] - a[i];
		s[i]= (float)(resta/fmax(a[i],b[i])); //        ---------------->CUIDADO AQUI<------------------
		calidad+=s[i];
	}
	
	
	// promedio y devolver
    return calidad/(double)(ngrupos);//        ---------------->CUIDADO AQUI, NS SI SE PUEDE HACER ESE CASTING<------------------
}

/*
//funcion pa comparar
#include <stdio.h>
#include <stdlib.h>
int comp (const void * elem1, const void * elem2) 
{
    int f = *((int*)elem1);
    int s = *((int*)elem2);
    if (f > s) return  1;
    if (f < s) return -1;
    return 0;
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
	float ar[MAXE], aux;
	int y=0, k=0;
	int indice,indice2, muro=0;
	//Para recorrer cada enfermedad
	for(int i=0; i<TENF; i++){
		muro=0;
		//Para recorrer cada cluster
		for(int j=0; j<ngrupos; j++){
			//Inicializamos el array a cero
			bzero(ar, MAXE);
			//Para obtener el indice de cada elemento asociado al cluster 'j' y asociar las probabilidades
			for(int x=0; x<listag[j].nelemg; x++){
				//Obtenemos el indice del elemento 'x' del cluster 'j'
				indice= listag[j].elemg[x];
				
				//Obtenemos la probabilidad de esa enfermedad en el elemento 'x'
				aux= enf[indice][i];

			//ar[x]= enf[indice][i];
				//Insertar ordenadamente en el array ar
				y=0;
				while(y<x && aux>ar[y]){
					y=y+1;
				}
				if(aux<ar[y]){
					k=listag[j].nelemg;
					while(y<k){
						ar[k]=ar[k-1];
						k--;
					}
				}
				ar[y]=aux;
				
			}
			indice2= listag[j].nelemg/2;

			if(listag[j].nelemg>0){
				if(muro==0){ //Se puede? -> &prob_enf[i].mmin==NULL
					prob_enf[i].mmin=ar[indice2];
					prob_enf[i].mmax=ar[indice2];
					prob_enf[i].gmin=j;
					prob_enf[i].gmax=j;
					muro++;
				}else{
					if(prob_enf[i].mmin>ar[indice2]){
						prob_enf[i].mmin=ar[indice2];
						prob_enf[i].gmin=j;
					}
					if(prob_enf[i].mmax<ar[indice2]){
						prob_enf[i].mmax=ar[indice2];
						prob_enf[i].gmax=j;
					}
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

