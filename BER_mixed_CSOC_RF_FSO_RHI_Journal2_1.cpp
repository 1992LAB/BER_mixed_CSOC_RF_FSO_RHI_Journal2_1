
  
// Version final
   
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include<conio.h>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */

//****** Code *******
//#include "R1_2_J4.h"
#include "R1_2_J6.h"
//#include "R1_2_J8.h"
//#include "R1_2_J10.h"
//#include "R1_2_J18.h"
//#include "R1_2_J22.h"

//#include "R2_3_J4.h"
//#include "R2_3_J6.h" 
//#include "R2_3_J10.h"

//#include "R4_5_J4.h"
//#include "R4_5_J6.h"
//#include "R4_5_J12.h"

//#include "R3_4_J11.h"
//#include "R3_4_J17.h"

//#include "R8_9_J3.h"
//#include "R8_9_J5.h"
//  -----------------

////////////// 
//#include "R1_2_J5_An.h"
//#include "R2_3_J5_An.h"
//#include "R1_2_J9_An.h"

//#include "R2_3_J3_4.h"
//#include "R2_3_J5_8.h"
// #include "R2_3_J9_16.h"
//#include "R2_3_J17_32.h"
//#include "R4_5_J3_8.h"
//#include "R4_5_J5_16.h"
//#include "R4_5_J9_32.h"
//#include "R8_9_J3_16.h"
//#include "R8_9_J5_32.h"
//#include "R16_17_J3_32.h"
//#include "R16_17_J4_64.h"
 
#define pi 3.14159265358979323846

// Allocation d'une matrice :----------------------------------------------------------
int **Alloc_Mat(int l, int c){  // n :nbre ligne , m: nbre colonne 
	int **M;
	M = (int**) malloc( l * sizeof(int*) );    
	for(int i= 0; i<l ; i++) {  M[i]=(int*) malloc( c * sizeof(int) );} 
	return M;
} 
   
// Allocation d'un tableau de int -----------------------------------------------------
int *Alloc_Tab(int l){  // n :nbre ligne , m: nbre colonne 
	int *T;
	T = (int*) malloc( l * sizeof(int) );    
	return T;
}
   
// Allocation d'une matrice de float :----------------------------------------------------------
double **Alloc_Mat_double(int l, int c){  // n :nbre ligne , m: nbre colonne 
	double **M;
	M = (double**) malloc( l * sizeof(double*) );    
	for(int i= 0; i<l ; i++) {  M[i]=(double*) malloc( c * sizeof(double) );} 
	return M;
}

// Allocation d'un tableau de float :----------------------------------------------------------
double *Alloc_Tab_double(int l){  // n :nbre ligne , m: nbre colonne 
	double *T;
	T = (double*) malloc( l * sizeof(double) );    
	return T;
}
   
// Libérer une mémoire du tableau -----------------------------------------------------
void free(void * bloc);   //free(entier); entier = NULL; (Forcer l"entier à être Null)

// Liberer une matrice de float (matrice) : ----------------------------------------------------
 void free_double_mtx(double ** mtx,int line){
    for(int i=0;i<line;i++)
   { free(mtx[i]);  mtx[i]=NULL ; }
    free(mtx);  mtx=NULL ;  
 } // ---------------------------------------------------------------------------------   
// Pour Liberer la mémoire de int (matrice) : ---------------------------------------------------
 void free_int_mtx(int** mtx,int line){
        for(int i=0;i<line;i++){
           free(mtx[i]);  mtx[i]=NULL ;
          }
        free(mtx); 	mtx=NULL ;  
 }// ----------------------------------------------------------------------------------
 // Pour Liberer la mémoire de int (tableau) : ----------------------------------------
 void free_int_tab(int* tab){
        free(tab); 	tab=NULL ;  
 } // ---------------------------------------------------------------------------------------------
 
// Pour Liberer la mémoire de float (tableau) : ---------------------------------------------------
 void free_double_tab(double* tab){
        free(tab); 	tab=NULL ;  
 }
 
// initialiser tous les éléments d'une matrice à 0 ------------------------------------
void initialiser_matrix(int **matrix,int l,int c){
	for(int i=0;i<l ;i++){ for(int j=0;j<c ;j++) {   matrix[i][j]=0;  } }
} //-----------------------------------------------------------------------------------

//initialiser tous les éléments d'une matrice float à 0.0 ----------------------------
void initialiser_matrix_double(double **matrix,int l,int c){
	for(int i=0;i<l ;i++){ for(int j=0;j<c ;j++) {   matrix[i][j]=0;  } }
}
//-----------------------------------------------------------------------------------

// initialiser tous les éléments d'un tableau à 0 -------------------------------------
void initialiser_Vect(int *vect,int c){   for(int i=0;i<c ;i++)  vect[i]=0;  }
//-------------------------------------------------------------------------------------	

// initialiser tous les éléments d'un tableau à 0 -------------------------------------
void initialiser_Vect_double(double *vect,int c){   for(int i=0;i<c ;i++)  vect[i]=0;  }
//-------------------------------------------------------------------------------------	

// Afficher une matrice :--------------------------------------------------------------
void affiche_matrix(int **matrix,int l,int c) {
    for(int i=0;i<l;i++) {
        for(int j=0;j<c;j++) {  printf("%3d",matrix[i][j]);  }
        printf("\n");
    }
}//------------------------------------------------------------------------------------

// Afficher un vecteur int :---------------------------------------------------------------
void affiche_vect(int *vect,int c) {
    for(int i=0;i<c;i++)  printf("%3d",vect[i]);
}// ----------------------------------------------------------------------------------

// Afficher un vecteur float :---------------------------------------------------------------
void affiche_vect_double(double *vect,int c) {
    for(int i=0;i<c;i++)  printf("%3e ",vect[i]);
}// ----------------------------------------------------------------------------------

// Max d'un tableau ------------------------------------------------------------------
int maxTab(int *T,int taille){
	int max=T[0];
	for (int i=0;i<taille;i++){ if(T[i]>max) {  max=T[i]; }	}
	return max;
} //----------------------------------------------------------------------------------

// Max d'une matrice :----------------------------------------------------------------- 
int maxMatrice(int **T,int l,int c){
	int max=T[0][0];
	for (int i=0;i<l;i++){ for (int j=0;j<c;j++){ if(T[i][j]>max){ max=T[i][j]; }  }  }
	return max;
} //----------------------------------------------------------------------------------

// fonction Décalage à gauche 'un tableau d'entiers  ---------------------------------
void decalerG_Tabl( int *T, int taille ){
	for(int i=1; i < taille ;i++){
		T[i-1]=T[i];
	}
	T[taille-1]=0;
}

// fonction Décalage à gauche d'une Matrice d'entiers  ------------------------------- 
void decalerG_Mat( int **T, int l , int c ){
	for(int i=0; i < l ;i++){
        for(int j=1; j < c ;j++){
        	T[i][j-1]=T[i][j];	
		}
    	T[i][c-1]=0;
	}
} // ----------------------------------------------------------------------------------------

// fonction XOR de deux bits : ---------------------------------------------------------------
 int  xxor(int pp, int ff)    // p et q sont des bits
 {   return ( pp ^ ff);  } // ------------------------------------------------------------------------------


// code_conv_1_n ----R=1/2-----Systématique-----Systématique-----Systématique-----------------------
int* code_conv_1_n(int *Info_In,int **g, int dim2) {  // retourne un tableau : Information codée 

    // dim2 : tailled Info_In
    
	// Tableau à retourner : code Word :
	int *Code_Out = (int*) malloc( 2*dim2 * sizeof(int) );
    for(int i=0; i<2*dim2; i++){  Code_Out[i]=0; } // initialisation
    
    int N=2 ; // (Toujours =2) (le nbre de polynome générateurs (Une g0=1, g1=[1 0 1 1 ......0 1] ) 
	//int m=dim2-1;  // la longueur de g
    
    // Initialization des registres mémoire à zeros
    int **mm = (int**)malloc(N * sizeof(int*));    
    for(int i=0; i<N; i++){  mm[i]=(int*)malloc( dim2 * sizeof(int)) ;  }	
	initialiser_matrix(mm ,N , dim2);
	
	/*printf("\nMatrice mm:\n");
	affiche_matrix(mm ,N , dim2);*/
	
    // Le 1 ér registre prend le premier bit de l'Information :
	for(int i=0; i < N ;i++){
    	mm[i][0]=Info_In[0];
	}
	
	/*printf("\nMatrice mm 1:\n");
	affiche_matrix(mm ,N , dim2);*/
     
    for(int kk=0; kk < dim2 ;kk++){   // * * Codage bit à bit
    // le codage : -----------------------------------------------------------
        int *s= (int*)malloc(N * sizeof(int));
	
    	for(int i=0; i < N ;i++){
        	s[i]=0;
        	for(int j=0; j < dim2 ;j++){
    	    	if(g[i][j]!=0){
    		    	s[i]=xxor(s[i],mm[i][j]);
		    	}
	    	}
	    }
    
    // Le signal du sortie ----------------------------------------------------
        for(int t=N-1; t>=0 ; t--){
        	Code_Out[((N)*(kk+1))-(t+1)]=s[(N-1)-t];
    	}
	
	// On doit faire le décalage des bits des registre vers la droite à chaque fois 
        	for(int i=0; i<N ; i++){
                for (int j=dim2-1; j>0; j--){
			        mm[i][j]=mm[i][j-1];
		        }
		    //Le 1 ér registre prend le prochain bit de l'Information
		        if(kk < (dim2-1)){
		            mm[i][0]=Info_In[kk+1];
		        }
            }
         free_int_tab(s) ;  
	}
        /*printf("\nCode_Out :\n");
	    affiche_vect(Code_Out , dim2);*/
    	return Code_Out;
    	
		free_int_mtx(mm, N );
    	free_int_tab(Code_Out);
}
// -----------------------------------------------------------------------------------------------------------------

// fonction codage_All_Rates ---------------------------------------------------------------------------------------
void codage_conv(int *Info_In, int *Code_Out, int **gg, int taille){
    
	// kk: Nbre de vecteurs d'entrées
    // nn: Nbre de vecteurs de sorties
    // Info_In = rand(1,kk*length(gg))>0.5 ;  // l'Information binnaire à coder
    // gg: Matrice contient tous les polynomes générateurs du code
    // taille : taille de l'information à coder  , taille length(Info_In)
    // long_gg : longueur de la matrice des polynômes générateurs 
    
	//le nombre des éléments d'un tableau = le rapport entre la taille du tableau et la taille d'un de ses éléments.
    int SizeInfo = taille;
    //printf("\nSizeInfo = %d \n", SizeInfo);
    
    //int *Code_Out;   // tableau à retourner // SizeCodeOut = SizeInfo * nn;
    int SizeCodeOut = SizeInfo + (SizeInfo/k) ;
    //printf("\nSizeCodeOut = %d \n", SizeCodeOut);
    
    //Code_Out = ( int* )malloc( SizeCodeOut * sizeof(int) );
    //printf("\nAlloc Code_Out");
    
	// Découper l'informations en kk séquences d'entrée (nbre d'entrées) ----------------------------------------
    // Allocation de u
    int **u = ( int** )malloc( n * sizeof(int*) );    
    for( int i=0; i<n; i++ ){  u[i] = ( int* )malloc( (SizeInfo/k) * sizeof(int) ) ;  }
    initialiser_matrix(u ,n , (SizeInfo/k));
	
	/*printf("\nMatrice u 1 :\n");
	affiche_matrix(u ,nn , (SizeInfo/kk));*/
	
	for(int i=k-1; i >=0 ;i--){
    	for(int j=0; j < (SizeInfo/k) ;j++){
    		u[(k-1)-i][j]=Info_In[(k)*(j+1)-(i+1)]; 
		}
	}
	
	/*printf("\nMatrice u 2 :\n");
	affiche_matrix(u ,nn , (SizeInfo/kk));*/
    // Codage :--------------------------------------------------------------------------------------------------
    // Coder chaque séquence d'entrée avec  code_conv_1_n
    
    int l=2;    // toujours =2
	int **c = ( int** )malloc( k * sizeof(int*) );    
    for( int i=0; i<k; i++ ){  c[i] = ( int* )malloc( l*(SizeInfo/k) * sizeof(int) ) ;  }
    initialiser_matrix(c ,k , l*(SizeInfo/k));
     
        /*printf("\nc :\n"); 
	affiche_matrix(c ,kk , l*(SizeInfo/kk));*/
	
	int **g_tab = ( int** )malloc( l * sizeof(int*) );    
    for( int i=0; i<l; i++ ){  g_tab[i] = ( int* )malloc( (SizeInfo/k) * sizeof(int) ) ;  }

	// Copier le poly num i+1 en g_tab 
//---------------------------------------------------------------------------	
	// printf("long_gg=%d", long_gg);
	for(int ii=0; ii<k ;ii++){
		
	    initialiser_matrix(g_tab ,l , (SizeInfo/k));
	    g_tab[0][0]=1;
	
	    /*printf("\ng_tab 0 :\n") ;
        affiche_matrix(g_tab,l,(SizeInfo/kk));*/
    
		for(int jj=0; jj<lg ;jj++){
			g_tab[1][jj]=gg[ii+1][jj];
		}
	   
	    /*printf("\n g_tab 1 : \n") ;
        affiche_matrix(g_tab,l,(SizeInfo/kk));*/
		
	    //printf("\n 	(SizeInfo/kk)=%d \n",(SizeInfo/kk));
	     
	    c[ii]=code_conv_1_n(u[ii],g_tab,(SizeInfo/k));
	    
	    
	//	}
    }
//---------------------------------------------------------------------------	
	
	/*printf("\nMatrice c: le mot de code :\n") ;  // le mot de code: comme si j'ai un code de rendement 1/2
    affiche_matrix(c,kk,(l*(SizeInfo/kk))); */
    
    // ---------------------------------------------------------------------------------------------------------
    // creer une matrice qui contient juste les bits de redandance :
    int **cc = ( int** )malloc( k * sizeof(int*) );    
    for( int i=0; i<k; i++ ){  cc[i] = ( int* )malloc( (SizeInfo/k) * sizeof(int) ) ;  }
    initialiser_matrix(cc ,k , (SizeInfo/k));
	
	/*printf("\ncc 0:\n") ;  
    affiche_matrix(cc,kk,(SizeInfo/kk)); */
    
	for(int i=0; i<k ;i++){
		for(int j=0; j<(SizeInfo/k) ;j++){
			cc[i][j]=c[i][2*(j+1)-1];
		}
	}
	
   /* printf("\ncc 1:\n") ;  
    affiche_matrix(cc,kk,(SizeInfo/kk));   */
    
	// La derniere séquence (de redandance) qui est égale au xor de toutes les séquences
    // codée (je l'ajoute à la matrice u)
    // calculer le xor de toutes les séquences de c (les séq de redandance) = séq
    // de parité du code : c_red et la mettre dans la derniére ligne de u

	int *c_red =  ( int* )malloc((SizeInfo/k) * sizeof(int) ); 
    // Initialisation du vecteur c_red par la 1 ére ligne de la matrice des vecteurs de redandance c :
	for (int i=0; i<(SizeInfo/k)  ;i++){
    	c_red[i] = cc[0][i];
	}
	
	// XOR de toutes les lignes des vecteurs de redandance de la matrice c :
    for (int i=1; i<k ;i++){
    	for (int j=0; j<(SizeInfo/k) ;j++){
    		c_red[j] = xxor(c_red[j],cc[i][j]);
		}
	}
   
    for(int j=0; j < (SizeInfo/k) ;j++){
    	u[(n-1)][j]=c_red[j];
    }
    
    
	//  en série
    //  Code_out : Mot de code 
	for(int i=0 ;i<n ;i++){
		Code_Out[i]=u[i][0];	
	}

	for(int j=1; j<(SizeInfo/k) ; j++){
	    for (int i=j*n ; i<(j+1)*n ; i++){ 
			Code_Out[i]=u[i-(j*n)][j];
		}
	}
    
	/*printf("\nCode_Out 1:\n");
	affiche_vect(Code_Out , SizeCodeOut);*/
    //return Code_Out; 
    
	free_int_tab(c_red ) ;
    free_int_mtx(cc, k );
    free_int_mtx( g_tab, l);
	free_int_mtx(c, k );
    free_int_mtx(u, n );
	//free_int_tab(Code_Out);  
} 
//---------------------------------------------------------------------------------------------------------

// SumEqua_all_rates : sommes (XOR) des équations Orthogonales----------------------------------------------------------------
int** SumEqua_all_rates( int **A , int d1 , int d2) {   // d1 : nbre de ligne , d2: nbre de colonne, Retourne matrice cc(d2/J, J)
	//[d1,d2]=size(A);  // d1=J*k , d2= J*k+1 
	 int *c = (int*) malloc( d1 * sizeof(int) );
     for(int i=0 ;i<d1 ;i++){
     	c[i]=0; 
     	for(int j=0 ;j<d2 ;j++){
     		c[i]=xxor(c[i],A[i][j]);
		 }
	 }
    
	// Déviser c (tableau de taille d2) ==> en matrice InfoBmsg de taille ( d2/J , J ) 
    int **cc = ( int** )malloc( (d1/J) * sizeof(int*) );    
    for( int i=0; i < (d1/J) ; i++ ){  cc[i] = ( int* )malloc( J * sizeof(int) ) ;  }
    
    initialiser_matrix(cc , (d1/J) , J );
	
    for( int j=0; j < (d1/J) ; j++ ){
	    for(int i = J-1; i >=0 ;i--){
    		cc[j][(J-1)-i]=c[J *(j+1)-(i+1)]; 
   		}
	}    
	
	return cc; 
	
	free_int_mtx(cc, (d1/J) );
	free_int_tab(c);
	
} // --------------------------------------------------------------------------------------------------------------------

// fonction sert à remplacer un vecteur dans les éauations Orthogonales  --------------------------------------------------
int **Remplacer_Vect_Equ( int** E , int** s,  int d1, int d2){
    // s2 : Matrice des Equations Orthogonales :
	// E : Les vecteurs à remplacer dans les équations Orthogonles 
    // d1 : Nbre de ligne de s2 
	// d2 : Nbre de colonne de s2
	// [d1,d2] = size(s2);     //  d1 : ligne , d2 : colonnes
    // E : Matrice des vecteurs à remplacer dans les équations Orthogonales (contient n vecteur de taille lg ) 
             
    // Remplacement du 1er vecteur : -------------------------------------------         
    int **B= ( int** )malloc( d1 * sizeof(int*) );
	for( int i=0; i < d1 ; i++ ){  B[i] = ( int* )malloc( d2* sizeof(int) ); }
    initialiser_matrix(B ,d1 , d2 );
   
    for( int i=0; i < d1 ; i++ ){
        for( int j=0; j < J ; j++ ){
		    if(s[i][j]!=0){
				B[i][j]= E[0][s[i][j]-1]; //E(1,(s2(i,j))); // 
			}
        }                         
    }
            		
    // remplacement des autres vecteurs sauf le dernier (qui -------------------
    // contient les bits de redandances)	
	for( int ii=1; ii < k ; ii++ ){
       	for( int i=0; i < d1 ; i++ ){
            for( int j=J*ii; j <J*ii+J ; j++ ){
             	if (s[i][j]!=0){
             		B[i][j]=E[ii][s[i][j]-1]  ;  // E(k+1,s2(i,j));
				}
			}
		}
	}
			
    // pour le vecteur des bits de redandances --------------------------------
    for( int i=0; i < d1 ; i++ ){
        B[i][d2-1]=E[k][s[i][d2-1]-1];  //E(kk+1,s2(i,j));
    }
        
	return B;
	
	free_int_mtx(B, d1 );          
}

// fonction de décodage : Decoder_MLD_HIHO_R_1_2_all_rates --------------------------------------------------------------------
void Decoder_MLD_HIHO(int *InfoDecod, int *Info_Demod_recu, int N, int **A,int **gg){
	// k : nbre d'entrées 
	// n : nbre de sortie 
	// Info_Demod_recu : vecteur reçu démodulé 
	// N : taille du vecteur reçu
	// lg : taille Matrice des polynomes générateurs
	// A :  matrice des équations de parité 
	// gg : matrice des polynômes générateurs
	// J : nbre des Equations orthogonales (dans chaque Ensemble d'equations)
	// tML : capacité de correction du code 
	// m : mémoire du code
	
	// Longueur des mots de codes ----------------------------------------------------------------
    int dim2 = n*((N/k)+m);
    //int nbreBloc = N/(k*lg) ;  // Nbre de blocs (dans l'information démodulé reçu)
	int dimInfoDemo = dim2 ;  // Taille de l'information démodulé reçu
    //printf("\ndim2=%d , dimInfoDemo=%d \n", dim2, dimInfoDemo);
	  
	// l'Information decode à retourner : -------------------------------------------------------
	/*int *InfoDecod=(int*)malloc( N * sizeof(int));
	initialiser_Vect(InfoDecod, N);*/
	/*printf("\nInfoDecode : \n");
    affiche_vect(InfoDecod, N );*/
      
	// La Matrice qui va contenir les séquence décodées ------------------------------------------
    /*int **u = Alloc_Mat( nbreBloc ,(k*lg) );
	initialiser_matrix(u , nbreBloc , (k*lg));*/
	// séquence décodé de chaque entrée
	int **u = Alloc_Mat( k ,(N/k) );
	initialiser_matrix( u ,k ,(N/k));
	
	/*printf("\nMatrice u:\n");
    affiche_matrix(u , k , (N/k) );*/
	     
    // ue : Matrice des vecteurs à remplacer dans les équ Orth ----------------------------------
    /*int **ue = Alloc_Mat(k ,lg );
    initialiser_matrix(ue , k ,lg );
    printf("\nMatrice ue:\n");
    affiche_matrix(ue , k ,lg );*/
      	                          
    // Diviser l'info démodulé en sous séquence (taille=taille d'un code) :-----------------------
    // Matrice contient les sous séquences
    int **Info_recuu = Alloc_Mat(n , dimInfoDemo/n);
    initialiser_matrix( Info_recuu, n, dimInfoDemo/n);
   /* printf("\nMatrice Info_recuu:\n");
    affiche_matrix(Info_recuu , n ,dimInfoDemo/n);*/
    	
	for(int i=n-1; i >=0 ;i--){
         	for(int j=0; j < (dimInfoDemo/n) ;j++){
    		    Info_recuu[(n-1)-i][j] = Info_Demod_recu[n*(j+1)-(i+1)]; 
		}
	}

	/*printf("\nMatrice Info_recuu:\n");
    affiche_matrix(Info_recuu , n ,dim2/n );*/
    
		// Sans la sequence de redandance (la sequence de redandance est la derniére dans r)--------------
		int** R1 =  Alloc_Mat( n-1 , dimInfoDemo/n );
        
        for(int i=0; i<n-1 ;i++){
        	for(int j=0; j<dimInfoDemo/n ;j++){
        		R1[i][j]=Info_recuu[i][j];
			}
		}
		/*printf("\nR1:\n");
        affiche_matrix(R1 , n-1 ,dim2/n);*/
        
        // vecteur de chaque séquence décodée-------------------------------------------------------------
        int **EE1 = Alloc_Mat(k , dimInfoDemo/n); 
		initialiser_matrix(EE1 , k , dimInfoDemo/n);
        /*printf("\nEE1 :\n");
        affiche_matrix(EE1 , k , dim2/n);*/
		         
        // vecteur d'erreur -------------------------------------------------------------------------------
		int **E = Alloc_Mat(k , dimInfoDemo/n); 
		initialiser_matrix(E , k , dimInfoDemo/n);
        /*printf("\nE :\n");
        affiche_matrix(E , k , dim2/n);*/
		   
        //Matrice ldes vecteurs à remplacre dans les équations Orthogonales  :----------------------------
		int** e =  Alloc_Mat(n ,lg );
        initialiser_matrix(e , n ,lg );
		for(int i=0; i<n ;i++){
        	for(int j=0; j<lg ;j++){
        		e[i][j]=Info_recuu[i][j];
			}
		}
		/*printf("\ne :\n");
        affiche_matrix(e , n ,lg);*/
     //getch();
        // check sums : Ap=A  -------------------------------------------------------------------------------------------
		// e : matrice contient les (nn-1) vecteurs d'erreurs à remplacer dancles équations de Parité Orthogonales--------
        //int **Ap = Alloc_Mat(k*J , k*J+1);  
        //initialiser_matrix(Ap, k*J , k*J+1) ;
		int **Ap=Remplacer_Vect_Equ( e , A,  J*k,  J*k+1);         
        /*printf("\nAp :\n");
        affiche_matrix(Ap , k*J , k*J+1);*/

		// sommation des bits des équ orth : -----------------------------------------------------------------------------
        int **c=SumEqua_all_rates( Ap , J*k, J*k+1); //--> la rendre Matrice chaque ligne contient c des equ Orth correspondantes 
        
	 /*printf("\nc :\n");
	 affiche_matrix(c , k ,J);*/
	 
		// les syndrome : initialisation ---------------------------------------------------------------------------------
		int *som_syndrome=(int*)malloc(lg*sizeof(int));
        initialiser_Vect(som_syndrome, lg) ;
        /*printf("\nsom_syndrome :\n");
        affiche_vect(som_syndrome , lg); */
		
		// Décodage :  ----------------------------------------------------------------------------------------------------       
	    for(int kk=0; kk<N/k ; kk++){
	    	
	    	//calculer la somme dans Majority logic gate ------------------------------------------------------------------
			int *som=(int*) malloc(k*sizeof(int));
	    	for(int i=0; i<k ;i++){
                som[i]=0;
                for(int j=0; j<J ;j++){       // on considére pas la premiére équation 
                     som[i]=som[i]+c[i][j] ;   // le nombre des équation pour faire
                }
            }
          /*printf("\nsom: \n");
          affiche_vect(som,k);
          getch(); */
            // estimation : ----------------------------------------------------------------------------------------------
             for(int i=0; i<k ;i++){       
				if(som[i]>tML){  // changer --> je dois corrigé le bit donc EE1(k)=1 --> bit corrigé : EE1(k)+r(k)
					E[i][kk]=1;
				}
				else{       // %if som<=tML  --> ne pas changer -->   pas d'erreur r(k) est correct donc EE1(k)=0
					E[i][kk]=0;
				}
			 }
			 
			 // the estimated information bit : correction  --------------------------------------------------------------
             for(int i=0; i<k ;i++){
             	EE1[i][kk]=xxor(R1[i][kk],E[i][kk]);
			 }      
			
		  /*printf("\nMatrice EE1 :\n");
          affiche_matrix(EE1 , k , dim2/n);	  
		  getch();*/
			 
			 // placer E dans som_syndrome -------------------------------------------------------------------------------
             // E(1) must be substructed from each syndrome bit it affects
              for(int j=0; j<k ;j++){
              	  for(int i=0; i<lg ;i++){
              	  	  if(gg[j+1][i]==1){
              	  	  	  som_syndrome[i]=xxor( som_syndrome[i] ,R1[j][kk] );
					  }
				  }
			  }      
		  /*printf("\nsom_syndrome  :\n");
          affiche_vect(som_syndrome ,lg);  
		  getch();*/		
			  // Décalage à gauche de som_syndrome: ----------------------------------------------------------------------
              decalerG_Tabl( som_syndrome, lg );
              
			  /*printf("\nsom_syndrome  :\n");
              affiche_vect(som_syndrome , lg);*/ 
			  
		  
			  // Prendre les valeurs som_syndrome associé a s2(i,1)=1 : --------------------------------------------------
			  // Utiliser un seul vecteur de syndrome (Attention!)
			  int **som_syn = Alloc_Mat(k,J);  //
			  
			  for(int j=0; j<k ;j++){
					int nn=0;
				    for (int i=0 ;i<lg ;i++){
				    	if (gg[j+1][i]==1 ){
				    		som_syn[j][nn]=som_syndrome[i];
				    		nn++;
						}
					}  
			  }
		      
			  int **som_syndrome_A=som_syn;
		  /*printf("\nsom_syndromeA  :\n");
          affiche_matrix(som_syndrome_A , k , J);   
          getch();*/
			  // Décalage à gauche de r: ----------------------------------------------------------------------
              decalerG_Mat( Info_recuu, n ,dimInfoDemo/n );
          /*printf("\nInfo_recuu :\n");
          affiche_matrix(Info_recuu , n ,dim2/n);*/
              
              // e prend les (1--> lg) vecteurs-----------------------------------------------------------------
              for(int i=0; i<n ;i++){
			      for(int j=0; j<lg ;j++){
        		     e[i][j]=Info_recuu[i][j];
			      }
	          }	
			   
			/*printf("\ne :\n");
			affiche_matrix(e , n ,lg);   
		    getch();*/
               // Calcul de Ap ------------------------------------------------------------------------------------
			  initialiser_matrix(Ap , J*k,  J*k+1); 
			  Ap=Remplacer_Vect_Equ( e , A,   J*k,  J*k+1);                 
              /*printf("\nAp :\n");
			  affiche_matrix(Ap , J*k,  J*k+1);*/
			  
			  // --------------------------------------------------------------------------------------------------
              initialiser_matrix(c , k , J);
			  c  = SumEqua_all_rates( Ap , J*k, J*k+1);
			  /*printf("\nc :\n");
			  affiche_matrix(c , k ,J);*/
			  
              
			  // XOR c avec som_syn ------------------------------------------------------------------------------
              for(int i=0; i<k ;i++){
              	  for(int j=0; j<J ;j++){
              	  	   c[i][j]=xxor(c[i][j],som_syn[i][j]);
			    	}
			  }
          /*printf("\n matrice c :\n");
		  affiche_matrix(c , k ,J); 
		  getch();*/
		  
		//free_int_mtx(som_syndrome_A,k);
		free_int_mtx(som_syn,k);
		//free_int_mtx(som_syndrome_A,k);
		free_int_tab(som) ;  
		  
		}
		
	
    // Information Decodée ---------------------------------------------------------------------------
    /*for(int i=0 ;i<N/k ;i++){
		InfoDecod[i]=EE1[0][i];	
	}
	   
    for (int i=1 ; i<k ; i++){
	    for(int j=i*N/k; j<(i+1)*N/k ; j++){
         	InfoDecod[j]=EE1[i][j-(i*N/k)];
        }
    } */
    
        /*int *uu_s = (int*)malloc(k*lg * sizeof(int));
		initialiser_Vect(uu_s, k*lg);*/
		
		for(int i=0 ;i<k ;i++){
		    InfoDecod[i]=EE1[i][0];	
	    }

	    for(int j=1; j<N/k ; j++){
	       for (int i=j*k ; i<(j+1)*k ; i++){ 
	    		InfoDecod[i]=EE1[i-(j*k)][j];
		   }
	    }
    
    /*printf("\nInfoDecod : \n");
	affiche_vect(InfoDecod,N);
	getch();*/
    
	//return InfoDecod;
	
	free_int_tab(som_syndrome) ;
	free_int_mtx(c,k);
	free_int_mtx(Ap,J*k);
	free_int_mtx(e,n);
	free_int_mtx(E,k);
	free_int_mtx(EE1,k);
	free_int_mtx(R1,n-1);
	free_int_mtx(Info_recuu,n);
	free_int_mtx(u, k );
	 
}
// ---------------------------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------------------------
// Nakagami Fading Channel :
//**********************************************************************************************
long   u[31]={0,119,1522,2314,8939,193,1856,183,419,93,11257,1592,549,25,1254,
    1783,458,3351,4578,65842,2541,5638,289,2583,3274,42,1001,2853,29,303,3125 };

double N0, mu, sigma ;
      //parametre_naka ;

//**********************************************************************************************
double   loi_uniforme(int i)
	{

      long M,c,a;
	    switch(i)
		{
            case 0 : M=2603125;a=596;c=549261;
			break;
			case 1 : M=324000;a=31;c=68363;
			break;
			case 2 : M=196000;a=71;c=41357;
			break;
			case 3 : M=893025;a=106;c=188428;
			break;
			case 4 : M=127008;a=43;c=26797;
			break;
			case 5 : M=209088;a=67;c=44117;
			break;
			case 6 : M=387200;a=111;c=81699;
			break;
			case 7 : M=735075;a=166;c=155101;
			break;
			case 8 : M=205821;a=232;c=43429;
			break;
			case 9 : M=415625;a=666;c=87697;
			break;
			case 10 : M=721875;a=1156;c=152318;
			break;
			case 11 : M=654368;a=287;c=138071;
			break;
			case 12 : M=451737;a=430;c=95317;
			break;
			case 13 : M=284375;a=456;c=60003;
			break;
			case 14 : M=446875;a=716;c=94291;
			break;
			case 15 : M=653125;a=1046;c=137809;
			break;
			case 16 : M=698139;a=664;c=147307;
			break;
			case 17 : M=772497;a=562;c=162997;
			break;
			case 18 : M=614061;a=400;c=129566;
			break;
			case 19 : M=957999;a=400;c=202138;
			break;
			case 20 : M=690625;a=1106;c=145721;
			break;
			case 21 : M=584375;a=936;c=123303;
			break;
			case 22 : M=771875;a=1236;c=162865;
			break;
			case 23 : M=655473;a=274;c=138304;
			break;
			case 24 : M=503125;a=806;c=106159;
			break;
			case 25 : M=790625;a=1266;c=166821;
			break;
			case 26 : M=934375;a=1496;c=197153;
			break;
		    case 27 : M=698625;a=346;c=147409;
			break;
			case 28 : M=676269;a=760;c=142693;
			break;
			case 29 : M=944541;a=898;c=199298;
			break;
			case 30 : M=955719;a=1312;c=201656;
			break;	};

    u[i]=((a*u[i])+c)%M;
    return (double)u[i]/(double)(M);
}

// -------------------- Gaussian  channel-------------------------
double loi_gauss(int loi1,int loi2)
{
	double gauss,v1,v2,vg, f1, f2;
	double u1, u2;
	do
	{
	   u1=loi_uniforme(loi1);
	   u2=loi_uniforme(loi2);
	   v1=(u1+u1)-1.0;
	   v2=(u2+u2)-1.0;
	   gauss=((v1*v1)+(v2*v2));
	}
	while((gauss>=1.0)||(gauss==0.0));

	f1= log(gauss);
	f2=(-f1) + (-f1);
	vg=f2/gauss;

	f2=sqrt(vg);
	gauss=(v1*f2);

	return gauss;
}


 // -------------------- Loi Beta -----------------------
double loi_beta(double pp, double bb)
{
	double s1, s2, beta1, u1, u2;
		do{     u1=loi_uniforme(28);
			u2=loi_uniforme(29);
			s1=pow(u1,1.0/pp);
			s2=pow(u2,1.0/bb);
		   } while ( (s1+s2) > 1.0 );
	    beta1=s1/(s1+s2);
	    return beta1;
}

// loi gamma --------------------------------------------------
// pour alpha=parametre_naka  et  beta=1 : Gamam(alpha,1)----------
double loi_gamma( double parametre_naka)  {
	  double gamma, pp, w, x , v;
	  double y,y1,temp;
	  int kk, v1, continu;
      
	  // m = parametre_naka
      // parametre_naka=0.5;   // k= alpha =parametre_naka  
      
      //v1=2;
	  //pp=1/parametre_naka;
	  
	  v1 = (int) floor(parametre_naka);  // v1 ??
	  pp = parametre_naka-v1;    // pp ??
	  
	  // printf("\n Partie entiere et fractionaire de m = %d,%f",v1,p);
	 etq:	y=1.0;
		for(kk=1;kk<=v1;kk++)
				    {
					temp=loi_uniforme(2+kk);
					y*=1.0-temp;
				    }
		            y1= y;
                    
		if(y<1.0e-30) 
				{
			  	// printf("\nerreur terme y = %e",y);
			  	goto etq;
				}
        
		if (pp==0.0)  	
		      gamma=(double)(-log(y1));
		else
		    {
			  w=loi_beta(pp, 1.0-pp);
// 			  printf("\n terme loi beta w = %e",w);

			  temp= loi_uniforme(30);
			  x=-log(1.0-temp);
//			  printf("\n terme loi exp x = %e",x);
              
			  gamma=(double)(-log(y1))+w*x;
		     }
	return gamma;
}

// Gamma (alpha,beta) :--------------------------------------------------
double loi_gamma_1( double alpha, double beta){
	double X ;
	double Y ;
	X = loi_gamma(alpha) ;   // X suit gamma(alpha,1)
	//Y = pow(beta,-1) * X ;
	Y = (1/beta) * X ;
	return Y ;   // Y suit Gamma(alpha,beta)
}

double Gamma_Gamma(double alpha , double beta){
	double GG;
	double x, y;
	
	x = loi_gamma_1(  alpha, alpha ) ;
	y = loi_gamma_1(  beta, beta ) ;
	GG = x*y ;
	
	return GG;  // Gamma Gamma
}

// Malaga Ia Random variable: ----------------------------------------------------
double Malaga(double alph , double bta, double omega, double b0, double rho){
	// Arguments: alph, bta, Omega, b0, rho
	double phia=0;
    double phib=-pi/2;
    //double omega_prime=omega+2*b0*rho;
    double G; // Ul, Usc, Usg;
    double Ia; // Ia: Malaga Random Variable
    
	G=loi_gamma_1(bta,bta); // G is a gamma random variable with mean 1
    //Ul=sqrt(G)*sqrt(omega)*exp(i*phia);
 //   Ul=sqrt(G)*sqrt(omega)*(cos(phia)+i*sin(phia));
    
    //Usc=sqrt(rho)*sqrt(G)*sqrt(2*b0)*exp(i*phib);
 //   Usc=sqrt(rho)*sqrt(G)*sqrt(2*b0)*(cos(phib)+i*sin(phib));
    
    //Usg=sqrt(1-rho)*(sqrt(b0))*(loi_gauss(1,2)+i*loi_gauss(1,2));
 //   Usg=sqrt(1-rho)*(sqrt(b0))*(loi_gauss(1,2)+i*loi_gauss(1,2));
    
    //Ia=pow(fabs(Ul+Usc+Usg),2)*loi_gamma_1(alph,1/alph);
   Ia= (pow((sqrt(G)*sqrt(omega)*cos(phia) + sqrt(rho)*sqrt(G)*sqrt(2*b0)*cos(phib) + sqrt(1-rho)*sqrt(b0)*loi_gauss(1,2)),2) + pow((sqrt(G)*sqrt(omega)*sin(phia) + sqrt(rho)*sqrt(G)*sqrt(2*b0)*sin(phib) + sqrt(1-rho)*sqrt(b0)*loi_gauss(1,2)),2)) *loi_gamma_1(alph,alph);
  //Ia= (pow((sqrt(G)*sqrt(omega)*cos(phia) + sqrt(rho)*sqrt(G)*sqrt(2*b0)*cos(phib) + sqrt(1-rho)*loi_gauss(1,2)),2) + pow((sqrt(G)*sqrt(omega)*sin(phia) + sqrt(rho)*sqrt(G)*sqrt(2*b0)*sin(phib) + sqrt(1-rho)*loi_gauss(1,2)),2))*exp(2*loi_gauss(1,2));
    
    return Ia;
}
//-------------------------------------------------------------------------------
// Pointing Errors Ip Random Variable : -----------------------------------------
// Loi de Rayleigh :
double loi_rayleigh(int loi1,int loi2, double sigma) // int i, int j = int loi1,int loi2
{
   double a,b,c,b_carre,c_carre;


   b = loi_gauss(loi1,loi1+1);
   b_carre =  b * b;
   c = loi_gauss(loi2,loi2+1);
   c_carre =  c * c;

   a = sqrt(pow(sigma,2))*sqrt(b_carre + c_carre);

   return a;

}

double Pointing_errors( double sigmaray, double A0, double wzeq){//, double rhomod){
	
	double Ip; // Ip: pointing error Random Variable
	double ray;
	
	//double g = 2*b0*(1-rho);
    //double omega_prime = omega + 2 * b0 * rho;
//	double b = sqrt(pi/2)*a/wz;  // b=nu dans la relation pdf du pointing errors 
	//double A0 = pow(erf(b),2);   
	//double wzeq = wz*sqrt(sqrt(pi)*sqrt(A0)/(2*b*exp(pow(-b,2))));
	// epsilon : c'est chssi (ou gamma petite)
//	double chssi = wzeq/(2*sigmaray);   // see gamma aprés l'éq (11) en [Ahmed farid]
	/*printf("\nchssi: =%e\n", chssi);
	getch();*/
    //double mom1 = chssi^2*A0*(g+omega_prime)/(chssi^2+1); 
    
    ray = loi_rayleigh( 1, 2, sigmaray);  // eq(10) in [Ahmed A. Farid] 
    Ip = A0*exp(-2*pow(ray,2)/pow(wzeq,2));    // eq (9) in [Ahmed A. Farid] 
    
    return Ip; 
}   
//-------------------------------------------------------------------

// Malaga channels with pointing errors :----------------------------
double Malaga_Pointing_errors( double Ia ,double Ip  ){
	double I;
	I=Ia*Ip; // I=Ia*Ip*Il , Il=1
	return I;
}
//-------------------------------------------------------------------
 
// Nakagami (alpha,beta) : beta=1----------------------------------------
double loi_nakagami(double parametre_naka) {       
       
	  double naka;
	  naka=sqrt(loi_gamma(parametre_naka)/parametre_naka);
	  return naka;
}

//--Nakagami (alpha,beta) : beta Diff 1------------------------------
double loi_nakagami_1(double parametre_naka,double omega) {       
       
	  double naka;
	  naka=sqrt(loi_gamma_1(parametre_naka,parametre_naka/omega));
	  return naka;
}

// -----  Bruit AWGN ------------------------------------------------

void bruit_Mixed_Naka_Malaga_Pointing_errors( double snr, double r, double Eta, int N, int Info_Mod[], double mot_bruite[], double alph , double bta, double omega, double b0, double rho, double sigmaray, double A0, double wzeq, double G , double N0_s, double EsN0_n1, double EsN0_n2, int Ns , int N_R, double T0, double T1, double eps, double teta, double Ps, double Kappa_D, double Kappa_R, double L,  double delta, double B_R , double Ps_l, double L_l){

	 double ecart_type, ecart_type1, P_R, sigma_E , Ia, Ip; // Naka ;
	 double val, Interf,  sigma_y1, sigma_y2, Z2;;
	 int i, j, k1;
	 
     double VFading;
	 double **H = Alloc_Mat_double(N, N_R);  initialiser_matrix_double(H, N, N_R);  // fading coeffitients hi
	 double **H_inter = Alloc_Mat_double(N, N_R);  initialiser_matrix_double(H_inter, N, N_R);  // fading coeffitients hi
	 double **Z = Alloc_Mat_double(N, N_R);  initialiser_matrix_double(Z, N, N_R);  // Awgn noise ni
	 double **Z1 = Alloc_Mat_double(N, N_R);  initialiser_matrix_double(Z1, N, N_R);  // Awgn noise ni

	 double **mot_bruite0 = Alloc_Mat_double(N, N_R);    initialiser_matrix_double(mot_bruite0, N, N_R);  // yi=hi*x+ni
	 double *ss = (double*)malloc(N*sizeof(double));   initialiser_Vect_double(ss,N);  // les normes ss[i]= sqrt(h1^2+h2^2+...+h_N_R^2)
	 double **W = Alloc_Mat_double(N, N_R);  initialiser_matrix_double(W, N, N_R);  // hi/norme
	 double *Ysum = (double*)malloc(N*sizeof(double));   initialiser_Vect_double(Ysum,N); 
	 
	 double *mot_bruite1 = (double*)malloc(N*sizeof(double));  initialiser_Vect_double(mot_bruite1,N);
     double *P_E = (double*)malloc(N*sizeof(double));  initialiser_Vect_double(P_E,N);
     
     double **Sources = Alloc_Mat_double(N, N_R);  initialiser_matrix_double(Sources, N, N_R);  // fading coeffitients hi
     double **gam_S1R_i = Alloc_Mat_double(N, N_R);  initialiser_matrix_double(gam_S1R_i, N, N_R);  // fading coeffitients hi

     double *Gam = (double*)malloc(N*sizeof(double));  initialiser_Vect_double(Gam,N);

          double parametre_naka[][N_R] = { {3 ,1 ,3 , 1, 1, 2},//, 5, 4}, //},//  // of source, %l=1, it s source, the others: interference
            						       {2 ,2 ,1 , 2, 2, 3}//, 5, 4}  //}//
            						        };//
            						       // {3 ,2 ,1 , 1, 1 ,2},
									        //{3 ,2 ,1 , 1, 1 ,2},
									        //{3 ,2 ,1 , 1, 1 ,2},
											//{3 ,2 ,1 , 1, 1 ,2},
									        //{3 ,2 ,1 , 1, 1 ,2}
									        //};//
									       
									       
									       
									   
                                       //2 1 2 1 1 1 5 4	};  
     double omega_i[][N_R] = {{ 1 ,1 , 2 ,1 , 1 ,1},//, 1, 1},  //},//  // of source
               				  { 2 ,1 , 1, 2 , 2 ,1}//, 1, 1}  // }//
               				   };//
               				  //{	1, 2 , 1, 1 , 1, 2},
							  //{	1, 2 , 1, 1 , 1, 2},
							  //{	1, 2 , 1, 1 , 1, 2},
							 // {	1, 2 , 1, 1 , 1, 2},
							 // {	1, 2 , 1, 1 , 1, 2}
							  // };//
							  
							 
							    
                                //1 2 1 2 2 1 1 1;};  // of interference 1


	 // EsN0_n = Es/N_0, avec Es=1 ,  ecart_type^2 = sigma^2 = N_0/2 ==>
	 // ===> ecart_type = sqrt(N_0/2.0);
     //ecart_type = sqrt(Ps/(2.0*EsN0_n1));    
    
	// Terminal A : Nakagami -->  y_Naka = h1 * Info_Mod + n1  : h1 is a fading amlitude of Nakagami distribution
	 
	 for(i=0;i<N;i++){
	    for(j=0;j<N_R;j++){
	 	 
	 	     H[i][j] = loi_nakagami_1(parametre_naka[0][j], omega_i[0][j]); //sqrt(loi_gamma_1( parametre_naka[0][j],parametre_naka[0][j]/omega_i[0][j]));  //loi_nakagami_1(parametre_naka[0][j], omega_i[0][j]);//  // k1=0
	 	     Z[i][j] = loi_gauss(1,2);   // noise n_i,k
	 	     
	 	     // RHI noise v_i,k^{Rr}
	 	     sigma_y1 = sqrt( L*Ps*pow(H[i][j],2)*pow(Kappa_R,2) ); // sqrt(variance of RHI)
	 	     Z1[i][j] = sigma_y1 * loi_gauss(1,2)  ; //sigma_y^2=L*Ps*(H[i][j])^2 * Kappa_R^2: variance, mu_y: mean=0
	     }
     }
     
      /*printf("\n H :\n");
      affiche_matrix_double(H,N,N_R) ;
      getch(); 
      
      printf("\n Z :\n");
      affiche_matrix_double(Z,N,N_R) ;
      getch(); 
	  
	  printf("\n Z1 :\n");
      affiche_matrix_double(Z1,N,N_R) ;
      getch(); */
     
    for(j=0;j<N_R;j++){
    	for(i=0;i<N;i++){  
    	    Interf=0;
    	    for(k1=1;k1<Ns;k1++){  // statrt from 1: interference , k1=0: c la source
    	    	Interf = Interf + sqrt(Ps_l*L_l) * loi_nakagami_1(parametre_naka[k1][j], omega_i[k1][j]) * Info_Mod[i];  //
    	    	//printf("\nInterf= %e \n", Interf);
			}
	 	    H_inter[i][j] = Interf;
	 	}
	}
	
	 /* printf("\n H_inter :\n");
      affiche_matrix_double(H_inter,N,N_R) ;
      getch(); */
     
     //_____________________Received signal at R ___________________________________
	 for(j=0;j<N_R;j++){
	     for(i=0;i<N;i++){ 
	 	      mot_bruite0[i][j] = sqrt(Ps*L) *  H[i][j] * Info_Mod[i] + H_inter[i][j] + Z1[i][j] + (sqrt(Ps/(2.0*EsN0_n1*omega_i[0][j])) *Z[i][j] ) ;
	     }
     }
      /*printf("\n mot_bruite0 :\n");
      affiche_matrix_double(mot_bruite0,N,N_R) ;
      getch();*/
     
     // _________________________MRC at Relay_____________________________________
	 // les normes 
     for(i=0;i<N;i++){ 
         for(j=0;j<N_R;j++){ 
     	     ss[i] = ss[i] + pow(H[i][j],2);
     	 }
     	 ss[i] = sqrt(ss[i]);   
     }
      
	      /*printf("\n les normes: ss :\n"); 
	      affiche_vect_double(ss,N);	
	      getch();*/
     
    for(i=0;i<N;i++){ 
         for(j=0;j<N_R;j++){
     	      W[i][j] = H[i][j]/ss[i];
     	      Ysum[i] = Ysum[i] + W[i][j] * mot_bruite0[i][j] ;
     	 } 
     } 
          /*printf("\n W :\n");
          affiche_matrix_double(W,N,N_R) ;
          getch();
             
          printf("\n Ysum :\n"); 
	      affiche_vect_double(Ysum,N);	
	      getch(); */
	      
	//__________________________gam_S1R (output of MRC)___________________________
    for(j=0;j<N_R;j++){
    	for(i=0;i<N;i++){ 
    	    for(k1=1;k1<Ns;k1++){
    	    	Sources[i][j] = Sources[i][j] + Ps_l*L_l/N0_s * pow( loi_nakagami_1(parametre_naka[k1][j], omega_i[k1][j]) , 2 );   //loi_gamma_1( parametre_naka[k1][j], parametre_naka[k1][j]/omega_i[k1][j]);    // ;
			}
			//gamma_S1R_i
	 	    gam_S1R_i[i][j] =  Ps*L/N0_s * pow(H[i][j],2) /( Sources[i][j]+  Ps*L/N0_s *pow(H[i][j],2) * pow(Kappa_R,2) + 1 ) ;
	 	}
	 }
	        /*printf("\n gam_S1R_i :\n");
            affiche_matrix_double(gam_S1R_i,N,N_R) ;
            getch(); */
	 
	 double s_S1R = 0 ;
	 double gam_bar_S1R = 0 ;
	 
	 for(i=0;i<N;i++){
	 	for(j=0;j<N_R;j++){
	 		Gam[i] = Gam[i] + gam_S1R_i[i][j] ;
	 	}
	 	s_S1R = s_S1R + Gam[i] ; //printf("\ns_S1R= %e \n", s_S1R);
	 }
	 gam_bar_S1R = s_S1R/N;  //printf("\ngam_bar_S1R= %e \n", gam_bar_S1R); getch();   // Average of sum_i=1^{N_R} gam_S1R_i
	 
	      /*printf("\n Gam :\n"); 
	      affiche_vect_double(Gam,N);	
	      getch(); */

    // _____________________Harvested Energy and Gain G (AF)________________________
	 for(i=0;i<N;i++){ 
	     //P_E = loi_gamma_1( parametre_naka*N_R, sigma_E) ;
	 	 //P_E = ss[i] * (teta*eps*T0*Ps)/pow(d_R,delta);
         for(k1=0;k1<Ns;k1++){ //Ns le nombre des intereferences, +1 : la source
         	 for(j=0;j<N_R;j++){
         	 	if(k1==0){
         	 		P_E[i] = P_E[i] + (teta*eps*T0*Ps*L)/T1 * loi_gamma_1(parametre_naka[k1][j], parametre_naka[k1][j]/omega_i[k1][j] ) ;  }
				else{
					P_E[i] = P_E[i] + (teta*eps*T0*Ps_l*L_l)/T1 * loi_gamma_1(parametre_naka[k1][j], parametre_naka[k1][j]/omega_i[k1][j] ) ;  }
			 }                                           
         }
	      
         if(P_E[i]<B_R){
	 	 	 P_R = P_E[i]/T1;
		 }
	 	 else{
	 	 	 P_R = B_R/T1;
		 } 
		 //P_R = 1 ;
	     mot_bruite1[i] = sqrt(P_R) * Ysum[i]  ;  // G : le gain 
     }
     
     //------------------------------
     /*FILE *fi1;
	 fi1=fopen("P_E.m","at");
	 for(i=0;i<N;i++){ 
         fprintf(fi1,"%e\n",P_E[i]);
	 }*/
     //------------------------------
     
	 // le gain
	 //G = sqrt(E[P_R]/(N0_s*C))) ;
	 double E_PE=0;   // average of P_E
	 for(k1=0;k1<Ns;k1++){
        for(j=0;j<N_R;j++){
        	if(k1==0){
        		E_PE = E_PE + (teta*eps*T0*Ps*L)/T1 *  omega_i[k1][j] ; }
			else{
			    E_PE = E_PE + (teta*eps*T0*Ps_l*L_l)/T1 *  omega_i[k1][j] ; }
		}                                           
     }  
      // E_PE = 1 ; 
      //printf("\nE_PE= %e \n", E_PE);    getch();
	  G = sqrt(E_PE/(N0_s*(gam_bar_S1R + 1))) ;  //printf("\nG= %e \n", G);   // semi blind
	  //G = sqrt(E_PE) * sqrt( 1 /(N0_s*(gam_bar_S1R + 1))) ;    // without EH
	  
	 // Terminal B , Relay -->  G * y_Naka : Amplifiying the signal
	 for(i=0;i<N;i++){
	     mot_bruite1[i] = G * mot_bruite1[i] ;  // G : le gain 
	 }
	      /*printf("\n mot_bruite1 :\n"); 
	      affiche_vect_double(mot_bruite1,N);	
	      getch();*/
	 
     // _______________________ Second hop: Malaga with MRC ________________________
	 ecart_type1 = sqrt(1.0/(2.0*EsN0_n2)); 
	 for(i=0;i<N;i++){
	 	  val = loi_gauss(1,2);
		  
		  Ia = Malaga( alph, bta, omega, b0, rho);
		  Ip = Pointing_errors(    sigmaray,    A0,   wzeq); // Ip: pointing error Random Variable
	 	  VFading = Malaga_Pointing_errors(Ia , Ip) ; 
	 	  
	 	  // RHI noise v_i,k^{Rr}
	 	  sigma_y2 = sqrt( pow(G,2) * P_R * pow(Eta*VFading, r) * pow(Kappa_D,2) ) ; // variance
	 	  Z2 = sigma_y2 * loi_gauss(1,2)  ; //sigma_y^2=L*Ps*(H[i][j])^2 * Kappa_R^2: variance, mu_y: mean=0

	 	  // For r=1 --> Heterodyne Detection, For r=2 --> IM/DD Detection
	 	  mot_bruite[i] =    pow(sqrt(Eta * VFading ),r) * mot_bruite1[i] + Z2 + ecart_type1 * val;     // mot bruité = mot modulé + bruit    
	 } 
	 
	     /*printf("\n mot_bruite :\n"); 
	      affiche_vect_double(mot_bruite,N);	
	      getch();*/
	      
	/*for(i=0;i<N;i++){
		mot_bruite[i] = mot_bruite[i] ; // Ysum[i] ; //mot_bruite0[i][0];
    }*/
	 
	 // free 
	 free_double_tab(mot_bruite1);
	 free_double_tab(Ysum);
	 free_double_mtx(W ,N );
	 free_double_tab(ss);
	 free_double_mtx(mot_bruite0,N );
	 free_double_mtx(Z ,N );
	 free_double_mtx(H ,N );
}
//**********************************************************************************************

//**********************************************************************************************
void generer_bloc( int mot_info[] ,int N)
{
   int i;
   //double val;

   for(i=0;i<N;i++)
   {
   	  if(loi_uniforme(0)>0.5) mot_info[i]=1;
      else mot_info[i]=0;
      //Info[i]=loi_uniforme(0); 
   }
   
   //free_int_tab(mot_info) ;
}
// ---------------------------------------------------------------------------------------------

//Modelise le modulateur, le canal gaussien et le d‚modulateur : HISO

// ------------------------------------------------------------------------------------------
// Déviser une info en bloc :-----------------------------------------------------------------
int **mtx_blocs( int nL,int nC, int *info,int N ){
         	 // N : taille info totale
			 //int nL = 2; // Nbre bloc 
             //int nC = taille/2; // Longueur de chaque bloc
	         int **Mtx_Bloc = Alloc_Mat ( nL, nC);
			 initialiser_matrix(Mtx_Bloc ,nL ,nC);
			 
			 for(int j=0; j < nL ;j++){
			     for(int i = nC-1; i >=0 ;i--){
    		        Mtx_Bloc[j][(nC-1)-i]=info[nC *(j+1)-(i+1)];
				}
			}
		
            return Mtx_Bloc;
          
		  
        free_int_tab(info) ;
		free_int_mtx(Mtx_Bloc, nL ); 
}
//---------------------------------------------------------------------- 

// Ajout des bits de terminaison : --------------------------------------         
 int **Ajout_bit_Term( int **m_bloc,int nL,int nC, int nCm   ){  
            //int nCm=(nC + k*m);
		
			int **m_bloc_bit_term = Alloc_Mat ( nL , nCm ); 
			initialiser_matrix(m_bloc_bit_term ,nL , nCm );
		
			for (int i=0; i < nL; i++){
    	        for(int j=0; j < nC ; j++){
    		        m_bloc_bit_term[i][j]=m_bloc[i][j];
				}
			}
			return m_bloc_bit_term;
			
		free_int_mtx(m_bloc, nL );
		free_int_mtx(m_bloc_bit_term, nL);
			
 } 
 //---------------------------------------------------------------------- 
int bits_differents(int mot_info[],int mot_decoded[],int N)
{
int i,nr=0; 
for(i=0;i<N;i++)
//if(mot_info[i]!=mot_decoded[i+r])
if(mot_info[i]!=mot_decoded[i])
          nr++;
  return nr;
  
  free_int_tab(mot_decoded) ;
  free_int_tab(mot_info) ;
  
} 
  
 
// Function BPSK Modulation  -------------------------------------------------------------------------------------------------
void BPSK_Mod( int *Code_Word , int *Info_Mod,  int Taille_Code_Word){
    
	for (int i=0 ; i<Taille_Code_Word ; i++){
    	if(Code_Word[i]==0){
    		Info_Mod[i]=-1;
		}
		else{// (Code_Word[i]==1)
			 Info_Mod[i]=1;
		}
	}
}//----------------------------------------------------------------------------------------------------------------------------

// Function BPSK Démodulation  -------------------------------------------------------------------------------------------------
void BPSK_Demod( int *Info_Demod, double *mot_bruite , int Taille_mot_bruite){
	
	for (int i=0 ; i<Taille_mot_bruite ; i++){
    	if(mot_bruite[i]>0){
    		Info_Demod[i]=1;
		}
		else{   // (Info_Demod[i]<0)
			Info_Demod[i]=0;
		}
	}
}//-------------------------------------------------------------------------------------------

// Calcul du Temps d'execution du Programme : **************************************
int frequency_of_primes (int nn){
  int i,j;
  int freq=nn-1;
  for (i=2; i<=nn; ++i) for(j=(int)sqrt(i);j>1;--j) if (i%j==0) {--freq; break;}
  return freq;
}
// *********************************************************************************        

// ------------*******------------fonction main --------------********-----------------------
int main(){
	
	//*****************************
	clock_t tt;
    tt = clock();
	//*****************************
	
	 // ********************** Equa_Orth ; **********************
	int **A=Alloc_Mat(k*J, k*J+1);
	initialiser_matrix(A,k*J, k*J+1);  
     
	 for(int i=0;i<k*J;i++){
    	for(int j=0;j<k*J+1;j++){
    		A[i][j]=Ap[i][j];
		}
     }
	// affiche_matrix(A,k*J,k*J+1)  ;  
	 
	// ***************** Les polynomes générateurs **************					
    int **gg=Alloc_Mat(n,lg);
	initialiser_matrix(gg,n,lg);  
 
	 for(int i=0;i<n;i++){
    	for(int j=0;j<lg;j++){
    		gg[i][j]=pg[i][j];
		}
     }
     
     // affiche_matrix(gg,n,lg);  
     // -----------------------------------------------------------------------------
     
     int NN=500 ;   // il faut choisir NN superieur à lg
	 int N=k*NN, N_bt=N+k*m,  dimInfoDemo=n*(NN+m) ; //
	 //printf("\n NN=%d , N=%d , N_bt =%d ,k=%d, dimInfoDemo=%d\n",NN,N,N_bt,k,dimInfoDemo);
	 //printf("\n k*NN+k*m = %d*%d+%d*%d \n",k,NN , k ,m);
      
	 //Parametres -------------------------------------------------------------------
	 int *mot_info = (int*) malloc( N * sizeof(int) );                        initialiser_Vect(mot_info,N); 
	 int *mot_info_bt = (int*) malloc(N_bt * sizeof(int) );                   initialiser_Vect(mot_info_bt , N_bt);
	 int *cod_Out = (int*) malloc( dimInfoDemo * sizeof(int) );               initialiser_Vect(cod_Out,dimInfoDemo);
	 int *InfoMod = (int*) malloc( dimInfoDemo * sizeof(int) );               initialiser_Vect(InfoMod,dimInfoDemo);	 
	 double *mot_bruite = (double*) malloc (dimInfoDemo * sizeof(double));    initialiser_Vect_double(mot_bruite,dimInfoDemo);
	 int *InfoDemo = (int*) malloc( dimInfoDemo * sizeof(int) );              initialiser_Vect(InfoDemo,dimInfoDemo);
	 int *InfoDecode = (int*) malloc( N * sizeof(int) );                      initialiser_Vect(InfoDecode,N);
	 
	 /*printf("init mot_info_bt : \n");	
	 affiche_vect(mot_info_bt,N_bt);*/
	 
	 // --------------------------------------------------------
	 long MinMot =2000;   // nbre minimale de mot (bloc)
	 long MinErreur=500;  // nbre minimale d'erreur résiduelle
     // --------------------------------------------------------
     
	 double snr , c ;
	 double rend_code,EsN0_n,EbN0_n, BER;
	 long mot, nber_erreur, nber_total; 
	 double EsN0_n1, EsN0_n2;
	 
	 // channel parameters :  _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 
	 // Nakagami:_____________________________________________________ 
     int N_R=6;// Relay Antennas
     int Ns=2; // sources nomber
     /*double parametre_naka[Ns][N_R] = { {3 ,1 ,3 ,1 },  // of source, %l=1, it s source, the others: interference
            							{2 ,2 ,1, 2 }
										};  
     double omega_i[Ns][N_R] = {{ 1 ,1, 2 ,1 },  // of source
               				 	{ 2 ,1 ,1, 2 }
								};  // of interference 1
     */
	           			
	 double Kappa_R = 0.2; // RHI
	 double Kappa_D = 0.1 ;
	 
	double Ps = 10 ;  //20
	double delta = 2.25 ; 
	double puis_G_s = 4 ;  // 30/10   // printf("\n puis =%e \n ",puis);  //double(21/10) 
	double puis_G_R = 2.5 ;     // printf("\n puis =%e \n ",puis);  //double(21/10) 
    double G_s = pow(10,puis_G_s); 
    double G_R = pow(10,puis_G_R); //21 dB
    double d_0 = 100;
    double d_sLr = 60;
    double Lamda = 28.6*pow(10,-3) ;  //=c/f
    double L = G_s*G_R* pow(Lamda/(4*pi*d_0),2) * pow(d_0/d_sLr,delta) ;   //printf("\n L =%e \n ",L);  getch();

    double Ps_l = 5 ; 
    double G_sl= pow(10,1) ;  //10/10
    double d_sl = 120 ;
    double Lamda_l=  28.6*pow(10,-3) ;  //printf("\n Lamda_l =%e \n ",Lamda_l);  getch();   //=c/f  
    double L_l = G_sl*G_R* pow(Lamda_l/(4*pi*d_0),2) * pow(d_0/d_sl,delta) ;   //printf("\n L_l =%e \n ",L_l);  getch();
    
     // Malaga :_____________________________________________________
     //double alph=2.296 ,  bta=2;  // Strong
	/* double alph=4.2 ,  bta=3;   // Moderate  
     double omega = 1.3265 ; // 1.3265(Malaga) , 0.78(GG case)  
	 double b0 = 0.1079 ; // 0.1079 (Malaga),  0.0000001(Log_Normal)   
	 double rho = 0.596 ;    // rho = 0.596(Malaga); rho=0.99(GG); rho=0 (LogNormal)  
     */
	 // GG case :____________________________________________________
     double alph=8,  bta=4;
     double omega = 0.78 ; // 1.3265(Malaga) , 0.78(GG case)  
	 double b0 = 0.1079 ; // 0.1079 (Malaga),  0.0000001(Log_Normal)   
	 double rho = 0.99 ;    // rho = 0.596(Malaga); rho=0.99(GG); rho=0 (LogNormal) 
    
	 // Pointing errors : _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
     // epsilon = 1 :________________________________________________
	 //double sigmaray=0.5 ,  wz=1,  a=0.20 ;       //is the Optical aperture radius// see values in table III of article FArid, dossier Malaga
	 //epsilon = 6.7 :______________________________________________
	 double sigmaray=0.5 ,  wz=6.7,  a=0.3 ;       //is the Optical aperture radius// see values in table III of article FArid, dossier Malaga
    
     // Energy Harvesting :__________________________________________
     double T0 = 1 ; 
	 double T1 = 1 ;
	 double eps = 0.7 ;
	 double teta = 0.7 ;
	 double B_R = 500;
	 //________________________________________________________________	 
	 // Detection technique : _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 
	 double r = 1;   // Heterodyne Detection
	 //double r = 2;   // IM/DD Detection
	 // _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 
	 double Eta = 0.9;
	
	 // epsilon du pointing errors: _________________________________
	 double b = sqrt(pi/2)*a/wz;  //   // c'est 'nu'
	 double A0 = pow(erf(b),2); 
	 double wzeq = wz*sqrt(sqrt(pi)*sqrt(A0)/(2*b*exp(-pow(b,2))));  
	 double epsilon = wzeq/(2*sigmaray);  // \xi
    
	 // Le 1er Moment E[I], 2ième Moment E[I^2] _______________________
	 double g = 2*b0*(1-rho);   
     double omega_prime = omega+2*b0*rho; 
     double mom; 
	 if(r==1)
     	mom = pow(epsilon,2)*A0*(g+omega_prime)/(pow(epsilon,2)+1); //mom1 : E[I] et I=Ia*Ip
     else // r==2 
		mom = (alph+1)/alph*(pow(omega_prime,2)*(1+1/bta)+4*omega_prime*g+2*pow(g,2))*pow(epsilon,2)*pow(A0,2)/(pow(epsilon,2)+2);  //mom2 : E[I^2]
     //________________________________________________________________
       
	   double N0_s, N_0; 
	   double G;   // le gain 
	 
	double snr_i = 0 ;    
    double pas = 1 ;
    double snr_max = 16 ;
	 //----------------------------------------------------------------
	 FILE *fi;
	 fi = fopen("BER_CSOC_RF_FSO_RHI_Results_Journal2.txt","at");
	 fprintf(fi,"\n\n__________________________________BER Results for the code CSOC(%d,%d,%d)__________________________________\n", n,k,m);
     fprintf(fi,"N_R = %d, Ns = %d, alph=%e,  bta=%e, r=%e , epsilon=%e\n",N_R,Ns, alph,bta,r,epsilon);
     fprintf(fi,"L=%e, L_l=%e, eps=%e, Ps=%e, Ps_l=%e, d_sLr=%e, d_sl=%e, Kappa_R=%e, Kappa_D=%e, puis_G_s=%e\n", L, L_l, eps, Ps,Ps_l, d_sLr,d_sl, Kappa_R, Kappa_D, puis_G_s);
     fprintf(fi,"NN=%d, MinMot=%ld, MinErreur=%ld\n", NN, MinMot, MinErreur);
     fprintf(fi,"SNR=%e:%e:%e;\n",snr_i,pas,snr_max);
	 fprintf(fi,"BER=[");
   	 // ---------------------------------------------------------------
       
   for (snr=snr_i; snr<=snr_max; snr+=pas){
   
        N0_s = Ps*pow(10,-snr/10);
        EbN0_n =  Ps/N0_s  ;  // N0(first hop) = pow(10,-snr/10) //exp((snr*log(10.0))/10.0) ;
       //EbN0_n = Ps/(pow(d_R,delta))*pow(10,snr/10); // en considérant d^delta dans la relation ri(the received information at Ri)
		
		rend_code=(double)k/(double)n ;
        
        // For a first hop _ _ _ _ _ _ _ _ _ _ _ _ _ _
	    EsN0_n1 = rend_code*EbN0_n;
	    
        // For a Relay _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 
        //	  G = sqrt( rend_code *(teta*eps*T0*Ps*N_R*OmegaE)/(pow(d_R,delta)*T1)/(Ps*pow(10,-snr/10))) ;
	    //  C=Es/(G^2 * N_0), C=1 ==> G^2=Es/N_0
        
		// For a second hop_ _ _ _ _ _ _ _ _ _ _ _ _ _ 
        N_0 = mom*pow(Eta,r)*(pow(10,-snr/10));    // N_0: c'est N_1 dans l'article journal
		EsN0_n2 = rend_code *1/N_0 ;              // 1/N_0 ;     Ps=Es: transmit power
        
        // _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 
        
   	    c=0 ;
	    nber_total=0 ;  // nber_total : nbere totale d'erreur dans chaque snr  , nbere erreur : nbre erreur dans chaque bloc
	    mot=0 ;         //nbre de mot dans chaque bloc
	    
        while(c==0){    // boucle infini car le nbre d'erreur est toujours égale à 0
		     
			 // **************** Info ********************
	          //source_info(mot_info, N);
	        generer_bloc(  mot_info , N );
	    /*printf("\nmot_info : \n");
	    affiche_vect(mot_info,N);
	    getch();*/
	    
	         // ajout bit de terminaison ------------------------
	        for(int i=0;i<N;i++){ mot_info_bt[i]=mot_info[i] ; }
		/*printf("\nmot_info_bt : \n");
	    affiche_vect(mot_info_bt,N_bt);
	    getch();*/
	         //--------------------------------------------------
	          
			 // *************** Codage : ****************** 
	          //cod_Out = code_conv_1_n(mot_info_bt,gg, N_bt);
	         codage_conv(mot_info_bt, cod_Out, gg, N_bt);
	    /*printf("\ncod_Out : \n");
	    affiche_vect(cod_Out,dimInfoDemo);
	    getch(); */
	          
			 // **************** Modulation ****************
		     BPSK_Mod(cod_Out, InfoMod, dimInfoDemo);
		/*printf("\nInfo_mod : \n");
	    affiche_vect(InfoMod,dimInfoDemo);
	    getch();*/
	    
			 // **************** Canal awgn: ***************
			 // bruit_Gaussien(  dimInfoDemo, Info_mod, mot_bruite,  EsN0_n );
			 //bruit_Rayleigh( dimInfoDemo, Info_mod, mot_bruite, EsN0_n);
			 //bruit_Nakagami( parametre_naka,dimInfoDemo, Info_mod, mot_bruite, EsN0_n);
	         bruit_Mixed_Naka_Malaga_Pointing_errors(snr, r , Eta , dimInfoDemo, InfoMod, mot_bruite, alph, bta, omega, b0, rho , sigmaray, A0, wzeq , G ,N0_s, EsN0_n1, EsN0_n2, Ns,N_R, T0, T1, eps, teta, Ps, Kappa_D, Kappa_R, L, delta, B_R, Ps_l, L_l); 
		    
		/*printf("\nMot bruité : \n");
	    affiche_vect_double(mot_bruite,dimInfoDemo);
	    getch();*/
	    
	         // ************** Démodulation :**************
	         BPSK_Demod( InfoDemo, mot_bruite , dimInfoDemo);
	         //InfoDemo = BPSK_Demod( Info_mod , dimInfoDemo);
	    /*printf("\nInfoDemo : \n");
	    affiche_vect(InfoDemo,dimInfoDemo);
	    getch(); */
	         // ************** Décodage : *****************
		     Decoder_MLD_HIHO(InfoDecode, InfoDemo, N,  A , gg);
	    /*printf("\nInfoDecode : \n");
	    affiche_vect(InfoDecode,N);
		getch();*/
		
	         // ************** Calcul d'erreur: ***********
	         nber_erreur = bits_differents(mot_info, InfoDecode, N);
	         //nber_erreur = bits_differents(cod_Out, InfoDemo, dimInfoDemo);
             nber_total = nber_total + nber_erreur;
             
             mot = mot + 1; 
             
        //printf("\nmot=%ld , MinMot=%ld , MinErreur=%ld \n",mot,MinMot,MinErreur);
        //printf("\nsnr=%d, mot=%ld , nber_erreur=%ld,  nber_total=%ld \n", snr, mot ,nber_erreur,nber_total);
        //getch();
        
             if(mot>=MinMot){
             	if(nber_total>=MinErreur){
                   	c++; 
				} 
			 }
			 			 
        } 
        
      BER=(double)nber_total/((double)mot*(double)N);  
      printf("\nsnr = %e ,mot=%ld , nber_total=%ld  << BER===%e \n",snr, mot ,nber_total, BER);	
	   //fprintf(fi,"\nsnr==%d,  mot =%ld , nber_total==%ld ,  BER===%e \n",snr, mot ,nber_total,BER);
      fprintf(fi," %e \t",BER);
	  	//getch(); 	   
   } 
   
    fprintf(fi," ]\n");
		
	// libérer la mémoire des tableaux alloués : 
      //free_double_tab(BER) ;
      free_int_tab(InfoDecode) ;
      free_int_tab(InfoDemo);
      free_double_tab(mot_bruite);
      free_int_tab(InfoMod);
	  free_int_tab(cod_Out );
	  free_int_tab(mot_info_bt); 
	  free_int_tab(mot_info);
	  free_int_mtx(A, k*J );
	  free_int_mtx(gg, n );
	
	//***************************************
	tt = clock() - tt;
    printf ("It took me %d clicks (%f seconds).\n",tt,((float)tt)/CLOCKS_PER_SEC);
    fprintf (fi,"It took me  (%f seconds).\n",((float)tt)/CLOCKS_PER_SEC);
    //***************************************
	   	
	getch();
	return 0;
} 
