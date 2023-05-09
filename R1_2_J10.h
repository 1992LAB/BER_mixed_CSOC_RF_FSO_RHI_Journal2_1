
//********Equation Orthogonales du code (55,1,2)********
/*le fichier.h du code(2,1,55)*/

#define k 1      // Le nbre d'entres du code
#define n 2      // Le nbre de sorties du code
#define lg 56     // La longueur du grand polynome generateur
#define m lg-1   // La memoire du code
#define J 10      // Le nbre des equations de parite orthogonales // le nbre des entiers positifs
#define tML J/2  // la capacite de correction du code

//#define dim n*(m+lg)    // la longueur du mot de code 
                       // avec L la longueur de l'information à coder (Systématique)
                       
// Les polynomes generateurs du code: 
int pg[][lg]={{  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,},
			  {  1 ,  0 ,  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  1 ,  0 ,  0 ,  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  1 ,  0 ,  0 ,  0 ,  1 ,  0 ,  0 ,  0 ,  0 ,  1 ,  1 ,},
			 };//0         2                                                          14                                 21                                      29             32                                                               45                  49                       54   55 

                 
// Les Equations de parite orthogonales :   // lignes : k*J+1, // colonne : k*J=8  
/*int Ap[][k*J+1]={{  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  1 ,},
				 {  1 ,  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  3 ,},
				 {  1 , 13 , 15 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 , 22 ,},
				 {  1 ,  9 , 16 , 28 , 30 ,  0 ,  0 ,  0 ,  0 ,  0 , 30 ,},
				 {  1 ,  4 , 12 , 19 , 31 , 33 ,  0 ,  0 ,  0 ,  0 , 33 ,},
				 {  1 , 14 , 17 , 25 , 32 , 44 , 46 ,  0 ,  0 ,  0 , 46 ,},
				 {  1 ,  5 , 18 , 21 , 29 , 36 , 48 , 50 ,  0 ,  0 , 50 ,},
				 {  1 ,  6 , 10 , 23 , 26 , 34 , 41 ,  0 ,  0 ,  0 , 15 ,},
				 {  1 ,  8 , 20 , 22 ,  0 ,  0 ,  0 , 53 , 55 ,  0 , 55 ,},
				 {  1 ,  2 ,  7 , 11 , 24 , 27 , 35 , 42 , 54 , 56 , 56 ,},
				 };*/
				 
				 
int Ap[][k*J+1]={{ 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1},
                 { 1 , 3 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 3},
                 { 1 , 13 , 15, 0 , 0 , 0 , 0 , 0 , 0 , 0 , 15},
                 { 1 , 8 , 20 , 22 , 0 , 0 , 0 , 0 , 0 , 0 , 22},
                 { 1 , 9 , 16 , 28 , 30 , 0 , 0 , 0 , 0 , 0 , 30},
                 { 1 , 4 , 12 , 19 , 31 , 33 , 0 , 0 , 0  ,0 , 33},
                 { 1 , 14 , 17 , 25 , 32 , 44 , 46 , 0 , 0 , 0 , 46},
                 { 1 , 5 , 18 , 21 , 29 , 36 , 48 , 50 , 0 , 0 , 50},
                 { 1 , 6 , 10 , 23 , 26 , 34 , 41 , 53 , 55 , 0 , 55},
                 { 1 , 2 , 7 , 11 , 24 , 27 , 35 , 42 , 54 , 56 , 56 },
                 };
