
//********************Equation Orthogonales du code (40,2,3)********************
/*le fichier.h du code(3,2,40)*/

#define k 2      // Le nbre d'entres du code
#define n 3      // Le nbre de sorties du code
#define lg 41     // La longueur du grand polynome generateur
#define m lg-1   // La memoire du code
#define J 6      // Le nbre des equations de parite orthogonales // le nbre des entiers positifs
#define tML J/2  // la capacite de correction du code


// Les polynomes generateurs du code: 
int pg[][lg]={{  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,},
			  {  1 ,  0 ,  1 ,  0 ,  0 ,  0 ,  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 , 1  ,  0 ,  0 ,  0 ,  0 , 1  ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 , 1 ,},
			  {  1 ,  0 ,  0 ,  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 , 1  ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 , 1  ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 , 1  , 1  ,  0 ,  0 ,  0 ,  0 ,},
			  };


// Les Equations de parite orthogonales :   // lignes : k*J+1, // colonne : k*J=8  
int Ap[][k*J+1]={{  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  1 },
				 {  1 ,  3 ,  0 ,  0 ,  0 ,  0 ,  3 ,  0 ,  0 ,  0 ,  0 ,  0 ,  3 },
				 {  1 ,  5 ,  7 ,  0 ,  0 ,  0 ,  4 ,  7 ,  0 ,  0 ,  0 ,  0 ,  7 },
				 {  1 , 19 , 23 , 25 ,  0 ,  0 , 10 , 22 , 25 ,  0 ,  0 ,  0 , 25 },
				 {  1 ,  6 , 24 , 28 , 30 ,  0 ,  2 , 15 , 27 , 30 ,  0 ,  0 , 30 },
				 {  1 , 12 , 17 , 35 , 39 , 41 ,  5 ,  6 , 13 , 26 , 38 , 41 , 41 },
				 {  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  1 },
				 {  2 ,  4 ,  0 ,  0 ,  0 ,  0 ,  1 ,  4 ,  0 ,  0 ,  0 ,  0 ,  4 },
				 { 10 , 14 , 16 ,  0 ,  0 ,  0 ,  1 , 13 , 16 ,  0 ,  0 ,  0 , 16 },
				 {  5 , 23 , 27 , 29 ,  0 ,  0 ,  1 , 14 , 26 , 29 ,  0 ,  0 , 29 },
				 {  7 , 12 , 30 , 34 , 36 ,  0 ,  1 ,  8 , 21 , 33 , 36 ,  0 , 36 },
				 {  8 , 13 , 31 , 35 , 37 ,  0 ,  1 ,  2 ,  9 , 22 , 34 , 37 , 37 }
				 };
