
	
//********Equation Orthogonales du code (35,1,2)********
/*le fichier.h du code(2,1,35)*/

#define k 1      // Le nbre d'entres du code
#define n 2      // Le nbre de sorties du code
#define lg 36     // La longueur du grand polynome generateur
#define m lg-1   // La memoire du code
#define J 8      // Le nbre des equations de parite orthogonales // le nbre des entiers positifs
#define tML J/2  // la capacite de correction du code

//#define dim n*(m+lg)    // la longueur du mot de code 
                       // avec L la longueur de l'information à coder (Systématique)
                       
// Les polynomes generateurs du code: 
int pg[][lg]={{  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,},
			  {  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  1 ,  0 ,  0 ,  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  1 ,  0 ,  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  1 ,  1 ,  0 ,  0 ,  0 ,  1 ,},
			  };


// Les Equations de parite orthogonales :   // lignes : k*J+1, // colonne : k*J=8  
int Ap[][k*J+1]={{  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  1 ,},
				 {  1 ,  8 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  8 ,},
				 {  1 ,  4 , 11 ,  0 ,  0 ,  0 ,  0 ,  0 , 11 ,},
				 {  1 ,  7 , 10 , 17 ,  0 ,  0 ,  0 ,  0 , 17 ,},
				 {  1 ,  3 ,  9 , 12 , 19 ,  0 ,  0 ,  0 , 19 ,},
				 {  1 , 13 , 15 , 21 , 24 , 31 ,  0 ,  0 , 31 ,},
				 {  1 ,  2 , 14 , 16 , 22 , 25 , 32 ,  0 , 32 ,},
				 {  1 ,  5 ,  6 , 18 , 20 , 26 , 29 , 36 , 36 ,},
				 };		
