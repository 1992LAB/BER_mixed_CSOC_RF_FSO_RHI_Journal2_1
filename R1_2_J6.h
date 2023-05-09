//********Equation Orthogonales du code (17,1,2)********
/*le fichier.h du code(2,1,17)*/

#define k 1      // Le nbre d'entres du code
#define n 2      // Le nbre de sorties du code
#define lg 18     // La longueur du grand polynome generateur
#define m lg-1   // La memoire du code
#define J 6      // Le nbre des equations de parite orthogonales // le nbre des entiers positifs
#define tML J/2  // la capacite de correction du code


//#define dim n*(m+lg)    // la longueur du mot de code 
                       // avec L la longueur de l'information à coder (Systématique)

// Les polynomes generateurs du code: 
int pg[][lg]={{  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 , 0 ,  0 ,  0 ,  0 ,  0 },
			  {  1 ,  0 ,  1 ,  0 ,  0 ,  0 ,  0 ,  1 ,  0 ,  0 ,  0 ,  0 ,  0 , 1 ,  0 ,  0 ,  1 ,  1 }
			  };


// Les Equations de parite orthogonales :   // lignes : k*J+1, // colonne : k*J=8  
int Ap[][k*J+1]={{  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  1 },
				 {  1 ,  3 ,  0 ,  0 ,  0 ,  0 ,  3 },
				 {  1 ,  6 ,  8 ,  0 ,  0 ,  0 ,  8 },
				 {  1 ,  7 , 12 , 14 ,  0 ,  0 , 14 },
				 {  1 ,  4 , 10 , 15 , 17 ,  0 , 17 },
				 {  1 ,  2 ,  5 , 11 , 16 , 18 , 18 },
				 };

