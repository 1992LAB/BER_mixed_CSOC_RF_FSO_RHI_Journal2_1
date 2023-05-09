 
 
// -----------------------------------------------------------------------------------------------
/*le fichier.h du code(3,2,13)*/

#define k 2        // Le nbre d'entrés du code
#define n 3        // Le nbre de sorties du code
#define lg 14      // La longueur du grand polynôme générateur
#define m lg-1     // La mémoire du code
#define J 4        // Le nbre des équations de parité orthogonales // le nbre des entiers positifs
#define tML J/2    // la capacité de correction du code

// Les polynomes générateurs du code: // lg=14;
int pg[][lg]={{1   ,  0  ,   0  ,   0   ,  0  ,   0   ,  0   ,   0   ,  0   ,   0  ,   0  ,  0  ,   0  ,   0},
              {1   ,  0  ,   0  ,   0   ,  0  ,   0   ,  1   ,   0   ,  0   ,   0  ,   0  ,  1  ,   0  ,   1},
              {1   ,  0  ,   0  ,   0   ,  0  ,   0   ,  0   ,   0   ,  1   ,   1  ,   0  ,  0  ,   1  ,   0},
			  }; 

// Les Equations de parité orthogonales :   // k*J+1=9, // k*J=8 
int Ap[][k*J+1]={{1  , 0 ,  0 ,  0 ,  1 ,  0 ,  0 ,  0 ,  1 },
                 {1  , 7 ,  0 ,  0 ,  7 ,  0 ,  0 ,  0 ,  7 },
                 {1  , 6 , 12 ,  0 ,  3 ,  4 , 12 ,  0 , 12 },
                 {1  , 3 ,  8 , 14 ,  2 ,  5 ,  6 , 14 , 14 },
                 {1  , 0 ,  0 ,  0 ,  1 ,  0 ,  0 ,  0 ,  1 },
                 {3  , 9 ,  0 ,  0 ,  1 ,  9 ,  0 ,  0 ,  9 },
                 {4  , 10,  0 ,  0 ,  1 ,  2 , 10 ,  0 , 10 },
                 {2  , 7 , 13 ,  0 ,  1 ,  4 ,  5 , 13 , 13 }
                };  
                
  
