# COMNUMFIP

![GitHub release](https://img.shields.io/github/v/release/vincmazet/comnumfip)


**Comnumfip** est un module Python pour les TP de communications numériques
à [Télécom Physique Strasbourg](http://www.telecom-physique.fr/).
Il est associé au cours de [Communications numériques](https://vincmazet.github.io/comnum/).

Il est associé aux fichiers nécessaires pour effectuer les TP de ce cours.


## Installation

Téléchargez le code en cliquant sur le bouton _Code_ ci-dessus, puis _Download ZIP_.
Décompressez les fichiers dans votre dossier de travail.


## Documentation


### eyediag

`eyediag(t, x, T, alpha=.5, color="tab:blue")`

Diagramme de l'oeil.

Entrées :
* **t** (array) : temps
* **x** (array) : signal
* **T** (scalar) : durée d'un symbole
* **alpha** (scalar) : transparence (0,5 par défaut)

Sortie :
* aucune


### rrc

`y = rrc(t, V, a)`
    
Impulsion en racine de cosinus surélevé (_root raised cosine_),
pour une durée de symbole égale à 1.

Entrées :
* **t** (array)  : temps
* **V** (scalar) : amplitude de l'impulsion
* **a** (scalar) : facteur de retombée (_roll-off factor_)
    
Sortie :
* **y** (array) : impulsion en racine de cosinus surélevé


### randmary

`c = randmary(N,p)`

Génération d'une séquence M-aire.

Entrées :
* **N** (scalar) : taille de la séquence (nombre de symboles)
* **P** (array)   : probabilité des symboles (sa taille correspond à la taille de l'alphabet)

Sortie :
* **c** (array) : séquence aléatoire M-aire où M = len(P). Le LSB est à gauche.

Exemples :

```
# séquence binaire de taille 1000, symboles équiprobables :
c1 = randmary(1000,[0.5, 0.5])

# séquence binaire de taille 100, p("0") = 0.3, p("1") = 0.7 :
c2 = randmary(100,[0.3, 0.7])

# séquence 4-aire de taille 10, symboles équiprobables :
c3 = randmary(10,np.ones(4)/4)
```


### bin2mary

`y = bin2mary(x,M)`

Convertit une séquence binaire en séquence M-aire.
Si la taille de x n'est pas multiple de log2(M), des "0" sont rajoutés à la fin de x.

Entrées :
* **x** (array)  : séquence binaire
* **M** (scalar) : taille de l'alphabet de la séquence traduite (M est une puissance de 2)
    
Sortie :
* **y** (array) : séquence M-aire



### mod_a, mob_b, mod_c, mod_d, mod_e

`t, x = mod_a(m, V, d)`

`t, x = mod_b(m, V, d)`

`t, x = mod_c(m, V, d)`

`t, x = mod_d(m, V, d)`

`t, x = mod_e(m, V, d)`

Modulations mystères A, B, C, D et E.

Entrées :
* **m** (array)    : séquence binaire (hexadécimale `mod_e`)
* **V** (scalaire) : amplitude de la forme d'onde
* **d** (scalaire) : durée de la forme d'onde
    
Sorties :
* **t** (array) : vecteur temps
* **x** (array) : signal modulé


### mod_rrc

`t, x = mod_rrc(m, V, T, a)`

Modulation NRZ avec une forme d'onde en racine de cosinus surélevé.

Entrées :
* **m** (array)  : séquence binaire
* **V** (scalar) : amplitude de la forme d'onde
* **T** (scalar) : durée de la forme d'onde
* **a** (scalar) : coefficient de retombée (_roll-off factor_)
    
Sorties :
* **t** (array) : vecteur temps
* **x** (array) : signal modulé


### channel

`y = channel(x,fc,s,T)`
    
Simule un canal de transmission en renvoyant le signal _y = x*g + b_ en sortie du canal,
où _g_ est le réponse impulsionnelle d'un filtre passe-bas, _b_ un bruit blanc gaussien et _*_ représente la convolution.

Entrées :
* **x** (array)   : signal émis
* **fc** (scalar) : fréquence de coupure du filtre g
* **s** (scalar)  : écart-type du bruit b
* **T** (scalar)  : durée d'un bit

Sortie :
* **y** (array) : signal transmis via le canal


### rleenc

`code = rleenc(msg)`

Compression RLE (_run length encoding_).

Entrée :
* **msg** (array) : séquence de symboles à compresser

Sortie :
* **code** (array) : séquence compressée en RLE

Exemple :

```
from skimage.io import imread
img = ioimread("image.png")      # Charge l'image image.png
code = rleenc(img.ravel())       # .ravel() permet de vectoriser l'image
                                 # pour en faire un tableau à une seule dimension
```

### rledec

`msg = rledec(code)`

Décompression RLE (_run length encoding_).

Entrée :
* **code** (array) : séquence compressée en RLE

Sortie :
* **msg** (array) : séquence de symboles décompressée

Exemple :

```
from numpy import reshape
msg = rledec(code)               # Effectue la décompression RLE
img = reshape(msg, (M,N))        # Si c'est une image qui est attendue,
                                 # transforme la séquence msg en image de taille M×N
```

### sample_and_threshold

`y = sample_and_threshold(x, T, S)`

Échantillonne à la période T et compare au seuil S le signal x, pour retourner une séquence binaire

Entrées :
* **x** (array)  : signal
* **T** (scalar) : période d'échantillonnage (= durée d'un bit)
* **S** (scalar) : seuil à appliquer

Sortie :
* **y** (array) : séquence binaire


## Licence

Ce programme est distribué sous licence CeCILL-B (www.cecill.info).
Copyright Université de Strasbourg 2021-2022.
Contributeur : vincent.mazet@unistra.fr.
