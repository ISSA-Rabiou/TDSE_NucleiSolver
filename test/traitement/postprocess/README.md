# Post-traitement d’une fonction d’autocorrélation

Ce dossier fournit une version corrigée et autonome de `spectre.f90`.

## Formats d’entrée

`auto_cor.txt` accepte une ligne numérique par point :

- `t Re(C)` : autocorrélation réelle ;
- `t Re(C) Im(C)` : autocorrélation complexe ;
- `label t Re(C) Im(C)` : ancien format à quatre colonnes.

Les lignes vides et les commentaires commençant par `#` ou `!` sont ignorés. La grille en temps doit être strictement croissante et régulière.

## Compilation et exécution

```bash
./script.sh
./bin/spectre < param.nml
```

Le résultat est écrit dans `file_spectrum.txt`, avec une ligne d’en-tête commentée et les colonnes :

```text
E  Re[S0(E)]  Re[S_cos(E)]  Re[S_cos2(E)]  Re[S_window(E)]
```

Le namelist `param.nml` permet de régler `dE`, `Emin`, `Emax`, les facteurs d’unité `conv_t`/`conv_E`, le fichier d’entrée/sortie et l’ordre du filtre exponentiel `n_filter`. Une valeur `Emax=-1` choisit automatiquement la limite de Nyquist ; `dE=-1` choisit un pas spectral automatique (résolution naturelle, sans grille démesurée).

## Tracé

Ouvrir `trace_autocorrelation.ipynb` dans Jupyter depuis ce dossier. Le notebook trace la partie réelle et imaginaire de l’autocorrélation, son module et les quatre spectres produits par Fortran.

## Corrections principales

- suppression des lectures fragiles par mot et de l’option implicite non documentée ;
- détection robuste des formats 2/3/4 colonnes ;
- validation du nombre de points, de l’ordre et de la régularité des temps ;
- protection contre les accès hors limites (`t(2)` si un seul point) ;
- correction de la quadrature : poids trapézoïdaux aux deux extrémités ;
- contrôle de la limite de Nyquist et calcul correct du nombre de points spectraux ;
- sortie documentée et remplacement sûr du fichier existant ;
- commentaires et messages d’erreur explicites.

## Détection des pics dans le notebook

Dans `trace_autocorrelation.ipynb`, modifier les paramètres :

```python
DATA_FILE = Path('/chemin/vers/auto_cor.txt')
SPECTRUM_FILE = Path('/chemin/vers/file_spectrum.txt')
SPECTRUM_COLUMN = 1  # colonne I(E) à analyser
N_PEAKS = 10         # nombre maximal de pics
I_MIN = 0.0          # seuil minimal
```

Le notebook parcourt les maxima locaux `I[i] >= I[i-1]` et `I[i] >= I[i+1]`. Pour chaque maximum admissible, il utilise exactement les trois couples autour du point central et calcule les coefficients `a`, `b`, `c` de `a E²+b E+c`. La position raffinée est `E_pic=-b/(2a)`. Les pics dont la parabole n’est pas concave ou dont le sommet sort de l’intervalle des trois points sont ignorés.


Le tableau final est un tableau Pandas. Le nombre de chiffres après la virgule est réglé par :

```python
DECIMALS = 8
```

La colonne `position_pic` contient la position du sommet de la parabole. Les colonnes `E_grille`, `intensite_pic`, `a`, `b` et `c` sont également affichées avec ce nombre de décimales.

## Calcul avec un fichier situé dans un autre répertoire

Le chemin de `auto_cor.txt` ne doit pas être défini en dur dans le notebook. Utiliser le script dédié :

```bash
./calcul_spectre.sh /chemin/vers/auto_cor.txt /chemin/vers/file_spectrum.txt
```

Le deuxième argument est facultatif. S’il est omis, le spectre est écrit à côté du fichier d’autocorrélation :

```bash
./calcul_spectre.sh /donnees/propagation/auto_cor.txt
```

Le script crée automatiquement `spectrum_run.json` à côté du spectre. Le notebook lit ce manifeste et retrouve les deux chemins sans qu’il soit nécessaire de modifier son code. Il suffit ensuite d’ouvrir le notebook depuis le dossier contenant `spectrum_run.json`.

Si le notebook est utilisé sans ce script, `DATA_FILE` vaut `None` et seul le spectre est analysé ; le chemin de l’autocorrélation peut alors être renseigné manuellement dans la cellule de configuration.

## Que signifie dE ?

`dE` est le **pas de discrétisation en énergie**. Si les énergies sont calculées comme

```text
E = Emin, Emin+dE, Emin+2*dE, ...
```

alors une petite valeur de `dE` donne plus de points et une meilleure résolution, mais demande davantage de calcul. `dE` doit être positif et utiliser les mêmes unités que l’énergie `E`.

Dans l’ancien fichier, `dE = -1.0` était une valeur spéciale : elle ne représentait pas un vrai pas d’énergie ; elle demandait au programme de calculer automatiquement un pas. Pour un débutant, le fichier `param.nml` conserve maintenant la convention `dE = -1.0`. Le programme calcule alors automatiquement une valeur positive et l’écrit dans `spectrum_parameters.txt`.

De même, `Emax = -1.0` demande encore une détermination automatique de la limite maximale autorisée par la fréquence de Nyquist. Il est possible de remplacer cette valeur par une énergie maximale explicite.

## Scripts simplifiés

1. Compiler une seule fois :

```bash
./script.sh
```

2. Calculer le spectre en indiquant le fichier d’entrée, le fichier de sortie et le pas `dE` :

```bash
./calcul_spectre.sh /donnees/auto_cor.txt /resultats/file_spectrum.txt 1.0e-5
```

Le troisième argument `dE` est facultatif ; s’il est absent, le script utilise `-1.0`, donc le programme calcule automatiquement le pas. Le script crée aussi `spectrum_run.json` et `spectrum_parameters.txt`.


## Rapport des paramètres utilisés

Après chaque calcul, le programme écrit `spectrum_parameters.txt`. Ce fichier indique les valeurs réellement utilisées, même lorsque `dE=-1` dans le namelist :

```text
dE              = ...
Emin            = ...
Emax            = ...
```

Il contient aussi `dt`, `t0`, `tmax`, `npts` et `nptE`.

## Un seul script : spectre.sh

Il n’y a maintenant qu’un seul script. Depuis le dossier qui contient `auto_cor.txt` :

```bash
./spectre.sh
```

Le programme comprend alors automatiquement :

```text
auto_cor.txt       = fichier d’entrée dans le dossier courant
file_spectrum.txt  = fichier de sortie dans le dossier courant
dE                 = -1, calcul automatique
```

On peut aussi donner les paramètres :

```bash
./spectre.sh /chemin/auto_cor.txt /chemin/file_spectrum.txt 1.0e-5
```

Le script compile automatiquement `spectre.f90`, lance le calcul et écrit également `spectrum_parameters.txt` et `spectrum_run.json` près du fichier spectral.

## Paramètres souples de spectre.sh

La commande simple reste :

```bash
./spectre.sh
```

Elle utilise `./auto_cor.txt`, `dE=-1`, `Emin=0`, `Emax=-1` et `n_filter=1`.

Tous les paramètres peuvent être donnés dans cet ordre :

```bash
./spectre.sh fichier_auto_cor fichier_spectrum dE Emin Emax n_filter
```

Exemple :

```bash
./spectre.sh auto_cor.txt file_spectrum.txt 1.0e-4 0.0 0.5 3
```

- `dE=-1` : calcul automatique du pas en énergie ;
- `Emin` : énergie minimale ;
- `Emax=-1` : calcul automatique de la limite maximale ;
- `n_filter` : ordre du filtre exponentiel utilisé pour la colonne `Re[S_window(E)]`.

Le filtre est calculé dans le programme Fortran par :

```text
C_window(t) = C(t) exp(-((t-t0)/(tmax-t0))^n_filter)
```

La valeur `n_filter` est donc utilisée pour produire la quatrième colonne du spectre, `Re[S_window(E)]`. Sa valeur effectivement utilisée est maintenant aussi écrite dans `spectrum_parameters.txt`.

## Tableau des pics

Le notebook n’affiche plus la liste brute dans la sortie texte. Il affiche un tableau Pandas avec les colonnes `position_pic`, `intensite_pic`, `a`, `b` et `c`, arrondies selon `DECIMALS`.
