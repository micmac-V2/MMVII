# Filtrage 1D par équation de la chaleur à conductivité variable

Récapitulatif de la conversation du 29 octobre 2026.

Fichiers produits :

- `HeatFilter1D.h` : classe template header-only (bibliothèque standard uniquement)
- `test_HeatFilter1D.cpp` : programme de test (`float` et `double`)

Compilation du test :

```bash
g++ -std=c++17 -O2 -Wall -Wextra test_HeatFilter1D.cpp -o test_HeatFilter1D
```

---

## 1. Spécifications retenues

| Point | Choix |
|---|---|
| Dimension | 1D |
| Conductivité k(x) | fournie par l'utilisateur (tableau de même taille que le signal) |
| Forme de l'équation | conservative : ∂u/∂t = ∂/∂x ( k(x) ∂u/∂x ) |
| Schéma temporel | Euler implicite |
| Conditions aux limites | Neumann (flux nul) |
| Grille | uniforme, pas h |
| Type | template sur le type flottant |
| Interface | classe ; entrée `const Type*`, sortie `Type*`, taille `size_t` |
| Dépendances | bibliothèque standard seule (Eigen autorisé, mais non nécessaire) |

---

## 2. Forme conservative ou non conservative

**Conservative** : ∂u/∂t = ∂x ( k ∂x u )

**Non conservative** : ∂u/∂t = k ∂²x u

Elles diffèrent d'un terme d'advection :

∂x(k ∂x u) = k ∂²x u + k'(x) ∂x u

Avantages de la forme conservative :

1. **Conservation** : avec Neumann, la somme des échantillons est exactement préservée.
2. **Sens physique** : loi de Fick/Fourier, le flux vaut −k ∂x u. Une zone à k ≈ 0 agit comme une barrière.
3. **Principe du maximum** : pour k ≥ 0, pas de nouveaux extrema ni d'oscillations parasites.
4. **Matrice symétrique** semi-définie négative, donc schéma stable et solveur efficace.

---

## 3. Discrétisation

Grille uniforme (pas h), k ramené aux interfaces k_{i+1/2} :

- moyenne **arithmétique** : (k_i + k_{i+1}) / 2
- moyenne **harmonique** (défaut) : 2 k_i k_{i+1} / (k_i + k_{i+1}), plus robuste aux forts contrastes

Opérateur :

(Lu)_i = [ k_{i+1/2} (u_{i+1} − u_i) − k_{i−1/2} (u_i − u_{i−1}) ] / h²

Neumann : on pose k_{−1/2} = k_{n−1/2} = 0 (flux nul aux extrémités, sans cellule fantôme).

Euler implicite, avec r = dt / h² :

(I − dt·L) u^{m+1} = u^m

La matrice est tridiagonale, symétrique et à diagonale strictement dominante :

- diagonale : 1 + r (a_{i−1} + a_i), avec a_i = k_{i+1/2}
- hors-diagonale : −r a_i

Elle est **constante** au cours des pas : factorisation de Thomas (sans pivot) faite **une seule fois** dans le constructeur, puis chaque pas coûte une descente et une remontée en O(n), sans allocation.

---

## 4. Utilisation

```cpp
HeatFilter1D<double>::Params p;
p.step    = 1.0;   // pas d'espace h
p.dt      = 5.0;   // pas de temps
p.nbSteps = 20;    // temps de filtrage = dt * nbSteps
p.mean    = HeatFilter1D<double>::InterfaceMean::Harmonic;

HeatFilter1D<double> filter(k, n, p);   // k : const double*, copié
filter.Apply(in, out);                  // in == out autorisé
```

Raccourci ponctuel : `HeatFilter1DApply(in, out, k, n, p)`.

Le constructeur lève `std::invalid_argument` si le pas est ≤ 0, si dt < 0, ou si k est négatif ou NaN.

---

## 5. Propriétés vérifiées par les tests

- conservation de la somme du signal
- principe du maximum (pas d'oscillations parasites)
- signal constant invariant
- résultat identique pour `in == out` et `in != out`
- barrière k = 0 : aucune fuite d'un côté à l'autre
- cas limite n = 1

---

## 6. Unité de la conductivité

Dans ce code, k est une **diffusivité**, de dimension (longueur)² / temps, c'est-à-dire de l'unité de h² / dt. Seul le produit

r = k · dt / h²

intervient dans le calcul.

---

## 7. Équivalence avec un filtre gaussien (k constant)

Pour k constant, l'équation de la chaleur revient à convoluer par une gaussienne de variance σ² = 2 k T, avec T = N · dt.

Dans le code, la variance totale est :

**σ² = 2 · N · k · dt / h²** (σ en échantillons)

Elle ne dépend pas du découpage en temps : chaque pas implicite apporte une variance 2r, et les variances s'additionnent.

Pour viser un σ donné avec N itérations (h = 1) :

**k · dt = σ² / (2N)**

Exemple : σ = 3 échantillons, N = 10 → k · dt = 0,45 (par exemple k = 1 et dt = 0,45).
Si σ est en unités physiques, on utilise σ / h.

### Qualité de l'approximation

- La variance est exacte, mais le noyau n'est que **approximativement** gaussien.
- Un seul pas implicite donne un noyau à queues exponentielles (type Laplace, excès de kurtosis égal à 3).
- Après N pas, l'excès de kurtosis est d'environ 3 / N (0 pour une gaussienne).
- En pratique, N entre 4 et 10 donne une très bonne approximation. Au-delà, le coût augmente (linéaire en N) sans gain visible.

### Limites

- Loin des bords (plus de 3σ environ) ; près des bords, Neumann se comporte comme une réflexion du signal.
- Avec k(x) variable : σ(x) ≈ √(2 k(x) T) localement, mais ce n'est plus un vrai filtre gaussien (le terme d'advection k' u' déforme légèrement le noyau).
