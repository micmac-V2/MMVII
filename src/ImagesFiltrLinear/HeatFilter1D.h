// GENERATE with CLAUDE see recapitulatif_filtre_chaleur_1D.md  for a summary of prompt

// HeatFilter1D.h
//
// Filtrage 1D par équation de la chaleur à conductivité variable.
//
//   du/dt = d/dx ( k(x) du/dx )            (forme conservative)
//
// Discrétisation :
//   - grille uniforme de pas h, n échantillons (cellules centrées en x_i = i*h)
//   - k fourni aux noeuds, ramené aux interfaces k_{i+1/2} (moyenne arithmétique
//     ou harmonique)
//   - conditions de Neumann (flux nul) : k_{-1/2} = k_{n-1/2} = 0
//   - Euler implicite :  (I - dt*L) u^{m+1} = u^m
//     La matrice est tridiagonale, symétrique, à diagonale strictement dominante
//     (donc factorisation de Thomas sans pivot, stable) et constante au cours des
//     pas : elle est factorisée une seule fois dans le constructeur.
//
// Propriétés (garanties par la construction) :
//   - inconditionnellement stable (aucune contrainte sur dt)
//   - conservation de la somme des échantillons (somme des lignes de la matrice = 1)
//   - principe du maximum : min(u) <= u_out <= max(u)
//
// Dépendances : bibliothèque standard uniquement.

#ifndef HEAT_FILTER_1D_H
#define HEAT_FILTER_1D_H

#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <vector>

// namespace Claude_HeatFilter1D

template <typename Type>
class cHeatFilter1D
{
    static_assert(std::is_floating_point<Type>::value,
                  "HeatFilter1D : Type doit etre un type flottant");

public:
    // Moyenne utilisée pour passer de k aux noeuds à k aux interfaces.
    enum class InterfaceMean
    {
        Arithmetic,  // (k_i + k_{i+1}) / 2
        Harmonic     // 2 k_i k_{i+1} / (k_i + k_{i+1}) : une zone a k~0 bloque la diffusion
    };

    struct Params
    {
        Type          step      = Type(1);   // pas d'espace h (> 0)
        Type          dt        = Type(1);   // pas de temps (>= 0)
        std::size_t   nbSteps   = 1;         // nombre de pas ; temps de filtrage = dt * nbSteps
        InterfaceMean mean      = InterfaceMean::Harmonic;
    };

    // conductivity : tableau de taille size, valeurs >= 0. Copié (peut être libéré après).
    cHeatFilter1D(const Type* conductivity, std::size_t size, const Params& params)
        : mSize(size), mParams(params)
    {
        if (size > 0 && conductivity == nullptr)
            throw std::invalid_argument("HeatFilter1D : conductivite nulle");
        if (!(params.step > Type(0)))
            throw std::invalid_argument("HeatFilter1D : le pas d'espace doit etre > 0");
        if (params.dt < Type(0))
            throw std::invalid_argument("HeatFilter1D : dt doit etre >= 0");
        for (std::size_t i = 0; i < size; ++i)
            if (!(conductivity[i] >= Type(0)))
                throw std::invalid_argument("HeatFilter1D : conductivite negative ou NaN");

        Factorize(conductivity);
    }

    // Filtre `in` vers `out` (n = size du constructeur). in == out est autorisé.
    void Apply(const Type* in, Type* out) const
    {
        if (mSize == 0)
            return;
        if (in == nullptr || out == nullptr)
            throw std::invalid_argument("HeatFilter1D::Apply : pointeur nul");

        if (in != out)
            for (std::size_t i = 0; i < mSize; ++i)
                out[i] = in[i];

        for (std::size_t step = 0; step < mParams.nbSteps; ++step)
            SolveInPlace(out);
    }

    std::size_t Size() const { return mSize; }

private:
    // Construit les interfaces puis factorise la matrice tridiagonale A = I - r*K,
    // avec r = dt/h^2 :
    //   A_ii       = 1 + r (a_{i-1} + a_i)      (a_{-1} = a_{n-1} = 0 : Neumann)
    //   A_i,i+1    = A_i+1,i = - r a_i
    void Factorize(const Type* k)
    {
        if (mSize == 0)
            return;

        const Type r = mParams.dt / (mParams.step * mParams.step);

        // a[i] = conductivite a l'interface entre i et i+1 (i = 0..n-2)
        std::vector<Type> a(mSize > 1 ? mSize - 1 : 0);
        for (std::size_t i = 0; i + 1 < mSize; ++i)
            a[i] = InterfaceValue(k[i], k[i + 1]);

        mLower.assign(mSize, Type(0));
        mInvPivot.assign(mSize, Type(0));
        mUpperScaled.assign(mSize, Type(0));

        // Factorisation de Thomas : on stocke l_i, 1/pivot_i, et c'_i = u_i / pivot_i
        for (std::size_t i = 0; i < mSize; ++i)
        {
            const Type aLeft  = (i > 0)            ? a[i - 1] : Type(0);
            const Type aRight = (i + 1 < mSize)    ? a[i]     : Type(0);

            const Type diag  = Type(1) + r * (aLeft + aRight);
            const Type lower = -r * aLeft;      // A_{i,i-1}
            const Type upper = -r * aRight;     // A_{i,i+1}

            const Type pivot = (i == 0) ? diag : diag - lower * mUpperScaled[i - 1];
            mLower[i]       = lower;
            mInvPivot[i]    = Type(1) / pivot;
            mUpperScaled[i] = upper * mInvPivot[i];
        }
    }

    Type InterfaceValue(Type k1, Type k2) const
    {
        if (mParams.mean == InterfaceMean::Arithmetic)
            return Type(0.5) * (k1 + k2);

        const Type s = k1 + k2;
        return (s > Type(0)) ? (Type(2) * k1 * k2) / s : Type(0);
    }

    // Résout A x = b ; b est dans x en entrée, la solution est dans x en sortie.
    void SolveInPlace(Type* x) const
    {
        // Descente
        x[0] *= mInvPivot[0];
        for (std::size_t i = 1; i < mSize; ++i)
            x[i] = (x[i] - mLower[i] * x[i - 1]) * mInvPivot[i];

        // Remontée
        for (std::size_t i = mSize - 1; i-- > 0;)
            x[i] -= mUpperScaled[i] * x[i + 1];
    }

    std::size_t       mSize;
    Params            mParams;
    std::vector<Type> mLower;        // sous-diagonale  A_{i,i-1}
    std::vector<Type> mInvPivot;     // 1 / pivot_i
    std::vector<Type> mUpperScaled;  // c'_i = A_{i,i+1} / pivot_i
};

// Raccourci pour un usage ponctuel (refactorise a chaque appel).
template <typename Type>
void HeatFilter1DApply(const Type* in, Type* out, const Type* conductivity, std::size_t size,
                       const typename cHeatFilter1D<Type>::Params& params)
{
    cHeatFilter1D<Type>(conductivity, size, params).Apply(in, out);
}

#endif // HEAT_FILTER_1D_H
