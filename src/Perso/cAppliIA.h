/**
 * @file convertisseur.cpp
 * @brief Convertisseur de devises interactif avec complétion automatique.
 *
 * Le fichier regroupe :
 *  - la classe Devise : un code de devise et son taux par rapport à l'euro ;
 *  - la classe BaseDevises : le catalogue des devises, construit soit avec des
 *    taux « en dur », soit à partir d'un fichier texte ;
 *  - la classe Convertisseur : conversion d'un montant et saisie interactive,
 *    avec complétion par la touche Tab si GNU Readline est disponible ;
 *  - la fonction CONVERTIR() : point d'entrée destiné à être appelé par un
 *    programme hôte (ce fichier ne définit pas de @c main).
 *
 * @par Format du fichier de taux
 * Une ligne par devise : <tt>CODE TAUX</tt>, séparés par des espaces ou une
 * tabulation. Les taux sont exprimés pour 1 EUR. Les lignes vides et celles
 * commençant par @c # sont ignorées ; la virgule décimale est acceptée.
 * @code
 * # Taux pour 1 EUR
 * EUR 1.0
 * USD 1.08
 * JPY 163
 * @endcode
 *
 * @par Compilation
 * C++17 est requis. GNU Readline est @b optionnel : il n'est utilisé que si la
 * macro @c HAVE_READLINE est définie à la compilation. Sans elle, la saisie se
 * fait avec @c std::getline (pas de complétion) et aucune dépendance externe
 * n'est nécessaire.
 * @code
 * # Sans Readline (par défaut)
 * g++ -std=c++17 -o monprogramme main.cpp convertisseur.cpp
 *
 * # Avec Readline (complétion par Tab)
 * g++ -std=c++17 -DHAVE_READLINE -o monprogramme main.cpp convertisseur.cpp -lreadline
 * @endcode
 *
 * @warning GNU Readline est distribué sous licence GPL. Distribuer un binaire
 *          lié à Readline impose de respecter les conditions de la GPL pour ce
 *          binaire. D'où le choix de ne l'activer que sur demande.
 *
 * @section regles Règles de conversion
 * Les règles ci-dessous sont repérées dans le code par leur numéro (@c [R7],
 * etc.). Celles marquées « limite connue » décrivent un comportement actuel
 * jugé perfectible.
 *
 * @par Calcul de la conversion
 *  - R1 : un taux est le nombre d'unités de la devise pour 1 EUR.
 *  - R2 : la conversion passe toujours par l'euro :
 *    <tt>résultat = montant / taux(départ) * taux(arrivée)</tt>.
 *  - R3 : l'EUR vaut 1 ; convertir depuis ou vers l'euro n'est qu'une seule
 *    opération.
 *  - R4 : convertir une devise vers elle-même redonne le montant (à l'arrondi
 *    flottant près).
 *  - R5 : un code absent de la base lève @c std::out_of_range.
 *  - R6 : aucun frais, commission ni écart achat/vente : un seul taux par
 *    devise.
 *
 * @par Validité d'une devise
 *  - R7 : le code est mis en majuscules à la construction.
 *  - R8 : un code vide est refusé.
 *  - R9 : un taux négatif ou nul est refusé (@c std::invalid_argument).
 *  - R10 : aucun contrôle de format sur le code (« 3 lettres » non vérifié) ;
 *    limite connue.
 *
 * @par Constitution de la base
 *  - R11 : constructeur par défaut : 7 devises en dur (EUR, USD, GBP, CHF,
 *    JPY, CAD, MAD).
 *  - R12 : constructeur depuis un fichier : une ligne égale un couple
 *    <tt>CODE TAUX</tt>, séparé par des espaces ou une tabulation.
 *  - R13 : blancs aux extrémités retirés ; lignes vides et lignes commençant
 *    par @c # ignorées.
 *  - R14 : la virgule décimale est acceptée dans le fichier.
 *  - R15 : une ligne doit contenir exactement deux champs (un code, un
 *    nombre), sinon erreur avec le nom du fichier et le numéro de ligne.
 *  - R16 : un code en double n'est pas signalé, la dernière ligne l'emporte ;
 *    limite connue.
 *  - R17 : fichier introuvable ou ne contenant aucune devise : erreur.
 *
 * @par Saisie du montant
 *  - R18 : les blancs autour de la saisie sont retirés.
 *  - R19 : toute la saisie doit être un nombre (@c 12abc est refusé).
 *  - R20 : un montant négatif est refusé ; zéro est accepté.
 *  - R21 : la virgule décimale n'est pas acceptée à la saisie ; limite connue.
 *  - R22 : @c nan et @c inf ne sont pas filtrés ; limite connue.
 *
 * @par Saisie des devises
 *  - R23 : codes nettoyés puis mis en majuscules : la casse est ignorée.
 *  - R24 : les deux codes doivent exister, sinon « Devise inconnue. » sans
 *    préciser lequel ; limite connue.
 *  - R25 : avec Readline, Tab complète d'après le préfixe, casse ignorée ;
 *    sans Readline, aucune complétion.
 *
 * @par Affichage et codes de retour
 *  - R26 : résultat affiché avec 2 décimales pour toutes les devises ; limite
 *    connue (inexact pour le JPY ou le KWD).
 *  - R27 : ordre des étapes : montant, devise de départ, devise d'arrivée ;
 *    une erreur arrête la session.
 *  - R28 : CONVERTIR() retourne 0 en cas de succès, 1 sinon (saisie invalide,
 *    devise inconnue, fin de saisie, erreur de fichier) ; les exceptions sont
 *    rattrapées et affichées sur @c std::cerr.
 *
 * @note Les classes sont dans un espace de noms anonyme, donc privées à ce
 *       fichier. Pour que Doxygen les documente, activer l'option
 *       @c EXTRACT_ANON_NSPACES (et @c EXTRACT_ALL ou @c EXTRACT_STATIC).
 */

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

//#define HAVE_READLINE  0

#ifdef HAVE_READLINE
#include <readline/history.h>
#include <readline/readline.h>
#endif

namespace {

// ============================================================
//  Utilitaires de chaîne
// ============================================================

/**
 * @brief Met une chaîne en majuscules.
 *
 * @note Règles R7 et R23.
 * @param s Chaîne à convertir (passée par copie).
 * @return La chaîne en majuscules.
 */
std::string EnMajuscules(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return std::toupper(c); });
    return s;
}

/**
 * @brief Retire les blancs (espace, tabulation, retours chariot) au début
 *        et à la fin d'une chaîne.
 * @param s Chaîne à nettoyer.
 * @return La chaîne sans blancs aux extrémités, ou une chaîne vide si @p s
 *         ne contient que des blancs.
 *
 * @note Règles R13, R18 et R23.
 */
std::string Nettoyer(const std::string& s) {
    const char* blancs = " \t\r\n";
    size_t debut = s.find_first_not_of(blancs);
    if (debut == std::string::npos) return "";
    size_t fin = s.find_last_not_of(blancs);
    return s.substr(debut, fin - debut + 1);
}

// ============================================================
//  Classe Devise
// ============================================================

/**
 * @class Devise
 * @brief Une devise : son code et son taux pour 1 EUR.
 *
 * L'euro sert de devise pivot : le taux est le nombre d'unités de la devise
 * obtenues pour 1 EUR (par exemple 1.08 pour l'USD). Le code est
 * automatiquement mis en majuscules.
 */
class Devise {
public:
    /**
     * @brief Construit une devise.
     * @param code Code de la devise (ex. "usd" ou "USD"), mis en majuscules.
     * @param tauxPourUnEuro Nombre d'unités de la devise pour 1 EUR.
     * @throw std::invalid_argument Si le code est vide ou si le taux est
     *        négatif ou nul.
     */
    Devise(std::string code, double tauxPourUnEuro)
        // [R7] code mis en majuscules
        : mCode(EnMajuscules(std::move(code))), mTaux(tauxPourUnEuro)
    {
        // [R10] aucun contrôle de format sur le code (limite connue)
        if (mCode.empty())  // [R8]
            throw std::invalid_argument("Code de devise vide");
        if (mTaux <= 0.0)   // [R9]
            throw std::invalid_argument("Taux invalide pour " + mCode);
    }

    /** @brief Code de la devise, en majuscules (ex. "USD"). */
    const std::string& Code() const { return mCode; }

    /** @brief Nombre d'unités de cette devise pour 1 EUR. */
    double Taux() const { return mTaux; }

    /**
     * @brief Convertit un montant exprimé dans cette devise en euros.
     * @param montant Montant dans cette devise.
     * @return Le montant équivalent en euros.
     */
    double VersEuro(double montant) const { return montant / mTaux; }  // [R1]

    /**
     * @brief Convertit un montant en euros vers cette devise.
     * @param montantEuro Montant en euros.
     * @return Le montant équivalent dans cette devise.
     */
    double DepuisEuro(double montantEuro) const { return montantEuro * mTaux; }  // [R1]

private:
    std::string mCode;  ///< Code de la devise, en majuscules.
    double      mTaux;  ///< Nombre d'unités de cette devise pour 1 EUR.
};

// ============================================================
//  Classe BaseDevises
// ============================================================

/**
 * @class BaseDevises
 * @brief Catalogue de devises, indexées par leur code.
 *
 * Peut être construit avec des taux « en dur » (constructeur par défaut) ou
 * à partir d'un fichier texte (voir BaseDevises(const std::string&)).
 */
class BaseDevises {
public:
    /**
     * @brief Constructeur par défaut : taux « en dur ».
     *
     * Remplit la base avec EUR, USD, GBP, CHF, JPY, CAD et MAD. Les taux
     * sont des valeurs d'exemple, à mettre à jour.
     *
     * @note Règle R11.
     */
    BaseDevises() {
        Ajouter(Devise("EUR", 1.0));
        Ajouter(Devise("USD", 1.08));
        Ajouter(Devise("GBP", 0.85));
        Ajouter(Devise("CHF", 0.96));
        Ajouter(Devise("JPY", 163.0));
        Ajouter(Devise("CAD", 1.47));
        Ajouter(Devise("MAD", 10.8));
    }

    /**
     * @brief Construit la base à partir d'un fichier texte.
     *
     * Chaque ligne contient un couple <tt>CODE TAUX</tt>. Les lignes vides et
     * les commentaires (commençant par @c #) sont ignorés, la virgule
     * décimale est acceptée. Si un code apparaît plusieurs fois, la dernière
     * ligne l'emporte.
     *
     * @note Règles R12 à R17.
     *
     * @param nomFichier Chemin du fichier (absolu, ou relatif au répertoire
     *        courant du programme).
     * @throw std::runtime_error Si le fichier est illisible, si une ligne est
     *        mal formée ou contient un taux invalide (le message indique le
     *        fichier et le numéro de ligne), ou si aucune devise n'est lue.
     */
    explicit BaseDevises(const std::string& nomFichier) {
        std::ifstream f(nomFichier);
        if (!f)  // [R17] fichier introuvable
            throw std::runtime_error("Impossible d'ouvrir " + nomFichier);

        std::string ligne;
        int numLigne = 0;
        while (std::getline(f, ligne)) {
            ++numLigne;
            ligne = Nettoyer(ligne);  // [R13] blancs retirés
            if (ligne.empty() || ligne[0] == '#') continue;  // [R13] vide ou commentaire

            std::replace(ligne.begin(), ligne.end(), ',', '.');  // [R14] virgule décimale

            std::istringstream iss(ligne);
            std::string code, reste;
            double taux;
            // [R15] Exactement deux champs : un code, puis un nombre.
            if (!(iss >> code >> taux) || (iss >> reste))
                throw std::runtime_error(nomFichier + ", ligne " +
                                         std::to_string(numLigne) +
                                         " : format attendu \"CODE TAUX\"");
            try {
                Ajouter(Devise(code, taux));  // [R9] taux invalide : erreur avec n° de ligne
            } catch (const std::invalid_argument& e) {
                throw std::runtime_error(nomFichier + ", ligne " +
                                         std::to_string(numLigne) + " : " + e.what());
            }
        }
        if (mDevises.empty())  // [R17] aucune devise lue
            throw std::runtime_error(nomFichier + " : aucune devise lue");
    }

    /**
     * @brief Ajoute une devise à la base, ou remplace celle de même code.
     * @param d Devise à ajouter.
     */
    void Ajouter(const Devise& d) { mDevises.insert_or_assign(d.Code(), d); }  // [R16] doublon écrasé sans avertissement

    /**
     * @brief Indique si un code de devise est présent dans la base.
     * @param code Code recherché, en majuscules.
     * @return @c true si la devise existe.
     */
    bool Existe(const std::string& code) const {
        return mDevises.count(code) != 0;
    }

    /**
     * @brief Retourne la devise correspondant à un code.
     * @param code Code recherché, en majuscules.
     * @return Référence vers la devise.
     * @throw std::out_of_range Si le code n'existe pas dans la base.
     */
    const Devise& Get(const std::string& code) const {
        auto it = mDevises.find(code);
        if (it == mDevises.end())
            throw std::out_of_range("Devise inconnue : " + code);
        return it->second;
    }

    /**
     * @brief Liste des codes de toutes les devises de la base.
     * @return Les codes, triés par ordre alphabétique.
     */
    std::vector<std::string> Codes() const {
        std::vector<std::string> codes;
        for (const auto& [code, devise] : mDevises) codes.push_back(code);
        return codes;
    }

private:
    std::map<std::string, Devise> mDevises;  ///< Devises indexées par leur code.
};

// ============================================================
//  Classe Convertisseur
// ============================================================

/**
 * @class Convertisseur
 * @brief Convertit des montants entre devises et gère la saisie interactive.
 *
 * La conversion passe par l'euro (devise pivot). Si la macro @c HAVE_READLINE
 * est définie, la saisie interactive utilise GNU Readline : la touche Tab
 * complète le nom des devises, en ignorant la casse. Sinon, la saisie se fait
 * avec @c std::getline, sans complétion.
 */
class Convertisseur {
public:
    /**
     * @brief Construit un convertisseur à partir d'une base de devises.
     * @param base Base utilisée (copiée). Par défaut, la base « en dur »
     *        construite par BaseDevises().
     */
    explicit Convertisseur(BaseDevises base = BaseDevises())
        : mBase(std::move(base)) {}

    /** @brief Base de devises utilisée par ce convertisseur. */
    const BaseDevises& Base() const { return mBase; }

    /**
     * @brief Convertit un montant d'une devise à une autre.
     *
     * Le montant est d'abord converti en euros, puis dans la devise cible.
     *
     * @param montant Montant à convertir.
     * @param de Code de la devise de départ, en majuscules.
     * @param vers Code de la devise d'arrivée, en majuscules.
     * @return Le montant converti.
     * @throw std::out_of_range Si l'un des codes est inconnu.
     *
     * @note Règles R2 à R6.
     */
    double Convertir(double montant, const std::string& de,
                     const std::string& vers) const {
        // [R2] via l'euro ; [R3] l'EUR vaut 1 ; [R4] de == vers redonne le montant ;
        // [R5] Get() lève out_of_range si un code est inconnu ; [R6] aucun frais.
        return mBase.Get(vers).DepuisEuro(mBase.Get(de).VersEuro(montant));
    }

    /**
     * @brief Lance une session interactive de conversion.
     *
     * Affiche les devises disponibles, demande un montant, la devise de
     * départ puis la devise d'arrivée (avec complétion par Tab si Readline
     * est activé), et affiche
     * le résultat avec deux décimales. Les erreurs de saisie sont signalées
     * sur @c std::cerr.
     *
     * @return 0 en cas de succès, 1 en cas de saisie invalide, de devise
     *         inconnue ou de fin de saisie (Ctrl+D).
     *
     * @note Règles R18 à R24 et R26 à R28.
     */
    int Executer() {
        std::string saisie, de, vers;

        std::cout << "Devises disponibles : ";
        for (const auto& code : mBase.Codes()) std::cout << code << " ";
        std::cout << "\n";
#ifdef HAVE_READLINE
        std::cout << "(Tab pour compléter le nom d'une devise)\n";
#endif
        std::cout << "\n";

        // [R27] ordre : montant, devise de départ, devise d'arrivée ;
        // une erreur (ou Ctrl+D) arrête la session avec le code 1 [R28].
        if (!LireLigne("Montant à convertir : ", saisie)) return 1;

        double montant = 0.0;
        try {
            size_t pos = 0;
            montant = std::stod(saisie, &pos);
            // [R19] pos != taille : du texte suit le nombre (ex. "12abc")
            // [R20] négatif refusé, zéro accepté
            // [R21] pas de virgule décimale ici (limite connue)
            // [R22] nan et inf ne sont pas filtrés (limite connue)
            if (pos != saisie.size() || montant < 0)
                throw std::invalid_argument("");
        } catch (...) {
            std::cerr << "Montant invalide.\n";
            return 1;
        }

        if (!LireLigne("Devise de départ : ", de, true))   return 1;
        if (!LireLigne("Devise d'arrivée : ", vers, true)) return 1;

        // [R23] casse ignorée
        de   = EnMajuscules(de);
        vers = EnMajuscules(vers);

        // [R24] message sans préciser quelle devise est en cause (limite connue)
        if (!mBase.Existe(de) || !mBase.Existe(vers)) {
            std::cerr << "Devise inconnue.\n";
            return 1;
        }

        // [R26] 2 décimales pour toutes les devises (limite connue)
        std::cout << std::fixed << std::setprecision(2)
                  << montant << " " << de << " = "
                  << Convertir(montant, de, vers) << " " << vers << "\n";
        return 0;
    }

private:
    BaseDevises mBase;  ///< Base de devises utilisée pour les conversions.

#ifdef HAVE_READLINE
    // ------------------------------------------------------------------
    //  Complétion Readline : comment ça marche
    //
    //  Quand l'utilisateur appuie sur Tab, Readline appelle la fonction
    //  enregistrée dans rl_attempted_completion_function (ici Completion).
    //  Celle-ci délègue à rl_completion_matches(), qui appelle Generateur
    //  en boucle pour obtenir les candidats un par un, puis construit
    //  elle-même le tableau final et décide de ce qu'il faut afficher :
    //    - 0 candidat  : rien ne se passe
    //    - 1 candidat  : il est inséré directement (suivi d'un espace)
    //    - n candidats : le préfixe commun est inséré, et un deuxième Tab
    //                    affiche la liste des candidats
    // ------------------------------------------------------------------

    /**
     * @brief Instance de Convertisseur en cours de saisie.
     *
     * Readline attend de simples fonctions C : pas de méthode non statique
     * ni de lambda capturante. Comme Generateur() et Completion() n'ont pas
     * accès à @c this, ce pointeur statique est positionné par LireLigne()
     * pendant la durée de la saisie.
     */
    static inline const Convertisseur* sCourant = nullptr;

    /**
     * @brief Générateur de candidats de complétion (callback Readline).
     *
     * Retourne, à chaque appel, la devise suivante dont le code commence par
     * le texte tapé (casse ignorée). Les variables statiques conservent l'état
     * entre deux appels d'une même série, car Readline rappelle la fonction
     * jusqu'à ce qu'elle retourne @c nullptr.
     *
     * @param texte Mot en cours de saisie (ce qui précède le curseur).
     * @param etat 0 au premier appel d'une série (réinitialise le parcours),
     *        non nul aux appels suivants.
     * @return Le candidat suivant, alloué avec @c strdup (Readline le libère
     *         avec @c free), ou @c nullptr quand il n'y en a plus.
     */
    static char* Generateur(const char* texte, int etat) {
        static std::vector<std::string> codes;  // codes à parcourir
        static size_t index;                    // position dans la liste
        static std::string prefixe;             // préfixe recherché, en majuscules

        // Premier appel de la série : on (ré)initialise tout.
        if (etat == 0) {
            codes   = sCourant ? sCourant->mBase.Codes()
                               : std::vector<std::string>{};
            index   = 0;
            prefixe = EnMajuscules(texte);  // [R25] "us" -> "US" : casse ignorée
        }

        // On reprend là où le précédent appel s'était arrêté et on s'arrête
        // au premier code qui commence par le préfixe.
        while (index < codes.size()) {
            const std::string& code = codes[index++];
            if (code.compare(0, prefixe.size(), prefixe) == 0)
                // strdup alloue avec malloc : Readline libérera avec free().
                // Il ne faut donc pas retourner c_str() directement.
                return strdup(code.c_str());
        }

        // Plus aucun candidat : fin de la série.
        return nullptr;
    }

    /**
     * @brief Fonction de complétion appelée par Readline sur Tab.
     *
     * Désactive la complétion de noms de fichiers (comportement par défaut
     * de Readline quand aucun candidat n'est trouvé) puis construit le
     * tableau de candidats à partir de Generateur().
     *
     * @param texte Mot à compléter.
     * @param debut Position de début du mot dans la ligne (inutilisée).
     * @param fin Position de fin du mot dans la ligne (inutilisée).
     * @return Le tableau de candidats construit par Readline, ou @c nullptr
     *         s'il n'y en a aucun.
     *
     * @note Règle R25.
     */
    static char** Completion(const char* texte, int debut, int fin) {
        (void)debut;
        (void)fin;

        // Sans ce drapeau, Readline se rabattrait sur la complétion de
        // noms de fichiers quand on retourne nullptr.
        rl_attempted_completion_over = 1;

        return rl_completion_matches(texte, Generateur);
    }
#endif  // HAVE_READLINE

    /**
     * @brief Affiche une invite et lit une ligne au clavier.
     *
     * Avec @c HAVE_READLINE, la saisie passe par Readline. Les réglages
     * globaux modifiés (fonction de complétion) sont sauvegardés puis
     * restaurés à la sortie, pour ne pas perturber un programme hôte qui
     * utiliserait lui aussi Readline. Sans @c HAVE_READLINE, la saisie se fait
     * avec @c std::getline et @p completer est ignoré.
     *
     * @param invite Texte affiché devant la saisie.
     * @param[out] sortie Reçoit la ligne saisie, sans blancs aux extrémités.
     * @param completer @c true pour activer la complétion des codes de devises
     *        par Tab, @c false sinon.
     * @return @c false si la saisie est terminée (Ctrl+D), @c true sinon.
     */
    bool LireLigne(const char* invite, std::string& sortie,
                   bool completer = false) {
#ifdef HAVE_READLINE
        // Sauvegarde des réglages actuels (programme hôte).
        rl_completion_func_t* ancienne = rl_attempted_completion_function;
        const Convertisseur*  ancien   = sCourant;

        // Active (ou non) notre complétion et désigne l'instance à consulter.
        rl_attempted_completion_function = completer ? Completion : nullptr;
        sCourant = this;

        // Chaîne allouée avec malloc, ou nullptr sur Ctrl+D.
        char* ligne = readline(invite);

        // Restauration, quel que soit le résultat de la saisie.
        rl_attempted_completion_function = ancienne;
        sCourant = ancien;

        if (!ligne) return false;

        sortie = Nettoyer(ligne);
        free(ligne);  // free et non delete : allocation faite par Readline
        return true;
#else
        // Repli sans Readline : saisie simple, sans complétion.
        (void)completer;
        std::cout << invite << std::flush;
        std::string ligne;
        if (!std::getline(std::cin, ligne)) return false;  // fin de saisie
        sortie = Nettoyer(ligne);
        return true;
#endif
    }
};

}  // namespace

// ============================================================
//  Point d'entrée
// ============================================================

/**
 * @brief Lance la conversion de devises interactive.
 *
 * Point d'entrée destiné au programme hôte. Les exceptions sont rattrapées
 * et affichées sur @c std::cerr.
 *
 * @param fichierTaux Chemin d'un fichier de taux (voir BaseDevises). Si la
 *        chaîne est vide (défaut), les taux « en dur » sont utilisés.
 * @return 0 en cas de succès, 1 en cas d'erreur (fichier invalide, saisie
 *         invalide, devise inconnue, fin de saisie).
 *
 * @note Règle R28.
 *
 * @note Si cette fonction est déclarée dans un header, la valeur par défaut
 *       de @p fichierTaux doit figurer dans la déclaration et non dans la
 *       définition.
 *
 * @code
 * CONVERTIR();            // taux en dur
 * CONVERTIR("taux.txt");  // taux lus dans le fichier
 * @endcode
 */
int CONVERTIR(const std::string& fichierTaux = "") {
    try {
        Convertisseur convertisseur = fichierTaux.empty()
            ? Convertisseur()
            : Convertisseur(BaseDevises(fichierTaux));
        return convertisseur.Executer();
    } catch (const std::exception& e) {
        std::cerr << "Erreur : " << e.what() << "\n";
        return 1;
    }
}
